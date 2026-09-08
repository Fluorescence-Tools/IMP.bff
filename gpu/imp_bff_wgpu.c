/*
 * The GPU backend for IMP::bff's diffusion solver: WGSL through wgpu-native,
 * loaded at run time.
 *
 * This library links nothing. `webgpu.h` is vendored beside it, and every
 * wgpu entry point is resolved with dlopen/GetProcAddress from a path handed
 * in by the caller -- so it builds on a machine with no wgpu present, IMP's
 * module build never sees wgpu, and neither conda recipe's dependency list
 * changes. That is the constraint PRD-110 defends and the reason this is a
 * plugin rather than a build option.
 *
 * The contract is IMP/bff/Compute.h, restated below because that header is
 * C++ and this file is C: a backend may be built by a different compiler than
 * the library that loads it, and a std::vector across that line is a promise
 * neither side can keep. `test/test_compute_abi.py` fails if the two drift.
 *
 * The shader is gpu/diffusion.wgsl, carried here as a string by
 * gpu/diffusion_wgsl.h. Read okf/validation/gpu_diffusion_is_worth_it.md
 * before changing it: every choice in it is a measurement, and the one that
 * looks most like a detail -- computing the self weight from the *rounded*
 * neighbour weights -- is worth three orders of magnitude of accuracy.
 */

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "webgpu.h"
#include "diffusion_wgsl.h"

#ifdef _WIN32
#  include <windows.h>
#  define IMPBFF_EXPORT __declspec(dllexport)
#  define DLOPEN(p) ((void*)LoadLibraryA(p))
#  define DLSYM(h, s) ((void*)GetProcAddress((HMODULE)(h), (s)))
#else
#  include <dlfcn.h>
#  define IMPBFF_EXPORT __attribute__((visibility("default")))
#  define DLOPEN(p) dlopen((p), RTLD_NOW | RTLD_LOCAL)
#  define DLSYM(h, s) dlsym((h), (s))
#endif

/* ---- the contract, mirroring IMP/bff/Compute.h -------------------------- */
#define IMPBFF_COMPUTE_BACKEND_ABI 1
typedef int (*ImpBffPropagateFn)(const double* cur, const double* d,
                                 const double* decay, const double* bounds,
                                 int ng, int flux_form, int n_steps, int n_out,
                                 double* out_fluorescence, int n_fluorescence,
                                 double* out_density);
struct ImpBffComputeBackend {
    int abi;
    const char* name;
    ImpBffPropagateFn propagate;
};

/* FLUX_SMOLUCHOWSKI = 0, FLUX_ITO = 1, as DiffusionSolver.h defines them. */
#define FLUX_SMOLUCHOWSKI 0

/* How many partial sums the population reduction uses. One workgroup for the
   whole grid cost 16% of the run; anything from 16 to 256 measures the same. */
#define NPART 64
/* Steps per command buffer. One buffer for four thousand dispatches is a lot
   of recorded state; five hundred keeps the queue busy without it. */
#define BATCH 500
/* Below these the transfer is worth more than the arithmetic, so decline and
   let the CPU have it. */
#define MIN_NG 12
#define MIN_STEPS 32

/* ---- the wgpu entry points we use --------------------------------------- */
#define WGPU_FUNCTIONS(X) X(CreateInstance) X(InstanceRequestAdapter)      \
    X(InstanceProcessEvents) X(AdapterRequestDevice) X(AdapterGetInfo)     \
    X(AdapterInfoFreeMembers) X(AdapterHasFeature) X(DeviceGetQueue)       \
    X(DeviceCreateShaderModule) X(DeviceCreateBuffer)                      \
    X(DeviceCreateBindGroupLayout) X(DeviceCreateBindGroup)                \
    X(DeviceCreatePipelineLayout) X(DeviceCreateComputePipeline)           \
    X(DeviceCreateCommandEncoder) X(CommandEncoderBeginComputePass)        \
    X(CommandEncoderCopyBufferToBuffer) X(CommandEncoderFinish)            \
    X(ComputePassEncoderSetPipeline) X(ComputePassEncoderSetBindGroup)     \
    X(ComputePassEncoderDispatchWorkgroups) X(ComputePassEncoderEnd)       \
    X(QueueSubmit) X(QueueWriteBuffer) X(BufferMapAsync)                   \
    X(BufferGetConstMappedRange) X(BufferUnmap) X(BufferRelease)           \
    X(ShaderModuleRelease) X(BindGroupRelease)                             \
    X(BindGroupLayoutRelease) X(PipelineLayoutRelease)                     \
    X(ComputePipelineRelease) X(CommandEncoderRelease)                     \
    X(CommandBufferRelease) X(ComputePassEncoderRelease)

/* The header spells its function-pointer types WGPUProc<Name>, and the
   entry points wgpu<Name>; the list above carries <Name> so both follow. */

/* wgpuDevicePoll is wgpu-native's own, not webgpu.h's: it is how a program
   with no event loop drives the queue to completion and drains the map
   callbacks. Declared here so that only webgpu.h has to be vendored. */
typedef WGPUBool (*WGPUProcDevicePoll)(WGPUDevice, WGPUBool, uint64_t const*);

struct Wgpu {
#define DECLARE(name) WGPUProc##name name;
    WGPU_FUNCTIONS(DECLARE)
    WGPUProcDevicePoll DevicePoll;
#undef DECLARE
    void* lib;
    WGPUInstance instance;
    WGPUAdapter adapter;
    WGPUDevice device;
    WGPUQueue queue;
    int f16;
    char name[192];
};
static struct Wgpu G;
static struct ImpBffComputeBackend gBackend;

/* ---- async replies, without an event loop ------------------------------- */
struct AdapterReply { WGPUAdapter adapter; int done; };
static void on_adapter(WGPURequestAdapterStatus status, WGPUAdapter adapter,
                       WGPUStringView message, void* u1, void* u2) {
    struct AdapterReply* r = (struct AdapterReply*)u1;
    (void)message; (void)u2;
    r->adapter = (status == WGPURequestAdapterStatus_Success) ? adapter : NULL;
    r->done = 1;
}
struct DeviceReply { WGPUDevice device; int done; };
static void on_device(WGPURequestDeviceStatus status, WGPUDevice device,
                      WGPUStringView message, void* u1, void* u2) {
    struct DeviceReply* r = (struct DeviceReply*)u1;
    (void)message; (void)u2;
    r->device = (status == WGPURequestDeviceStatus_Success) ? device : NULL;
    r->done = 1;
}
struct MapReply { int done; int ok; };
static void on_map(WGPUMapAsyncStatus status, WGPUStringView message,
                   void* u1, void* u2) {
    struct MapReply* r = (struct MapReply*)u1;
    (void)message; (void)u2;
    r->ok = (status == WGPUMapAsyncStatus_Success);
    r->done = 1;
}

static WGPUStringView sv(const char* s) {
    WGPUStringView v;
    v.data = s;
    v.length = s ? strlen(s) : 0;
    return v;
}

/* ---- host-side precomputation ------------------------------------------- */
/* IEEE binary32 -> binary16, round to nearest even. The weights are small and
   well scaled (they are products of a mobility and a decay factor), so the
   subnormal and overflow paths are for safety rather than for this data. */
static uint16_t f16_of(float f) {
    uint32_t x;
    uint32_t sign, exp, man;
    memcpy(&x, &f, 4);
    sign = (x >> 16) & 0x8000u;
    exp = (x >> 23) & 0xFFu;
    man = x & 0x7FFFFFu;
    if (exp == 0xFF) return (uint16_t)(sign | 0x7C00u | (man ? 0x200u : 0u));
    if (exp > 112) {                      /* normal in half */
        int e = (int)exp - 112;
        if (e >= 0x1F) return (uint16_t)(sign | 0x7C00u);   /* overflow */
        {
            uint32_t h = (uint32_t)(e << 10) | (man >> 13);
            uint32_t rest = man & 0x1FFFu;
            if (rest > 0x1000u || (rest == 0x1000u && (h & 1u))) h++;
            return (uint16_t)(sign | h);
        }
    }
    if (exp >= 103) {                     /* subnormal in half */
        uint32_t m = man | 0x800000u;
        int shift = 126 - (int)exp;
        uint32_t h = m >> (shift + 1);
        uint32_t rest = m & ((1u << (shift + 1)) - 1u);
        uint32_t half = 1u << shift;
        if (rest > half || (rest == half && (h & 1u))) h++;
        return (uint16_t)(sign | h);
    }
    return (uint16_t)sign;
}

/* ---- loading ------------------------------------------------------------ */
static int resolve(const char* path, char* err, size_t nerr) {
    void* h = DLOPEN(path);
    if (!h) {
        snprintf(err, nerr, "could not open the wgpu library at '%s'", path);
        return 0;
    }
#define RESOLVE(name)                                                         \
    G.name = (WGPUProc##name)DLSYM(h, "wgpu" #name);                          \
    if (!G.name) {                                                            \
        snprintf(err, nerr, "'%s' has no symbol wgpu%s -- is it wgpu-native?",\
                 path, #name);                                                \
        return 0;                                                             \
    }
    WGPU_FUNCTIONS(RESOLVE)
#undef RESOLVE
    G.DevicePoll = (WGPUProcDevicePoll)DLSYM(h, "wgpuDevicePoll");
    if (!G.DevicePoll) {
        snprintf(err, nerr, "'%s' has no wgpuDevicePoll -- is it wgpu-native?",
                 path);
        return 0;
    }
    G.lib = h;
    return 1;
}

static int open_device(char* err, size_t nerr) {
    struct AdapterReply ar = {NULL, 0};
    struct DeviceReply dr = {NULL, 0};
    WGPURequestAdapterOptions opts;
    WGPURequestAdapterCallbackInfo aci;
    WGPURequestDeviceCallbackInfo dci;
    WGPUDeviceDescriptor dd;
    WGPUAdapterInfo info;
    WGPUFeatureName want = WGPUFeatureName_ShaderF16;
    int guard;

    G.instance = G.CreateInstance(NULL);
    if (!G.instance) { snprintf(err, nerr, "wgpuCreateInstance returned null"); return 0; }

    memset(&opts, 0, sizeof opts);
    opts.powerPreference = WGPUPowerPreference_HighPerformance;
    memset(&aci, 0, sizeof aci);
    aci.mode = WGPUCallbackMode_AllowProcessEvents;
    aci.callback = on_adapter;
    aci.userdata1 = &ar;
    G.InstanceRequestAdapter(G.instance, &opts, aci);
    for (guard = 0; !ar.done && guard < 100000; ++guard)
        G.InstanceProcessEvents(G.instance);
    if (!ar.adapter) { snprintf(err, nerr, "no GPU adapter"); return 0; }
    G.adapter = ar.adapter;

    G.f16 = G.AdapterHasFeature(G.adapter, want) ? 1 : 0;

    memset(&dd, 0, sizeof dd);
    if (G.f16) {
        dd.requiredFeatureCount = 1;
        dd.requiredFeatures = &want;
    }
    memset(&dci, 0, sizeof dci);
    dci.mode = WGPUCallbackMode_AllowProcessEvents;
    dci.callback = on_device;
    dci.userdata1 = &dr;
    G.AdapterRequestDevice(G.adapter, &dd, dci);
    for (guard = 0; !dr.done && guard < 100000; ++guard)
        G.InstanceProcessEvents(G.instance);
    if (!dr.device) { snprintf(err, nerr, "no GPU device"); return 0; }
    G.device = dr.device;
    G.queue = G.DeviceGetQueue(G.device);
    if (!G.queue) { snprintf(err, nerr, "no queue"); return 0; }

    /* The name is what get_compute_backend_name() reports, so that nobody
       measures a CPU and believes they measured a GPU. */
    info = (WGPUAdapterInfo)WGPU_ADAPTER_INFO_INIT;
    if (G.AdapterGetInfo(G.adapter, &info) == WGPUStatus_Success) {
        const char* backend = "gpu";
        /* A string view out of the API carries an explicit length and need
           not be NUL-terminated, so the precision matters. */
        size_t len = info.device.length;
        if (len == WGPU_STRLEN) len = info.device.data ? strlen(info.device.data) : 0;
        if (len > 96) len = 96;
        switch (info.backendType) {
            case WGPUBackendType_Metal: backend = "Metal"; break;
            case WGPUBackendType_Vulkan: backend = "Vulkan"; break;
            case WGPUBackendType_D3D12: backend = "D3D12"; break;
            case WGPUBackendType_OpenGL: backend = "OpenGL"; break;
            default: break;
        }
        snprintf(G.name, sizeof G.name, "wgpu:%s:%.*s%s", backend, (int)len,
                 info.device.data ? info.device.data : "?", G.f16 ? ":f16" : "");
        /* The four string views are allocations; freeing takes the struct by
           value, and G.name already has its own copy. */
        G.AdapterInfoFreeMembers(info);
    } else {
        snprintf(G.name, sizeof G.name, "wgpu:unknown%s", G.f16 ? ":f16" : "");
    }
    return 1;
}


/* ---- pipelines, built once per weight width ----------------------------- */
/* None of this depends on the grid, so it is built on the first call and
   kept. Calibration makes 100-1000 solves per site; recompiling the shader
   for each would be the dominant cost. Two sets, because f16 weights are
   right for one flux form and not for the other -- see `use_f16_for`. */
static struct Pipes {
    int ready;
    WGPUShaderModule module;
    WGPUBindGroupLayout bgl;
    WGPUPipelineLayout layout;
    WGPUComputePipeline sweep, reduce_a, reduce_b;
} P[2];

static WGPUComputePipeline make_pipeline(struct Pipes* q, const char* entry) {
    WGPUComputePipelineDescriptor cpd;
    memset(&cpd, 0, sizeof cpd);
    cpd.layout = q->layout;
    cpd.compute.module = q->module;
    cpd.compute.entryPoint = sv(entry);
    return G.DeviceCreateComputePipeline(G.device, &cpd);
}

/* Whether the weights may be f16 for this flux form.

   They may for FLUX_SMOLUCHOWSKI and they may not for FLUX_ITO, and the
   reason is the physics rather than the arithmetic. Rounding the six
   neighbour weights and rebuilding the self term from the rounded ones keeps
   the row sum exact, so the step is `(self + sum w) * p0` plus a term in the
   *differences* between neighbours -- the rounding cancels to first order
   exactly when neighbouring voxels hold nearly the same density.

   Smoluchowski's equilibrium is uniform whatever the mobility, and its
   density varies by 0.06% between neighbours on the benchmark field. Ito's
   equilibrium is `p ~ 1/D`, so its density inherits the mobility's own
   voxel-to-voxel variation -- measured at 12%, against the mobility's 11.9%.
   There is nothing for the compensation to cancel against, and the raw f16
   error survives: 2.3e-4 after 32 steps, where Smoluchowski is at 4.5e-7.

   `IMP_BFF_GPU_F16=off` turns it off everywhere, for comparing. */
static int use_f16_for(int smoluchowski) {
    static int checked = 0, allowed = 1;
    if (!checked) {
        const char* e = getenv("IMP_BFF_GPU_F16");
        checked = 1;
        allowed = !(e && (!strcmp(e, "off") || !strcmp(e, "0") || !strcmp(e, "no")));
    }
    return G.f16 && allowed && smoluchowski;
}

static int build_pipelines(int use16) {
    struct Pipes* q = &P[use16 ? 1 : 0];
    WGPUShaderSourceWGSL wgsl;
    WGPUShaderModuleDescriptor smd;
    WGPUBindGroupLayoutEntry entries[8];
    WGPUBindGroupLayoutDescriptor bgld;
    WGPUPipelineLayoutDescriptor pld;
    char* source;
    const char* prefix;
    size_t len;
    int i;

    if (q->ready) return 1;

    /* The shader is written once and typed here: `WT` is the weight type, and
       f16 halves the biggest read in the kernel. */
    prefix = use16 ? "enable f16;\nalias WT = f16;\n" : "alias WT = f32;\n";
    len = strlen(prefix) + sizeof kDiffusionWgsl;
    source = (char*)malloc(len);
    if (!source) return 0;
    memcpy(source, prefix, strlen(prefix));
    memcpy(source + strlen(prefix), kDiffusionWgsl, sizeof kDiffusionWgsl);

    memset(&wgsl, 0, sizeof wgsl);
    wgsl.chain.sType = WGPUSType_ShaderSourceWGSL;
    wgsl.code = sv(source);
    memset(&smd, 0, sizeof smd);
    smd.nextInChain = (WGPUChainedStruct*)&wgsl;
    q->module = G.DeviceCreateShaderModule(G.device, &smd);
    free(source);
    if (!q->module) return 0;

    memset(entries, 0, sizeof entries);
    for (i = 0; i < 8; ++i) {
        entries[i].binding = (uint32_t)i;
        entries[i].visibility = WGPUShaderStage_Compute;
        if (i == 5) {
            entries[i].buffer.type = WGPUBufferBindingType_Uniform;
        } else if (i == 1 || i == 4 || i == 6) {
            entries[i].buffer.type = WGPUBufferBindingType_Storage;
        } else {
            entries[i].buffer.type = WGPUBufferBindingType_ReadOnlyStorage;
        }
    }
    memset(&bgld, 0, sizeof bgld);
    bgld.entryCount = 8;
    bgld.entries = entries;
    q->bgl = G.DeviceCreateBindGroupLayout(G.device, &bgld);
    if (!q->bgl) return 0;

    memset(&pld, 0, sizeof pld);
    pld.bindGroupLayoutCount = 1;
    pld.bindGroupLayouts = &q->bgl;
    q->layout = G.DeviceCreatePipelineLayout(G.device, &pld);
    if (!q->layout) return 0;

    q->sweep = make_pipeline(q, "sweep");
    q->reduce_a = make_pipeline(q, "reduce_a");
    q->reduce_b = make_pipeline(q, "reduce_b");
    if (!q->sweep || !q->reduce_a || !q->reduce_b) return 0;
    q->ready = 1;
    return 1;
}

/* ---- the propagation ---------------------------------------------------- */
static WGPUBuffer make_buffer(uint64_t size, WGPUBufferUsage usage) {
    WGPUBufferDescriptor bd;
    memset(&bd, 0, sizeof bd);
    bd.size = (size + 3u) & ~(uint64_t)3u;    /* copies want a multiple of 4 */
    bd.usage = usage;
    return G.DeviceCreateBuffer(G.device, &bd);
}

static int read_back(WGPUBuffer staging, size_t bytes, void* dst) {
    struct MapReply mr = {0, 0};
    WGPUBufferMapCallbackInfo mci;
    const void* mapped;
    int guard;
    memset(&mci, 0, sizeof mci);
    mci.mode = WGPUCallbackMode_AllowProcessEvents;
    mci.callback = on_map;
    mci.userdata1 = &mr;
    G.BufferMapAsync(staging, WGPUMapMode_Read, 0, bytes, mci);
    for (guard = 0; !mr.done && guard < 100000; ++guard) {
        G.DevicePoll(G.device, 1, NULL);
        G.InstanceProcessEvents(G.instance);
    }
    if (!mr.ok) return 0;
    mapped = G.BufferGetConstMappedRange(staging, 0, bytes);
    if (!mapped) return 0;
    memcpy(dst, mapped, bytes);
    G.BufferUnmap(staging);
    return 1;
}

static int propagate(const double* cur, const double* d, const double* decay,
                     const double* bounds, int ng, int flux_form,
                     int n_steps, int n_out, double* out_fluorescence,
                     int n_fluorescence, double* out_density) {
    const int smol = (flux_form == FLUX_SMOLUCHOWSKI);
    const size_t n = (size_t)ng * (size_t)ng * (size_t)ng;
    long str[6];
    uint32_t* idx = NULL;
    float* self_f = NULL;
    float* w32 = NULL;
    uint16_t* w16 = NULL;
    float* cur_f = NULL;
    float* trace_f = NULL;
    size_t n_act = 0, t, wbytes;
    int rc = 1, q, step, which = 0, slot = 0, use16;
    struct Pipes* pipes;
    uint32_t params[4];
    WGPUBuffer bA = NULL, bB = NULL, bW = NULL, bS = NULL, bPart = NULL;
    WGPUBuffer bParams = NULL, bTrace = NULL, bIdx = NULL, bAll = NULL;
    WGPUBuffer stTrace = NULL, stDens = NULL;
    WGPUBindGroup g0 = NULL, g1 = NULL;

    if (!G.device || ng < MIN_NG || n_steps < MIN_STEPS || n_out <= 0) return 1;
    if (n_fluorescence != n_steps / n_out + 1) return 1;
    use16 = use_f16_for(smol);
    if (!build_pipelines(use16)) return 1;
    pipes = &P[use16 ? 1 : 0];

    str[0] = -(long)ng * ng; str[1] = (long)ng * ng;
    str[2] = -(long)ng;      str[3] = (long)ng;
    str[4] = -1;             str[5] = 1;

    idx = (uint32_t*)malloc(n * sizeof(uint32_t));
    cur_f = (float*)calloc(n, sizeof(float));
    if (!idx || !cur_f) goto done;
    {
        long ix, iy, iz;
        for (ix = 1; ix + 1 < ng; ++ix)
            for (iy = 1; iy + 1 < ng; ++iy)
                for (iz = 1; iz + 1 < ng; ++iz) {
                    const size_t c = (size_t)((ix * ng + iy) * ng + iz);
                    if (bounds[c] == 0.0) continue;
                    idx[n_act++] = (uint32_t)c;
                    /* Outside the domain the density is zero and stays zero:
                       the CPU sweep clears the whole grid each step, and the
                       ping-pong here writes only active voxels, so anything
                       left outside would live for ever. */
                    cur_f[c] = (float)cur[c];
                }
    }
    if (n_act == 0) goto done;         /* nothing to propagate; let the CPU say so */

    self_f = (float*)malloc(n_act * sizeof(float));
    wbytes = use16 ? 6 * n_act * sizeof(uint16_t) : 6 * n_act * sizeof(float);
    if (use16) w16 = (uint16_t*)malloc(wbytes + 4);
    else       w32 = (float*)malloc(wbytes + 4);
    trace_f = (float*)malloc((size_t)n_fluorescence * sizeof(float));
    if (!self_f || (!w16 && !w32) || !trace_f) goto done;

    for (t = 0; t < n_act; ++t) {
        const size_t c = idx[t];
        const double d0 = d[c], k0 = decay[c];
        double out_sum = 0.0, in_sum = 0.0, rowsum, rounded = 0.0;
        double a[6];
        for (q = 0; q < 6; ++q) {
            const size_t m = (size_t)((long)c + str[q]);
            const double bm = bounds[m];
            /* Smoluchowski takes the interface mobility, Ito the neighbour's;
               the self term is the sum that leaves, which is not the same sum
               for the two forms. */
            a[q] = smol ? 0.5 * (d0 + d[m]) * bm : d[m] * bm;
            out_sum += smol ? a[q] : d0 * bm;
            in_sum += a[q];
        }
        rowsum = k0 * (1.0 - out_sum) + k0 * in_sum;
        for (q = 0; q < 6; ++q) {
            const double wq = k0 * a[q];
            if (use16) {
                const uint16_t h = f16_of((float)wq);
                float back;
                w16[(size_t)q * n_act + t] = h;
                /* Read the rounded value back and build the self term from
                   *it*: the row sum is what a long propagation conserves, and
                   rounding the seven numbers independently is worth three
                   orders of magnitude of error. */
                {
                    uint32_t sign = (uint32_t)(h & 0x8000u) << 16;
                    uint32_t e = (h >> 10) & 0x1Fu, mant = h & 0x3FFu, bits;
                    if (e == 0) {
                        if (mant == 0) bits = sign;
                        else {
                            int sh = 0;
                            while (!(mant & 0x400u)) { mant <<= 1; ++sh; }
                            mant &= 0x3FFu;
                            bits = sign | ((uint32_t)(113 - sh) << 23) | (mant << 13);
                        }
                    } else if (e == 0x1F) {
                        bits = sign | 0x7F800000u | (mant << 13);
                    } else {
                        bits = sign | ((e + 112u) << 23) | (mant << 13);
                    }
                    memcpy(&back, &bits, 4);
                }
                rounded += (double)back;
            } else {
                const float f = (float)wq;
                w32[(size_t)q * n_act + t] = f;
                rounded += (double)f;
            }
        }
        self_f[t] = (float)(rowsum - rounded);
    }

    params[0] = (uint32_t)ng;
    params[1] = (uint32_t)n;
    params[2] = (uint32_t)n_act;
    params[3] = (uint32_t)NPART;

    bA = make_buffer(n * 4, WGPUBufferUsage_Storage | WGPUBufferUsage_CopySrc |
                            WGPUBufferUsage_CopyDst);
    bB = make_buffer(n * 4, WGPUBufferUsage_Storage | WGPUBufferUsage_CopySrc |
                            WGPUBufferUsage_CopyDst);
    bW = make_buffer(wbytes, WGPUBufferUsage_Storage | WGPUBufferUsage_CopyDst);
    bS = make_buffer(n_act * 4, WGPUBufferUsage_Storage | WGPUBufferUsage_CopyDst);
    bIdx = make_buffer(n_act * 4, WGPUBufferUsage_Storage | WGPUBufferUsage_CopyDst);
    bPart = make_buffer(NPART * 4, WGPUBufferUsage_Storage);
    bTrace = make_buffer(4, WGPUBufferUsage_Storage | WGPUBufferUsage_CopySrc);
    bAll = make_buffer((uint64_t)n_fluorescence * 4,
                       WGPUBufferUsage_Storage | WGPUBufferUsage_CopySrc |
                       WGPUBufferUsage_CopyDst);
    bParams = make_buffer(16, WGPUBufferUsage_Uniform | WGPUBufferUsage_CopyDst);
    stTrace = make_buffer((uint64_t)n_fluorescence * 4,
                          WGPUBufferUsage_MapRead | WGPUBufferUsage_CopyDst);
    stDens = make_buffer(n * 4, WGPUBufferUsage_MapRead | WGPUBufferUsage_CopyDst);
    if (!bA || !bB || !bW || !bS || !bIdx || !bPart || !bTrace || !bAll ||
        !bParams || !stTrace || !stDens) goto done;

    G.QueueWriteBuffer(G.queue, bA, 0, cur_f, n * 4);
    G.QueueWriteBuffer(G.queue, bW, 0, use16 ? (void*)w16 : (void*)w32, wbytes);
    G.QueueWriteBuffer(G.queue, bS, 0, self_f, n_act * 4);
    G.QueueWriteBuffer(G.queue, bIdx, 0, idx, n_act * 4);
    G.QueueWriteBuffer(G.queue, bParams, 0, params, 16);

    {
        WGPUBindGroupEntry e[8];
        WGPUBindGroupDescriptor bgd;
        WGPUBuffer order0[8], order1[8];
        int k;
        order0[0] = bA; order0[1] = bB;
        order1[0] = bB; order1[1] = bA;
        order0[2] = order1[2] = bW;
        order0[3] = order1[3] = bS;
        order0[4] = order1[4] = bPart;
        order0[5] = order1[5] = bParams;
        order0[6] = order1[6] = bTrace;
        order0[7] = order1[7] = bIdx;
        memset(&bgd, 0, sizeof bgd);
        bgd.layout = pipes->bgl;
        bgd.entryCount = 8;
        bgd.entries = e;
        memset(e, 0, sizeof e);
        for (k = 0; k < 8; ++k) {
            e[k].binding = (uint32_t)k;
            e[k].buffer = order0[k];
            e[k].size = (k == 5) ? 16 : WGPU_WHOLE_SIZE;
        }
        g0 = G.DeviceCreateBindGroup(G.device, &bgd);
        for (k = 0; k < 8; ++k) e[k].buffer = order1[k];
        g1 = G.DeviceCreateBindGroup(G.device, &bgd);
        if (!g0 || !g1) goto done;
    }

    {
        const uint32_t groups = (uint32_t)((n_act + 255) / 256);
        WGPUBindGroup g[2];
        g[0] = g0; g[1] = g1;
        for (step = 0; step < n_steps; ) {
            WGPUCommandEncoder enc = G.DeviceCreateCommandEncoder(G.device, NULL);
            const int here = (n_steps - step < BATCH) ? (n_steps - step) : BATCH;
            int k = 0;
            if (!enc) goto done;
            while (k < here) {
                int run;
                if ((step + k) % n_out == 0) {
                    /* The report is its own pass: a copy orders between
                       passes, not inside one, and a queue write would run
                       immediately rather than in dispatch order. */
                    WGPUComputePassEncoder cp =
                        G.CommandEncoderBeginComputePass(enc, NULL);
                    G.ComputePassEncoderSetPipeline(cp, pipes->reduce_a);
                    G.ComputePassEncoderSetBindGroup(cp, 0, g[which], 0, NULL);
                    G.ComputePassEncoderDispatchWorkgroups(cp, NPART, 1, 1);
                    G.ComputePassEncoderSetPipeline(cp, pipes->reduce_b);
                    G.ComputePassEncoderDispatchWorkgroups(cp, 1, 1, 1);
                    G.ComputePassEncoderEnd(cp);
                    G.ComputePassEncoderRelease(cp);
                    G.CommandEncoderCopyBufferToBuffer(enc, bTrace, 0, bAll,
                                                           (uint64_t)slot * 4, 4);
                    ++slot;
                }
                run = here - k;
                if (run > n_out - ((step + k) % n_out))
                    run = n_out - ((step + k) % n_out);
                {
                    WGPUComputePassEncoder cp =
                        G.CommandEncoderBeginComputePass(enc, NULL);
                    int s;
                    G.ComputePassEncoderSetPipeline(cp, pipes->sweep);
                    for (s = 0; s < run; ++s) {
                        G.ComputePassEncoderSetBindGroup(cp, 0, g[which], 0, NULL);
                        G.ComputePassEncoderDispatchWorkgroups(cp, groups, 1, 1);
                        which ^= 1;
                    }
                    G.ComputePassEncoderEnd(cp);
                    G.ComputePassEncoderRelease(cp);
                }
                k += run;
            }
            {
                WGPUCommandBuffer cb = G.CommandEncoderFinish(enc, NULL);
                G.QueueSubmit(G.queue, 1, &cb);
                G.CommandBufferRelease(cb);
                G.CommandEncoderRelease(enc);
            }
            step += here;
        }
        {
            WGPUCommandEncoder enc = G.DeviceCreateCommandEncoder(G.device, NULL);
            WGPUComputePassEncoder cp;
            WGPUCommandBuffer cb;
            if (!enc) goto done;
            cp = G.CommandEncoderBeginComputePass(enc, NULL);
            G.ComputePassEncoderSetPipeline(cp, pipes->reduce_a);
            G.ComputePassEncoderSetBindGroup(cp, 0, g[which], 0, NULL);
            G.ComputePassEncoderDispatchWorkgroups(cp, NPART, 1, 1);
            G.ComputePassEncoderSetPipeline(cp, pipes->reduce_b);
            G.ComputePassEncoderDispatchWorkgroups(cp, 1, 1, 1);
            G.ComputePassEncoderEnd(cp);
            G.ComputePassEncoderRelease(cp);
            G.CommandEncoderCopyBufferToBuffer(enc, bTrace, 0, bAll,
                                                   (uint64_t)slot * 4, 4);
            G.CommandEncoderCopyBufferToBuffer(enc, bAll, 0, stTrace, 0,
                                                   (uint64_t)n_fluorescence * 4);
            G.CommandEncoderCopyBufferToBuffer(enc, which ? bB : bA, 0,
                                                   stDens, 0, (uint64_t)n * 4);
            cb = G.CommandEncoderFinish(enc, NULL);
            G.QueueSubmit(G.queue, 1, &cb);
            G.CommandBufferRelease(cb);
            G.CommandEncoderRelease(enc);
        }
    }

    G.DevicePoll(G.device, 1, NULL);
    if (!read_back(stTrace, (size_t)n_fluorescence * 4, trace_f)) goto done;
    {
        int r;
        for (r = 0; r < n_fluorescence; ++r)
            out_fluorescence[r] = (double)trace_f[r];
    }
    if (out_density) {
        float* dens = (float*)malloc(n * 4);
        if (!dens) goto done;
        if (read_back(stDens, n * 4, dens)) {
            size_t c;
            for (c = 0; c < n; ++c) out_density[c] = (double)dens[c];
        } else {
            free(dens);
            goto done;
        }
        free(dens);
    }
    rc = 0;

done:
    if (g0) G.BindGroupRelease(g0);
    if (g1) G.BindGroupRelease(g1);
    if (bA) G.BufferRelease(bA);
    if (bB) G.BufferRelease(bB);
    if (bW) G.BufferRelease(bW);
    if (bS) G.BufferRelease(bS);
    if (bIdx) G.BufferRelease(bIdx);
    if (bPart) G.BufferRelease(bPart);
    if (bTrace) G.BufferRelease(bTrace);
    if (bAll) G.BufferRelease(bAll);
    if (bParams) G.BufferRelease(bParams);
    if (stTrace) G.BufferRelease(stTrace);
    if (stDens) G.BufferRelease(stDens);
    free(idx); free(cur_f); free(self_f); free(w16); free(w32); free(trace_f);
    return rc;
}

/* ---- the door ----------------------------------------------------------- */
IMPBFF_EXPORT const struct ImpBffComputeBackend* imp_bff_compute_backend(
        const char* argument) {
    char err[256];
    err[0] = '\0';
    if (G.device) return &gBackend;                   /* already open */
    if (!argument || !*argument) return NULL;         /* no wgpu to run on */
    if (!resolve(argument, err, sizeof err)) return NULL;
    if (!open_device(err, sizeof err)) return NULL;
    gBackend.abi = IMPBFF_COMPUTE_BACKEND_ABI;
    gBackend.name = G.name;
    gBackend.propagate = propagate;
    return &gBackend;
}
