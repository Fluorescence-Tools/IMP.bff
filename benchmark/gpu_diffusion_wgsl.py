"""The diffusion stencil in WGSL, held against the CPU kernel.

This is the shader the GPU backend will run and the harness that says whether
it is worth running: same fields, same steps, the fluorescence trace compared
value by value and the wall clock taken as a minimum over repetitions.

It is a *benchmark*, not a backend. It drives wgpu through `wgpu-py` because
that package already carries the wgpu-native this machine has, which made the
number available today rather than after the C plugin was written. The WGSL
below is what the plugin will use; the plumbing around it is what the plugin
will replace.

Run it with a build that has `IMP.bff` importable and `wgpu` installed:

    python benchmark/gpu_diffusion_wgsl.py

Two traps found while writing it, both worth knowing before writing the same
thing in C:

* **`queue.write_buffer` runs immediately; dispatches run when the encoder is
  submitted.** Updating a uniform between recorded dispatches therefore does
  not interleave with them -- every dispatch sees the last value written. The
  report index has to travel through something the encoder orders, which here
  is `copy_buffer_to_buffer` between passes.
* **A dynamic storage-buffer offset did not take effect in this wgpu.** Every
  reduction wrote to slot 0, with a constant as well as with the sum. The
  copies above are the way around it, and they cost nothing: forty of them
  against four thousand sweeps.
"""
import time, numpy as np, wgpu, IMP.bff as b

SHADER = """
struct Params { ng: u32, n: u32, slot: u32, pad: u32 };
@group(0) @binding(0) var<storage, read>       cur:    array<f32>;
@group(0) @binding(1) var<storage, read_write> nxt:    array<f32>;
@group(0) @binding(2) var<storage, read>       d:      array<f32>;
@group(0) @binding(3) var<storage, read>       decay:  array<f32>;
@group(0) @binding(4) var<storage, read>       bnds:   array<f32>;
@group(0) @binding(5) var<uniform>             p:      Params;
@group(0) @binding(6) var<storage, read_write> trace:  array<f32>;

@compute @workgroup_size(64)
fn sweep(@builtin(global_invocation_id) gid: vec3<u32>) {
    let c = gid.x;
    if (c >= p.n) { return; }
    let ng = p.ng;
    let iz = c % ng;
    let iy = (c / ng) % ng;
    let ix = c / (ng * ng);
    if (ix == 0u || iy == 0u || iz == 0u ||
        ix + 1u >= ng || iy + 1u >= ng || iz + 1u >= ng) { nxt[c] = 0.0; return; }
    if (bnds[c] == 0.0) { nxt[c] = 0.0; return; }
    let p0 = cur[c];
    let d0 = d[c];
    var flux = 0.0;
    let nb = array<u32, 6>(c - ng * ng, c + ng * ng, c - ng, c + ng, c - 1u, c + 1u);
    for (var q = 0u; q < 6u; q = q + 1u) {
        let m = nb[q];
        flux = flux + 0.5 * (d0 + d[m]) * (p0 - cur[m]) * bnds[m];
    }
    nxt[c] = (p0 - flux) * decay[c];
}

var<workgroup> scratch: array<f32, 256>;
@compute @workgroup_size(256)
fn reduce(@builtin(local_invocation_id) lid: vec3<u32>) {
    var acc = 0.0;
    var i = lid.x;
    loop {
        if (i >= p.n) { break; }
        acc = acc + cur[i];
        i = i + 256u;
    }
    scratch[lid.x] = acc;
    workgroupBarrier();
    var s = 128u;
    loop {
        if (s == 0u) { break; }
        if (lid.x < s) { scratch[lid.x] = scratch[lid.x] + scratch[lid.x + s]; }
        workgroupBarrier();
        s = s / 2u;
    }
    if (lid.x == 0u) { trace[0] = scratch[0]; }
}
"""

def fields(ng, seed=3):
    rng = np.random.default_rng(seed)
    i = np.arange(ng) - ng // 2
    X, Y, Z = np.meshgrid(i, i, i, indexing="ij")
    inside = (X**2 + Y**2 + Z**2) <= (ng // 2 - 2) ** 2
    cur = (inside / inside.sum()).astype(np.float64)
    d = (0.10 + 0.05 * rng.random((ng, ng, ng))) * inside
    decay = np.exp(-(0.02 + 0.20 * (X**2 + Y**2 + Z**2) / (ng**2)) * 0.005) * inside + (~inside)
    return (cur.ravel(), d.ravel(), decay.ravel(), inside.ravel().astype(np.float64))


def gpu_propagate(dev, ng, cur, d, decay, bnds, n_steps, n_out, batch=500, wgs=64):
    """Der ganze Lauf auf dem Gerät: ein Pass je Stapel, ein dynamischer
    Offset für den Bericht.

    Der Bericht kann nicht mit einem Kopierbefehl aus dem Pass heraus
    gerettet werden -- Kopien stehen zwischen Pässen, nicht in ihnen -- und
    ein `queue.write_buffer` liefe sofort statt in der Reihenfolge der
    Dispatches. Also schreibt die Reduktion an einen dynamischen Offset, den
    `set_bind_group` beim Aufzeichnen mitgibt.
    """
    n = ng ** 3
    n_rep = n_steps // n_out + 1
    U = wgpu.BufferUsage
    ALIGN = 256                       # Mindestausrichtung für dynamische Offsets

    def buf(a):
        return dev.create_buffer_with_data(
            data=np.ascontiguousarray(a, dtype=np.float32),
            usage=U.STORAGE | U.COPY_SRC)

    b_a, b_b = buf(cur), buf(np.zeros(n))
    b_d, b_de, b_bo = buf(d), buf(decay), buf(bnds)
    # Die Reduktion schreibt immer nach b_tr[0]; ein Kopierbefehl *zwischen*
    # den Pässen rettet den Wert an seinen Platz. Dynamische Offsets wären
    # der elegantere Weg und wurden zuerst versucht: in dieser wgpu-Fassung
    # kam nur der erste Bericht an, auch mit einer Konstante statt der Summe.
    b_tr = dev.create_buffer(size=4, usage=U.STORAGE | U.COPY_SRC)
    b_all = dev.create_buffer(size=4 * n_rep, usage=U.STORAGE | U.COPY_SRC | U.COPY_DST)
    b_pa = dev.create_buffer(size=16, usage=U.UNIFORM | U.COPY_DST)
    dev.queue.write_buffer(b_pa, 0, np.array([ng, n, 0, 0], dtype=np.uint32).tobytes())

    mod = dev.create_shader_module(code=SHADER.replace("@workgroup_size(64)",
                                                       "@workgroup_size(%d)" % wgs))
    ro = {"type": wgpu.BufferBindingType.read_only_storage}
    rw = {"type": wgpu.BufferBindingType.storage}
    un = {"type": wgpu.BufferBindingType.uniform}
    dyn = {"type": wgpu.BufferBindingType.storage}
    bgl = dev.create_bind_group_layout(entries=[
        {"binding": i, "visibility": wgpu.ShaderStage.COMPUTE, "buffer": t}
        for i, t in enumerate((ro, rw, ro, ro, ro, un, dyn))])
    play = dev.create_pipeline_layout(bind_group_layouts=[bgl])
    p_sweep = dev.create_compute_pipeline(layout=play,
                                          compute={"module": mod, "entry_point": "sweep"})
    p_red = dev.create_compute_pipeline(layout=play,
                                        compute={"module": mod, "entry_point": "reduce"})

    def group(a, b_):
        return dev.create_bind_group(layout=bgl, entries=[
            {"binding": 0, "resource": {"buffer": a, "offset": 0, "size": a.size}},
            {"binding": 1, "resource": {"buffer": b_, "offset": 0, "size": b_.size}},
            {"binding": 2, "resource": {"buffer": b_d, "offset": 0, "size": b_d.size}},
            {"binding": 3, "resource": {"buffer": b_de, "offset": 0, "size": b_de.size}},
            {"binding": 4, "resource": {"buffer": b_bo, "offset": 0, "size": b_bo.size}},
            {"binding": 5, "resource": {"buffer": b_pa, "offset": 0, "size": 16}},
            {"binding": 6, "resource": {"buffer": b_tr, "offset": 0, "size": 4}},
        ])
    g = (group(b_a, b_b), group(b_b, b_a))
    wg = (n + wgs - 1) // wgs

    t0 = time.perf_counter()
    which, slot, step = 0, 0, 0

    def report(enc, which_, slot_):
        # Ein eigener Pass je Bericht: in einem gemeinsamen Pass landete nur
        # der erste Schreibvorgang, obwohl jeder Dispatch seinen eigenen
        # dynamischen Offset bekam. Vierzig kleine Pässe kosten nichts gegen
        # viertausend Sweeps.
        cp = enc.begin_compute_pass()
        cp.set_pipeline(p_red)
        cp.set_bind_group(0, g[which_])
        cp.dispatch_workgroups(1)
        cp.end()
        enc.copy_buffer_to_buffer(b_tr, 0, b_all, 4 * slot_, 4)

    while step < n_steps:
        enc = dev.create_command_encoder()
        n_here = min(batch, n_steps - step)
        k = 0
        while k < n_here:
            if (step + k) % n_out == 0:
                report(enc, which, slot)
                slot += 1
            # so viele Sweeps am Stück, wie bis zum nächsten Bericht passen
            run = min(n_here - k, n_out - ((step + k) % n_out))
            cp = enc.begin_compute_pass()
            cp.set_pipeline(p_sweep)
            for _ in range(run):
                cp.set_bind_group(0, g[which])
                cp.dispatch_workgroups(wg)
                which ^= 1
            cp.end()
            k += run
        dev.queue.submit([enc.finish()])
        step += n_here
    enc = dev.create_command_encoder()
    report(enc, which, slot)
    dev.queue.submit([enc.finish()])
    trace = np.frombuffer(dev.queue.read_buffer(b_all), dtype=np.float32)[:slot + 1].copy()
    return trace, time.perf_counter() - t0


dev = wgpu.utils.get_default_device()
print("Adapter:", dev.adapter.info["device"], "|", dev.adapter.info["backend_type"])
print("%-24s %10s %10s %8s %12s" % ("Gitter / Schritte", "CPU (ms)", "GPU (ms)", "Faktor", "Abweichung"))
REPS = 3
for ng, n_steps in ((41, 4000), (61, 4000), (81, 4000)):
    cur, d, decay, bnds = fields(ng)
    n_out = 100
    cpu_t, gpu_t, dev_max = [], [], 0.0
    for _ in range(REPS):
        t0 = time.perf_counter()
        fl_cpu, _ = b.diffusion_propagate(cur.tolist(), d.tolist(), decay.tolist(),
                                          bnds.tolist(), ng, 0, n_steps, n_out)
        cpu_t.append(time.perf_counter() - t0)
        fl_gpu, gt = gpu_propagate(dev, ng, cur, d, decay, bnds, n_steps, n_out)
        gpu_t.append(gt)
        fl_cpu = np.asarray(fl_cpu)
        k = min(len(fl_cpu), len(fl_gpu))
        rel = np.abs(fl_gpu[:k] - fl_cpu[:k]) / np.maximum(np.abs(fl_cpu[:k]), 1e-30)
        dev_max = max(dev_max, rel.max())
    # Minima: unter Last ist das Minimum der belastbare Schätzer, der
    # Mittelwert misst die Nachbarprozesse mit.
    c, g_ = min(cpu_t) * 1e3, min(gpu_t) * 1e3
    print("ng=%-3d %6d Schritte %10.1f %10.1f %7.1fx %12.2e"
          % (ng, n_steps, c, g_, c / g_, dev_max))
