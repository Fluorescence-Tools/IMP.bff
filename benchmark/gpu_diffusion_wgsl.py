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



# ---------------------------------------------------------------------------
# The fast variant, and why it is faster
#
# Five things, measured one at a time (okf/validation/gpu_diffusion_is_worth_it.md):
#
# * **Only the voxels that compute.** A third to two fifths of the cube is
#   inside the volume; the rest was being dispatched and then discarded by a
#   branch. An index list of the active voxels is worth 1.3-1.5x -- less than
#   the 2.5x the counts suggest, because the neighbour gathers are scattered
#   and skipping a voxel does not skip a cache line.
# * **A workgroup of 256** rather than 64, worth about 4%.
# * **Precomputed weights.** `d`, `decay` and `bounds` do not change over a
#   propagation, so neither do the stencil coefficients:
#       nxt = [decay*(1 - sum a_q)] * p0 + sum [decay*a_q] * cur[m]
#   Computing the six `a_q` once turns eighteen scattered reads per voxel into
#   six, and the six weights are read in a coalesced layout. Worth 1.6x.
# * **No neighbour table.** The six neighbours of `c` are `c +/- 1`,
#   `c +/- ng` and `c +/- ng*ng` -- the compaction changes which voxels are
#   visited, not how the grid is numbered. The table of six `u32` per voxel
#   was therefore twenty-four bytes of pure redundancy, and it was the largest
#   single read in the kernel. Worth 1.34-1.59x, and exact.
# * **f16 weights, with the self term compensated.** With the table gone the
#   weights are most of what is left to read, and halving them is worth
#   1.09-1.23x. Naively it costs three orders of magnitude of accuracy
#   (7e-7 -> 1e-3), which is why it looked like a dead end; see
#   `half_weights` for the one line that gives all of it back.
# * **A two-stage population sum.** The report ran in one workgroup over every
#   active voxel -- 897 serial iterations per thread at ng = 81, while the rest
#   of the device waited. Forty of them cost 23 ms of 119. Splitting it into
#   NPART partial sums and a second tiny pass costs 3.7 ms, worth 16% of the
#   whole run -- and, unlooked for, two orders of magnitude of accuracy: the
#   deviation at ng = 81 falls from 1.1e-5 to 1.3e-7, because the old
#   accumulator was summing a quarter of a million f32 values in series. Most
#   of what was being blamed on f32 was the reduction, not the stencil.
#
# Together: 13.3x / 13.2x / 14.0x against the eight-threaded CPU at
# ng = 41 / 61 / 81, where the naive kernel managed 5.3x / 4.1x / 3.8x. Note
# which way those two rows run: the naive kernel loses ground as the grid
# grows, the tuned one gains.
#
# Where the remaining time goes, measured with a no-op sweep and with the
# reports switched off (`us` per step, 4 000 steps):
#
#   ng    total   dispatch floor   reports   kernel
#   41   26.5 ms      10.8 (41%)    2.9 ms   ~13 ms
#   81   98.3 ms      11.5 (12%)    3.7 ms   ~83 ms
#
# So the two ends want different things. At ng = 41 the kernel is already
# cheaper than the command stream that launches it, and only fewer dispatches
# (temporal blocking) can help. At ng = 81 the kernel is the cost, roughly
# half of it the six neighbour gathers, which is what workgroup-memory tiling
# would attack.
# ---------------------------------------------------------------------------

SHADER_W = """
ENABLE
struct Params { ng: u32, n: u32, n_act: u32, pad: u32 };
@group(0) @binding(0) var<storage, read>       cur:   array<f32>;
@group(0) @binding(1) var<storage, read_write> nxt:   array<f32>;
@group(0) @binding(2) var<storage, read>       w:     array<WT>;   // 6 je Voxel
@group(0) @binding(3) var<storage, read>       self_: array<f32>;  // decay*(1-Sum w)
@group(0) @binding(4) var<storage, read_write> part:  array<f32>;  // Teilsummen
@group(0) @binding(8) var<storage, read>       nbi:   array<u32>;  // nur GATHER_TABLE
@group(0) @binding(5) var<uniform>             p:     Params;
@group(0) @binding(6) var<storage, read_write> trace: array<f32>;
@group(0) @binding(7) var<storage, read>       idx:   array<u32>;

@compute @workgroup_size(WGS)
fn sweep(@builtin(global_invocation_id) gid: vec3<u32>) {
    let t = gid.x;
    if (t >= p.n_act) { return; }
    let c = idx[t];
    let na = p.n_act;
    var acc = self_[t] * cur[c];
NEIGHBOURS
    nxt[c] = acc;
}

// Die Populationssumme, zweistufig. Einstufig -- eine Arbeitsgruppe ueber
// alle aktiven Voxel -- kostete bei ng=81 23 ms von 119, also 16% des ganzen
// Laufs fuer vierzig Zahlen: 897 Durchlaeufe je Thread auf einem Kern, waehrend
// der Rest des Geraets wartet. Mit NPART Gruppen sind es 3.7 ms.
var<workgroup> scratch: array<f32, 256>;

fn fold(lid: u32) {
    workgroupBarrier();
    var s = 128u;
    loop {
        if (s == 0u) { break; }
        if (lid < s) { scratch[lid] = scratch[lid] + scratch[lid + s]; }
        workgroupBarrier(); s = s / 2u;
    }
}

@compute @workgroup_size(256)
fn reduce_a(@builtin(local_invocation_id) lid: vec3<u32>,
            @builtin(workgroup_id) wid: vec3<u32>) {
    var acc = 0.0;
    var i = wid.x * 256u + lid.x;
    loop { if (i >= p.n_act) { break; } acc = acc + cur[idx[i]]; i = i + NPART * 256u; }
    scratch[lid.x] = acc; fold(lid.x);
    if (lid.x == 0u) { part[wid.x] = scratch[0]; }
}

@compute @workgroup_size(256)
fn reduce_b(@builtin(local_invocation_id) lid: vec3<u32>) {
    var acc = 0.0;
    var i = lid.x;
    loop { if (i >= NPART) { break; } acc = acc + part[i]; i = i + 256u; }
    scratch[lid.x] = acc; fold(lid.x);
    if (lid.x == 0u) { trace[0] = scratch[0]; }
}
"""

# Wie viele Teilsummen. Zwischen 16 und 256 ist der Unterschied Rauschen
# (1.6-4.4 ms bei ng=81); was zaehlt, ist ueberhaupt mehr als eine.
NPART = 64

# Die Nachbarn aus einer Tabelle -- was der erste schnelle Kernel tat.
GATHER_TABLE = """
    for (var q = 0u; q < 6u; q = q + 1u) {
        acc = acc + f32(w[q * na + t]) * cur[nbi[q * na + t]];
    }
"""
# Die Nachbarn aus den Schrittweiten. Die Reihenfolge ist die, in der
# `coefficients` sie ablegt: (-1,0,0) (1,0,0) (0,-1,0) (0,1,0) (0,0,-1) (0,0,1).
GATHER_STRIDE = """
    let s2 = p.ng * p.ng;
    acc = acc + f32(w[t])         * cur[c - s2];
    acc = acc + f32(w[na + t])    * cur[c + s2];
    acc = acc + f32(w[2u*na + t]) * cur[c - p.ng];
    acc = acc + f32(w[3u*na + t]) * cur[c + p.ng];
    acc = acc + f32(w[4u*na + t]) * cur[c - 1u];
    acc = acc + f32(w[5u*na + t]) * cur[c + 1u];
"""


def coefficients(ng, d, decay, bnds):
    """nxt[c] = decay*(p0 - Sum 0.5(d0+dm)(p0-cm)bm)
              = [decay*(1 - Sum a_q)]*p0 + Sum [decay*a_q]*c_m,  a_q = 0.5(d0+dm)bm
    Die Felder stehen fest, also stehen auch die Gewichte fest.

    `nbi` ist die Nachbartabelle. Der Kernel braucht sie nicht -- sie steht
    hier, weil `GATHER_TABLE` misst, was sie kostet."""
    D = d.reshape(ng, ng, ng); B = bnds.reshape(ng, ng, ng); K = decay.reshape(ng, ng, ng)
    interior = np.zeros((ng, ng, ng), bool); interior[1:-1, 1:-1, 1:-1] = True
    act = (B != 0) & interior
    idx = np.flatnonzero(act).astype(np.uint32)
    ix, iy, iz = np.unravel_index(idx, (ng, ng, ng))
    offs = [(-1,0,0),(1,0,0),(0,-1,0),(0,1,0),(0,0,-1),(0,0,1)]
    w = np.zeros((6, len(idx)), np.float64); nbi = np.zeros((6, len(idx)), np.uint32)
    a_sum = np.zeros(len(idx), np.float64)
    d0 = D[ix, iy, iz]; k0 = K[ix, iy, iz]
    for q, (ox, oy, oz) in enumerate(offs):
        jx, jy, jz = ix + ox, iy + oy, iz + oz
        a = 0.5 * (d0 + D[jx, jy, jz]) * B[jx, jy, jz]
        a_sum += a
        w[q] = k0 * a
        nbi[q] = (jx * ng * ng + jy * ng + jz).astype(np.uint32)
    self_ = k0 * (1.0 - a_sum)
    return idx, w, self_, nbi.ravel(), k0


def half_weights(w, self_, k0):
    """Die Gewichte auf f16 runden, ohne die Genauigkeit zu verlieren.

    Was der Lauf ueber tausende Schritte erhaelt, ist die Zeilensumme:
    `self_ + sum w == decay`. Rundet man die sieben Zahlen unabhaengig, ist
    die Summe um ~5e-4 falsch, jeder Schritt legt denselben Fehler nach, und
    die Spur weicht um 1e-3 ab -- tausendmal mehr als in f32.

    Also nicht unabhaengig runden: erst die sechs Gewichte auf f16 bringen,
    dann den Selbstterm *aus den gerundeten* Werten bilden und ihn in f32
    lassen. Die Zeilensumme stimmt dann bis auf die f32-Rundung, und die
    Abweichung der Spur ist die von f32 (gemessen: 1.0e-6 / 3.0e-6 / 1.1e-5
    bei ng = 41 / 61 / 81, gegen 7.2e-7 / 3.0e-6 / 1.1e-5).

    Ein f16-Selbstterm macht es wieder kaputt, und zwar schlimmer als gar
    keine Kompensation (1.8e-3 bei ng=41): die Kompensation schiebt den
    ganzen Zeilenfehler in genau die Zahl, die dann gerundet wird."""
    wh = w.astype(np.float16)
    return wh, (k0 - wh.astype(np.float64).sum(0)).astype(np.float32)


def gpu_weights(dev, ng, cur, d, decay, bnds, n_steps, n_out,
                stride=True, half=False, batch=500, wgs=256):
    n = ng ** 3
    n_rep = n_steps // n_out + 1
    U = wgpu.BufferUsage
    idx, w, self_, nbi, k0 = coefficients(ng, d, decay, bnds)
    n_act = len(idx)
    if half:
        w, self_ = half_weights(w, self_, k0)
    else:
        w, self_ = w.astype(np.float32), self_.astype(np.float32)
    w = w.ravel()

    def buf(a, dt=np.float32):
        a = np.ascontiguousarray(a, dtype=dt)
        if a.nbytes % 4:                       # Pufferlaenge auf 4 Byte runden
            a = np.concatenate([a, np.zeros(1, dt)])
        return dev.create_buffer_with_data(data=a, usage=U.STORAGE | U.COPY_SRC)

    b_a, b_b = buf(cur), buf(np.zeros(n))
    b_w = buf(w, np.float16 if half else np.float32)
    b_s, b_ix = buf(self_), buf(idx, np.uint32)
    b_nb = buf(np.zeros(1) if stride else nbi, np.uint32)
    b_pt = dev.create_buffer(size=4 * NPART, usage=U.STORAGE | U.COPY_SRC)
    b_tr = dev.create_buffer(size=4, usage=U.STORAGE | U.COPY_SRC)
    b_all = dev.create_buffer(size=4 * n_rep, usage=U.STORAGE | U.COPY_SRC | U.COPY_DST)
    b_pa = dev.create_buffer(size=16, usage=U.UNIFORM | U.COPY_DST)
    dev.queue.write_buffer(b_pa, 0, np.array([ng, n, n_act, 0], np.uint32).tobytes())
    code = (SHADER_W.replace("ENABLE", "enable f16;" if half else "")
                    .replace("WT", "f16" if half else "f32")
                    .replace("NEIGHBOURS", GATHER_STRIDE if stride else GATHER_TABLE)
                    .replace("NPART", "%du" % NPART)
                    .replace("WGS", str(wgs)))
    mod = dev.create_shader_module(code=code)
    ro = {"type": wgpu.BufferBindingType.read_only_storage}
    rw = {"type": wgpu.BufferBindingType.storage}
    un = {"type": wgpu.BufferBindingType.uniform}
    bgl = dev.create_bind_group_layout(entries=[
        {"binding": i, "visibility": wgpu.ShaderStage.COMPUTE, "buffer": t}
        for i, t in enumerate((ro, rw, ro, ro, rw, un, rw, ro, ro))])
    play = dev.create_pipeline_layout(bind_group_layouts=[bgl])
    p_sw = dev.create_compute_pipeline(layout=play, compute={"module": mod, "entry_point": "sweep"})
    p_ra = dev.create_compute_pipeline(layout=play, compute={"module": mod, "entry_point": "reduce_a"})
    p_rb = dev.create_compute_pipeline(layout=play, compute={"module": mod, "entry_point": "reduce_b"})
    def group(a, b_):
        r = [a, b_, b_w, b_s, b_pt, b_pa, b_tr, b_ix, b_nb]
        return dev.create_bind_group(layout=bgl, entries=[
            {"binding": i, "resource": {"buffer": x, "offset": 0, "size": x.size}}
            for i, x in enumerate(r)])
    g = (group(b_a, b_b), group(b_b, b_a))
    wg = (n_act + wgs - 1) // wgs
    t0 = time.perf_counter()
    which, slot, step = 0, 0, 0
    while step < n_steps:
        enc = dev.create_command_encoder(); n_here = min(batch, n_steps - step); k = 0
        while k < n_here:
            if (step + k) % n_out == 0:
                cp = enc.begin_compute_pass()
                cp.set_pipeline(p_ra); cp.set_bind_group(0, g[which])
                cp.dispatch_workgroups(NPART)
                cp.set_pipeline(p_rb); cp.dispatch_workgroups(1)
                cp.end()
                enc.copy_buffer_to_buffer(b_tr, 0, b_all, 4 * slot, 4); slot += 1
            run = min(n_here - k, n_out - ((step + k) % n_out))
            cp = enc.begin_compute_pass(); cp.set_pipeline(p_sw)
            for _ in range(run):
                cp.set_bind_group(0, g[which]); cp.dispatch_workgroups(wg); which ^= 1
            cp.end(); k += run
        dev.queue.submit([enc.finish()]); step += n_here
    enc = dev.create_command_encoder()
    cp = enc.begin_compute_pass()
    cp.set_pipeline(p_ra); cp.set_bind_group(0, g[which]); cp.dispatch_workgroups(NPART)
    cp.set_pipeline(p_rb); cp.dispatch_workgroups(1)
    cp.end()
    enc.copy_buffer_to_buffer(b_tr, 0, b_all, 4 * slot, 4)
    dev.queue.submit([enc.finish()])
    tr = np.frombuffer(dev.queue.read_buffer(b_all), np.float32)[:slot + 1].copy()
    return tr, time.perf_counter() - t0


VARIANTS = [
    ("naive",           gpu_propagate, {}),
    ("weights, table",  gpu_weights,   dict(stride=False)),
    ("weights, stride", gpu_weights,   dict(stride=True)),
    ("stride, f16",     gpu_weights,   dict(stride=True, half=True)),
]


def main(reps=5, sizes=((41, 4000), (61, 4000), (81, 4000)), n_out=100):
    ad = wgpu.gpu.request_adapter_sync(power_preference="high-performance")
    half_ok = "shader-f16" in ad.features
    dev = ad.request_device_sync(required_features=["shader-f16"] if half_ok else [])
    print("Adapter:", ad.info["device"], "|", ad.info["backend_type"],
          "| shader-f16:", half_ok)
    # Ohne diese Zeile ist die CPU-Spalte nicht lesbar: derselbe Aufruf ist
    # ein- oder achtthreadig, je nachdem, welcher Build importiert wurde.
    print("CPU: %d thread(s), OpenMP %s"
          % (b.openmp_thread_count(), "on" if b.built_with_openmp() else "off"))
    print("%-8s %-16s %10s %9s %10s %11s"
          % ("Grid", "variant", "time (ms)", "vs naive", "vs CPU", "deviation"))
    for ng, n_steps in sizes:
        cur, d, decay, bnds = fields(ng)
        cpu_t = []
        for _ in range(reps):
            t0 = time.perf_counter()
            fl_cpu, _ = b.diffusion_propagate(cur.tolist(), d.tolist(), decay.tolist(),
                                              bnds.tolist(), ng, 0, n_steps, n_out)
            cpu_t.append(time.perf_counter() - t0)
        fl_cpu = np.asarray(fl_cpu)
        # Minima: under load the minimum is the estimator that means
        # something, the mean measures the neighbours. Interleaved, so that a
        # burst of load does not land on one variant alone.
        cpu = min(cpu_t) * 1e3
        want = [v for v in VARIANTS if half_ok or "f16" not in v[0]]
        ts = {n: [] for n, _, _ in want}; trs = {}
        for _ in range(reps):
            for name, fn, kw in want:
                tr, t = fn(dev, ng, cur, d, decay, bnds, n_steps, n_out, **kw)
                ts[name].append(t); trs[name] = tr
        base = None
        for name, _, _ in want:
            ms = min(ts[name]) * 1e3
            if base is None: base = ms
            tr = trs[name]; k = min(len(tr), len(fl_cpu))
            rel = np.abs(tr[:k] - fl_cpu[:k]) / np.maximum(np.abs(fl_cpu[:k]), 1e-30)
            print("ng=%-5d %-16s %10.1f %8.2fx %9.1fx %11.1e"
                  % (ng, name, ms, base / ms, cpu / ms, rel.max()))


if __name__ == "__main__":
    main()
