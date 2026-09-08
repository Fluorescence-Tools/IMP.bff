// The diffusion stencil, in the form measured fastest on an M1 Pro.
//
// Read okf/validation/gpu_diffusion_is_worth_it.md before changing anything
// here; every line below is the result of a measurement.
//
//  * Only the active voxels are dispatched. `idx` is their list, and a third
//    to two fifths of the cube is outside the domain.
//  * The stencil weights are precomputed on the host. `d`, `decay` and
//    `bounds` do not change over a propagation, so neither do the
//    coefficients: `nxt = self*p0 + sum w_q * cur[m]`. Eighteen scattered
//    reads per voxel become six.
//  * There is no neighbour table. The six neighbours of `c` are `c +/- 1`,
//    `c +/- ng`, `c +/- ng*ng`: compaction changes which voxels are visited,
//    not how the grid is numbered.
//  * `w` is f16 where the device has `shader-f16` (the host writes `alias WT`
//    accordingly). The self term stays f32 *and is computed from the rounded
//    weights*, because what the propagation conserves is the row sum
//    `self + sum w == decay`; rounding the seven numbers independently costs
//    three orders of magnitude of accuracy.
//  * The population sum is two-stage. One workgroup over every active voxel
//    was 16% of the whole run.

struct Params { ng: u32, n: u32, n_act: u32, npart: u32 };

@group(0) @binding(0) var<storage, read>       cur:   array<f32>;
@group(0) @binding(1) var<storage, read_write> nxt:   array<f32>;
@group(0) @binding(2) var<storage, read>       w:     array<WT>;
@group(0) @binding(3) var<storage, read>       self_: array<f32>;
@group(0) @binding(4) var<storage, read_write> part:  array<f32>;
@group(0) @binding(5) var<uniform>             p:     Params;
@group(0) @binding(6) var<storage, read_write> trace: array<f32>;
@group(0) @binding(7) var<storage, read>       idx:   array<u32>;

@compute @workgroup_size(256)
fn sweep(@builtin(global_invocation_id) gid: vec3<u32>) {
    let t = gid.x;
    if (t >= p.n_act) { return; }
    let c = idx[t];
    let na = p.n_act;
    let s2 = p.ng * p.ng;
    var acc = self_[t] * cur[c];
    acc = acc + f32(w[t])         * cur[c - s2];
    acc = acc + f32(w[na + t])    * cur[c + s2];
    acc = acc + f32(w[2u * na + t]) * cur[c - p.ng];
    acc = acc + f32(w[3u * na + t]) * cur[c + p.ng];
    acc = acc + f32(w[4u * na + t]) * cur[c - 1u];
    acc = acc + f32(w[5u * na + t]) * cur[c + 1u];
    nxt[c] = acc;
}

var<workgroup> scratch: array<f32, 256>;

fn fold(lid: u32) {
    workgroupBarrier();
    var s = 128u;
    loop {
        if (s == 0u) { break; }
        if (lid < s) { scratch[lid] = scratch[lid] + scratch[lid + s]; }
        workgroupBarrier();
        s = s / 2u;
    }
}

@compute @workgroup_size(256)
fn reduce_a(@builtin(local_invocation_id) lid: vec3<u32>,
            @builtin(workgroup_id) wid: vec3<u32>) {
    var acc = 0.0;
    var i = wid.x * 256u + lid.x;
    let stride = p.npart * 256u;
    loop {
        if (i >= p.n_act) { break; }
        acc = acc + cur[idx[i]];
        i = i + stride;
    }
    scratch[lid.x] = acc;
    fold(lid.x);
    if (lid.x == 0u) { part[wid.x] = scratch[0]; }
}

@compute @workgroup_size(256)
fn reduce_b(@builtin(local_invocation_id) lid: vec3<u32>) {
    var acc = 0.0;
    var i = lid.x;
    loop {
        if (i >= p.npart) { break; }
        acc = acc + part[i];
        i = i + 256u;
    }
    scratch[lid.x] = acc;
    fold(lid.x);
    if (lid.x == 0u) { trace[0] = scratch[0]; }
}
