// One dense layer, y = activation(W x + b), a thread per output element.
//
// One dispatch per layer, ping-ponging two activation buffers, so a forward
// pass crosses the plugin boundary once however deep the network is. The
// uniform changes between layers and a queue write runs immediately rather
// than in dispatch order, so each layer is its own submission -- there are
// only a handful of them, against tens of thousands of rows.
//
// The reduction over n_in is a plain loop in a register. For the shapes this
// is for -- many rows, a few dozen units -- the parallelism is in the rows,
// and a workgroup reduction per output would be less of it, not more.

struct Params {
    n_rows: u32,
    n_in: u32,
    n_out: u32,
    act: u32,
    w_off: u32,
    b_off: u32,
    stride: u32,
    pad0: u32,
};

@group(0) @binding(0) var<storage, read>       xin:  array<f32>;
@group(0) @binding(1) var<storage, read_write> xout: array<f32>;
// Transposed by the host to n_in x n_out. With the weights stored n_out x
// n_in as the model holds them, adjacent threads -- which differ in the
// output unit -- read addresses n_in apart, and the load is scattered. Stored
// this way they read adjacent addresses. Measured worth ~3x.
@group(0) @binding(2) var<storage, read>       wgt:  array<f32>;
@group(0) @binding(3) var<storage, read>       bia:  array<f32>;
@group(0) @binding(4) var<uniform>             p:    Params;

// The codes are IMP::bff::internal::Activation in declaration order.
fn activate(v: f32, act: u32) -> f32 {
    if (act == 0u) { return v; }                                  // identity
    if (act == 1u) { return max(v, 0.0); }                        // relu
    if (act == 2u) { return tanh(v); }                            // tanh
    if (act == 3u) { return 1.0 / (1.0 + exp(-v)); }              // logistic
    if (act == 4u) {
        // log1p(exp(v)), written so that a large v does not overflow before
        // the logarithm brings it back.
        if (v > 20.0) { return v; }
        return log(1.0 + exp(v));                                 // softplus
    }
    if (act == 5u) { return v / (1.0 + exp(-v)); }                // silu
    if (act == 6u) { return sin(v); }                             // sin
    return v;
}

// A grid-stride loop rather than one thread per element: a dispatch may not
// exceed 65 535 workgroups in a dimension, and a batch of 100 000 rows through
// a 64-unit layer wants 100 000 of them. Exceeding it is a validation error
// that aborts the process rather than returning, so the host caps the count
// and each thread walks the rest.
@compute @workgroup_size(64)
fn layer(@builtin(global_invocation_id) gid: vec3<u32>) {
    let total = p.n_rows * p.n_out;
    var t = gid.x;
    loop {
        if (t >= total) { break; }
        let r = t / p.n_out;
        let o = t % p.n_out;
        let xrow = r * p.n_in;
        var acc = bia[p.b_off + o];
        for (var i = 0u; i < p.n_in; i = i + 1u) {
            acc = acc + wgt[p.w_off + i * p.n_out + o] * xin[xrow + i];
        }
        xout[t] = activate(acc, p.act);
        t = t + p.stride;
    }
}
