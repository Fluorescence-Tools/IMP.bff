---
type: validation
title: "The network through the same door: 12x, and 6x of it was one transpose"
description: MlpCore's forward pass routed through the compute door as a whole network rather than per GEMM, with a WGSL kernel behind it. 12.4x the eight-threaded CPU on a 4-128x3-1 net at 400k rows, 5.4x on a narrow one. The first working version managed 2.4x; transposing the weights on the host so that adjacent threads read adjacent addresses was worth 5.9x of the rest, and one command buffer per batch instead of per layer fixed the small-batch case. Deviation 3e-7 to 9e-7, which is f32. Three validation errors found on the way, each of which aborts the process rather than returning.
resource: /Users/tpeulen/dev/imp.bff
tags: [validation, imp.bff, performance, gpu, webgpu, neural-net]
timestamp: '2026-09-08T00:00:00Z'
---
# The network on the GPU

`include/internal/MlpCore.h` has been vendored and unused: header-only,
std-only, with a GEMM template policy written so the batch products could be
routed elsewhere. `IMP::bff::NeuralNet` is the face that uses it, and
`gpu/mlp.wgsl` is the elsewhere.

## The door takes the whole network, not one GEMM

The GEMM policy is the obvious place to hook an accelerator and it is the
wrong one. A door at GEMM granularity moves the activations across it once per
layer, and for the shape this is for -- a network evaluated per voxel, so
hundreds of thousands of rows against a few dozen units -- that is over a
hundred megabytes each way against about a millisecond of arithmetic. The
library already learned this once: porting the diffusion kernel to C++ and
leaving the loop in Python made the suite 40× slower.

So `ImpBffMlpForwardFn` takes the layer shapes, the weights and the whole
batch, and gives back the outputs. One crossing per `predict()`, however deep
the network. The compute-door ABI went to 2 for it; a plugin built against 1
is declined rather than read one field short.

## What it does

| net | rows | multiply-accumulates | GPU | CPU, 8 threads | |
|---|---|---|---|---|---|
| 4-64-64-2 | 1 000 | 4.5 M | — | 1.0 ms | declined, the CPU is quicker |
| 4-64-64-2 | 10 000 | 45 M | 3.6 ms | 10.0 ms | 2.8× |
| 4-64-64-2 | 100 000 | 448 M | 23.7 ms | 101 ms | 4.3× |
| 4-64-64-2 | 400 000 | 1.8 G | 75.7 ms | 411 ms | **5.4×** |
| 4-128×3-1 | 1 000 | 33 M | 2.8 ms | 6.6 ms | 2.3× |
| 4-128×3-1 | 10 000 | 334 M | 9.1 ms | 65 ms | 7.2× |
| 4-128×3-1 | 100 000 | 3.3 G | 55.9 ms | 663 ms | 11.9× |
| 4-128×3-1 | 400 000 | 13.4 G | 219 ms | 2 711 ms | **12.4×** |

Deviation against the CPU's f64: 3 × 10⁻⁷ to 9 × 10⁻⁷, which is f32 and what
the diffusion kernel sees. The backend declines below 16 M
multiply-accumulates, which is below every case that wins here and above the
one that does not.

## Two changes worth 5× between them

The first working version managed 2.4× at the largest size — 13.4 GMAC in
1 287 ms is about 21 GFLOP/s on a device that does thousands. Two changes:

**Transposing the weights, worth 5.9×.** A thread computes one output element,
so adjacent threads differ in the output unit. With the weights stored
`n_out × n_in`, as the model holds them, those threads read addresses `n_in`
apart and every load is scattered. Transposed to `n_in × n_out` on the host —
once per call, against thousands of reads per row — adjacent threads read
adjacent addresses. Nothing else changed.

**One command buffer per batch instead of per layer, worth 2.7× on small
batches** (1 000 rows: 4.25 ms → 1.58 ms) and nothing on large ones. The
uniform that carries a layer's shape is written by `wgpuQueueWriteBuffer`,
which runs when it is called rather than in the order the dispatches were
recorded — the same trap the diffusion kernel hit. The first version worked
around it by submitting and waiting per layer, which is a device round trip
per layer. Giving each layer its own 32-byte uniform buffer and its own bind
group makes the whole network one submission.

It is still only ~122 GFLOP/s, so the kernel is far from the machine and
workgroup-memory tiling is the obvious next thing. It has not been tried.

## Three ways to abort the process, all found here

wgpu validation failures are **not** return codes. They reach a Rust panic in
a callback that cannot unwind, and the process aborts — no exception, no
chance for the CPU path to take over. Each of these was found by hitting it:

- **A dispatch may not exceed 65 535 workgroups in a dimension.** One thread
  per output element and 100 000 rows through a 64-unit layer wants 100 000.
  The shader now uses a grid-stride loop and the host caps the count. The
  diffusion kernel declines rather than exceeding it; at 256 voxels per
  workgroup it would take a grid of ng > 250 to get there.
- **A storage binding has a ceiling**, 128 MiB by default, which 400 000 rows
  through a 128-unit layer passes. The rows of a forward pass do not see each
  other, so the batch is split into as many as fit rather than refused. The
  limit is read from the device, not assumed.
- **`WGPUBindGroupEntry.size` defaults to `WGPU_WHOLE_SIZE`, not 0**, so a
  `memset`-ed entry is a zero-length binding. Several fields in this API
  revision are like that; `WGPUBindGroupLayoutEntry` also grew a
  `bindingArraySize` between `visibility` and `buffer` that shifts everything
  after it if missed.

The lesson for the plugin generally: **a backend cannot decline what it has
already recorded.** Anything the device might reject has to be checked on the
host before it is put in a command buffer.

## What is tested

`test/test_neural_net.py` checks each of the seven activations against an
independent numpy forward pass, that the accelerated and CPU answers agree,
and that a small batch stays on the CPU. `test/test_wgsl_header.py` fails if
`MlpCore.h` grows an eighth activation, because the shader dispatches on seven
and a missing one would silently return the identity.
