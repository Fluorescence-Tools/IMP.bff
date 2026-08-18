# Is the quenching model identifiable from one decay? (PRD-110 stage 0)

Recorded by `python benchmark/quenching_identifiability.py --fit`.
T4 lysozyme 3GUN, chain A residue 132 CB, AV1 (20 Å / 0.5 Å / 3.5 Å), grid
29³ at **1.5 Å**, 2 180 accessible voxels.

> **Correction, 2026-08-18.** This page originally said 2.5 Å. It was 1.5 Å --
> the builder's default. `compute_av` *writes* `simulation_grid_resolution` into
> `source_info` from its `disc_step` argument (`fret/av.py:167`) rather than
> reading it back, so the benchmark's `--resolution` flag was silently ignored
> and every run used the default. The numbers below are unaffected and
> reproduce exactly -- only the label was wrong. The flag now reaches
> `disc_step` and the grid spacing is asserted against it, so this cannot recur. Field model, τ₀ = 4 ns, decay sampled at
64 points over 25 ns, Poisson noise for a decay peaking at 10⁴ counts.

**Answer: no.** Two of five directions carry essentially all the information.
This is [PRD-110](../prds/prd-110.md)'s stage-0 gate, and it fails.

## Fisher information of F(t) with respect to θ

Relative parameter changes, so units do not decide the answer. Eigen-directions
strongest first:

| λ | direction | reading |
|---|---|---|
| 2.46 × 10⁵ | `+0.93 slow_factor − 0.36 rC` | measured very well |
| 3.75 × 10² | `+0.91 rC + 0.37 slow_factor` | measured |
| 6.2 × 10⁻² | `−0.98 kQ_scale` | ~invisible |
| 3.2 × 10⁻³ | `−0.76 free_diffusion + 0.63 contact_distance` | invisible |
| 1.7 × 10⁻⁵ | `−0.77 contact_distance − 0.62 free_diffusion` | invisible |

**Condition number 1.4 × 10¹⁰.** Every individual parameter is poorly
determined — the smallest one-parameter uncertainty is `slow_factor` at ±168 %,
and `contact_distance` reaches ±1.9 × 10⁴ %.

The ordering is not an artefact of the parameter point. Four different
parameterisations and starting vectors were tried during this work; the
condition number ranged over 1.4 × 10¹⁰ – 3.1 × 10¹⁵ and **never fell below
seven orders of magnitude**.

## The fit says the same thing, out loud

A least-squares fit through the solver, to synthetic data generated at known θ
with a deliberately displaced start:

| | evaluations | wall clock |
|---|---|---|
| `scipy.least_squares`, 5 parameters | **676** | **271 s** |

| parameter | true | found | error |
|---|---|---|---|
| `slow_factor` | 0.985 | 0.945 | **−4.0 %** |
| `rC` | 1.500 | 1.583 | **+5.5 %** |
| `contact_distance` | 6.500 | 7.472 | +15.0 % |
| `free_diffusion` | 8.000 | 3.813 | **−52.3 %** |
| `kQ_scale` | 1.000 | 0.307 | **−69.3 %** |

**It fits the data without recovering the parameters**, and it recovers exactly
the two the Fisher analysis calls measured. That is the practical face of a
10¹⁰ condition number: the optimiser slides along the flat directions and stops
wherever it happens to be.

The 676 evaluations / 271 s is also the baseline any surrogate has to beat —
and beating it is worth nothing while the answer is this.

## Three findings about the model, found on the way

Each was a surprise, and each is a property of the model rather than of this
measurement.

### 1. The equilibrium occupancy has a closed form: `p ∝ 1/D`

The solver discretises the flux as `d[i]·p[i] − d[j]·p[j]`, which is
`∂p/∂t = ∇²(Dp)` — the **Itô** form, not `∇·(D∇p)`. Their stationary states
differ exactly when `D` varies in space: `∇·(D∇p)` relaxes to a *uniform*
distribution, `∇²(Dp)` to `D p = const`, i.e. `p ∝ 1/D`.

Verified against the iterative solver to a maximum relative deviation of
**1.1 × 10⁻¹³**. So `GridDiffusionSolver.equilibrium()` was spending tens of
thousands of iterations converging to something available in closed form — and
on this site it *did not converge at all*: still drifting after 40 000
iterations, peak-to-mean climbing 13.6 → 25.1 → 37.8 → 60.4 → 75.1. Every
gradient measured through it was the gradient of an unconverged field.

`equilibrium_occupancy()` now returns it directly, and
`DynamicAccessibleVolume.occupancy` uses it.

**This is a convention, not a derivation.** Itô against Stratonovich against
the isothermal convention is a real modelling choice for diffusion in a mobility
gradient, it was inherited from ChiSurf, and it is not small: it is the whole
reason the dye accumulates in the contact shell.

### 2. `slow_factor` is only meaningful within a whisker of 1.0

It is applied **once per contacting atom**, so it compounds. At a 6.5 Å contact
distance on this site a voxel sees up to 45 atoms:

| `slow_factor` | 0.985 | 0.95 | 0.9 | 0.6 |
|---|---|---|---|---|
| `0.985ⁿ` at n = 45 | 0.50 | 0.10 | 9 × 10⁻³ | **1 × 10⁻¹⁰** |

The median voxel has **zero** contacting atoms, so the field is bimodal: free
at `D` almost everywhere, and exponentially frozen in a thin shell. ChiSurf's
default of 0.985 is the sane end; anything much below it is a different model.

### 3. At a small `slow_factor` the decay becomes *exactly* independent of `D`

Measured before the range was corrected: over an 80× span of `free_diffusion`
(0.5 → 40 Å²/ns) the decay moved by **9.5 × 10⁻⁹**, including at `D = 0.5`
where the diffusion length (5 Å) is far below the AV extent (43.5 Å) — so this
is not the well-mixed limit.

The mechanism is findings 1 and 2 together: with the mobility spanning ten
orders of magnitude, `p ∝ 1/D` puts essentially all the population in the
frozen shell, where it cannot redistribute; and scaling `D` scales `p_eq` and
the diffusion coefficient in a way that cancels. Worth knowing before anyone
tries to fit `D` to a decay.

### 4. Gradients need a θ-independent time step

The obvious choice — set `t_step` from the stability limit `dg²/(6 D_max)` —
makes the integrator's discretisation error a function of the parameter being
differentiated, and that error lands in the derivative:

| finite-difference step | 2 % | 5 % | 10 % | 20 % |
|---|---|---|---|---|
| `‖dF/dD‖`, adaptive `t_step` | 0.0273 | 0.0255 | 0.0193 | 0.0173 |
| `‖dF/dD‖`, fixed `t_step` | 0.01452 | 0.01453 | 0.01458 | 0.01479 |

58 % drift against 1.9 %, and the adaptive step **inflated the sensitivity to
`D` by about a factor of two**. Any gradient-based calibration — surrogate or
not — has to hold the step fixed across θ.

## What this means for PRD-110

The gate says stop, and it says why: a faster forward model fits an
under-determined problem faster, which is not progress. The question that
replaces it is not *how do we accelerate this fit* but **what else has to be
measured** for the parameters to exist as quantities at all — decays at several
sites jointly, anisotropy for the mobility, or FCS. Reducing the model to the
two combinations one decay does determine is the other honest option.
