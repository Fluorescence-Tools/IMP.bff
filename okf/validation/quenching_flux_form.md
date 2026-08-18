# Where does the dye sit at equilibrium? (the flux discretisation)

Recorded 2026-08-18, prompted by the owner pointing at the **Haas–Steinberg**
equation. Measurements by `benchmark/kq_sensitivity_analysis.py --flux-form
{smoluchowski,ito}` and `test/quenching/test_quenching_field.py::FluxFormTests`.

**Answer: the inherited discretisation was wrong, and it mattered more than
anything else measured in PRD-110 or PRD-111.** `flux_form="smoluchowski"` is
now the default; `"ito"` reproduces the old behaviour.

## The two operators

The kernel folded in from ChiSurf propagates the flux between neighbouring
voxels as `d[i]·p[i] − d[j]·p[j]`, which discretises

```
∂p/∂t = ∇²(D p)                     stationary state:  p ∝ 1/D
```

The alternative writes the flux as `D_ij·(p[i] − p[j])`:

```
∂p/∂t = ∇·(D ∇p)                    stationary state:  p uniform on the domain
```

They coincide when `D` is uniform and disagree completely when it is not.
PRD-110 recorded this as *"a convention, not a derivation — Itô against
Stratonovich against the isothermal convention is a real choice."* **That was
too generous.** Two things decide it.

### 1. Equilibrium is thermodynamics; mobility is kinetics

A dye slowed by friction near the protein surface, with no attractive
interaction, must still be found uniformly across its accessible volume at
equilibrium — it merely takes longer to get around. Letting a friction field
set the distribution asserts a potential that was never specified. Under
`∇²(Dp)`, halving the local mobility doubles the local population, which is a
statement about free energy made by a parameter that only describes drag.

### 2. Haas–Steinberg is written to avoid exactly this

The canonical treatment of diffusion-modulated FRET on a flexible linker:

```
∂N(r,t)/∂t = −[1/τ_D + k_T(r)]·N(r,t) + D ∂/∂r [ p(r) ∂/∂r ( N(r,t) / p(r) ) ]
```

The diffusion operator is deliberately of the form `D ∂/∂r [ p ∂/∂r (N/p) ]`, so
its stationary state is the **given** `p(r)` for any `D`. `p(r)` comes from the
chain statistics — here, from the accessible volume — and `D` is a *separate*
kinetic parameter fitted against it. The structure is the whole point: the
equilibrium distribution is an input, not an output of the mobility.

QuEst's own `notebooks/04_diffusion_modulated_fret.ipynb` integrates this
equation, but with a constant `D` and a uniform `p(r)` — the one case where both
discretisations coincide, so it never exercised the difference.

## Both closed forms, verified against the kernel

A domain with a 16× mobility contrast down one side, propagated to equilibrium
with no decay:

| | iterated vs closed form | peak / min occupancy |
|---|---|---|
| `smoluchowski` | **0.00e+00** | **1.0000** |
| `ito` | 1.2e-11 | **16.0000** — exactly the `D` ratio |

So `p ∝ 1/D` is correct algebra for the operator that was implemented. The
operator was the problem, not the algebra.

## What actually changes: `slow_factor` was a disguised attraction

Under `∇²(Dp)`, reducing `D` near an atom *concentrates* the dye there — and the
quenchers are exactly where the mobility is reduced. So `slow_factor` was not
modelling friction; it was modelling an attraction to the quenchers, and the
compounding `slow_factor^n_contacts` made that attraction exponentially strong.

Mean donor lifetime on T4L A132, `D = 0.5 Å²/ns`, over the full sensible range of
`slow_factor`:

| | `slow_factor` 0.985 | 0.90 | change |
|---|---|---|---|
| `ito` | 3.4253 ns | 2.3221 ns | **−32.2 %** |
| `smoluchowski` | 3.4920 ns | 3.4988 ns | **+0.2 %** |

**A factor of ~160 in how much the parameter matters.** Under the correct
operator `slow_factor` does what its name says — it slows redistribution — and
because quenching here is slower than transport, that barely changes the decay.

## What this does to the identifiability results

Same six T4L sites, 2.0 Å, same θ₀ and noise as PRD-111 stage 0.

### Joint Fisher information, six sites

| | `ito` | `smoluchowski` |
|---|---|---|
| eigenvalues | 2.76e7, 3.74e4, 1275, 537, 90.6 | 8.21e5, 2.99e4, 735, 128, **2.66** |
| condition number | 3.04 × 10⁵ | 3.09 × 10⁵ |
| `free_diffusion` | ±3.7 % | ±6.7 % |
| `slow_factor` | ±0.2 % | ±2.9 % |
| `contact_distance` | ±4.5 % | **±61.1 %** |
| `kQ_scale` | ±9.9 % | ±8.1 % |
| `rC` | ±2.3 % | ±2.1 % |
| rotation gain from multi-site | 5 623× | 396× |
| largest principal angle | 88.4° | 49.1° |

The leading eigenvalue falls **34×**: a third of the apparent information in the
old model was `slow_factor` acting as an attraction. Note the condition numbers
are nearly identical and mean different things — another reason not to read a
condition number alone.

### The fit, which is where it gets better

| parameter | true | `ito` | `smoluchowski` |
|---|---|---|---|
| `kQ_scale` | 1.000 | +0.5 % | **−0.0 %** |
| `rC` | 1.500 | +0.2 % | **−0.1 %** |
| `slow_factor` | 0.985 | −2.5 % | **+1.1 %** |
| `free_diffusion` | 8.000 | +2.7 % | +10.0 % |
| `contact_distance` | 6.500 | −30.7 % | +69.2 % (runs to its bound) |
| forward solves | | 516 | **312** |
| reduced χ² | | 1.012 | 1.023 |

**The PET chemistry comes back essentially exactly** — `kQ_scale` to −0.0 % and
`rC` to −0.1 %, better than the old model ever managed — in 40 % fewer solves.

## Four earlier findings are now known to be artifacts

* **PRD-110 finding 1 (`p ∝ 1/D` in closed form).** The algebra is right and the
  1.1e-13 verification stands, but it is the equilibrium of the wrong operator.
  Under the default it is uniform on the accessible domain. The *practical*
  half of that finding survives intact: `GridDiffusionSolver.equilibrium()` was
  iterating tens of thousands of steps toward something available in closed
  form, and on a real site was not converging at all.
* **PRD-110 finding 2 (`slow_factor` only meaningful within a whisker of 1.0).**
  It compounds per contacting atom either way, but that only mattered because
  compounding drove an attraction. Under the default, 0.90 changes the lifetime
  by 0.2 %.
* **PRD-110 finding 3 (the decay becomes *exactly* independent of `D`).** That
  was `p ∝ 1/D` freezing the whole population into the contact shell at a
  pathological `slow_factor`. It does not happen under the default.
* **PRD-111 stage 0 finding 1 (`contact_distance` and `slow_factor` are one
  parameter).** They traded exactly because `slow_factor^n(contact_distance)`
  set the strength of the spurious attraction. Under the default they do not
  trade: `slow_factor` is recovered to +1.1 % and `contact_distance` is simply
  **uninformative** (±61 %, its own eigenvector at λ = 2.66, running to its
  bound in the fit). The practical advice is unchanged and better founded —
  **fix `contact_distance`, do not fit it** — but the reason is different.

## What survives, and what has not been re-measured

**Survives.** The multi-site argument itself: the blind directions still rotate
with geometry and the joint matrix is still far better than the sum of its parts
(396× — smaller than 5 623×, still decisive). `kQ_scale` is still the parameter
to be most careful about. `contact_distance` still must not be fitted.

**Not re-measured.** PRD-111 **stage 1 was run entirely under `ito`** — the
free-dye nuisance, the 9.4 σ detectability, the benign 5 % geometry error. Those
numbers are not carried over, and the free-fraction result in particular could
move a lot, since it traded against `free_diffusion` and `free_diffusion`'s role
changes here. Re-running stage 1 under the default is the next thing.

## The capability this gives up, and how to get it back

Real dyes *do* stick to protein surfaces, and occupancy differing from the AV
density is a real effect. ChiSurf's model was reaching for it. **The objection is
not to stickiness but to expressing it as a mobility**, which conflates a free
energy with a friction and makes one parameter do both jobs — badly, since it
forces the two to scale together.

The correct form takes a separate equilibrium `p_eq`:

```
flux_ij = D_ij · p_eq,ij · ( p_i/p_eq,i − p_j/p_eq,j )
```

which reduces to the implemented `D_ij·(p_i − p_j)` when `p_eq` is uniform. So a
sticky dye is representable — as an attractive potential supplied alongside the
mobility, with its own parameter. That is not implemented yet, and until it is,
**the model has no stickiness at all** rather than the wrong kind. Whether the
data can support a separate stickiness parameter is exactly the identifiability
question this line of work is set up to answer.
