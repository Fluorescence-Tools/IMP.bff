# Is the quenching model identifiable across sites? (PRD-111 stage 0)

Recorded by `python benchmark/kq_sensitivity_analysis.py --resolution 1.5 --fit`.
Six chain-A sites on T4 lysozyme 3GUN, AV1 (20 Å / 0.5 Å / 3.5 Å), grid 29³ at
**1.5 Å** at site A132 (the AV sets the extent, so the other sites differ). Field model, τ₀ = 4 ns, decay sampled at 64 points over 25 ns, Poisson
noise for a decay peaking at 10⁴ counts — **the same θ₀, noise model and time
axis as [PRD-110 stage 0](quenching_identifiability.md)**, so the two are
directly comparable. Site A132 reproduces PRD-110's condition number
`1.425 × 10¹⁰` to four figures, which is what makes the comparison one.

**Answer: yes.** [PRD-110](../prds/prd-110.md) found θ unrecoverable from one
decay. Six decays sharing one θ recover four of the five parameters to within a
few percent. This is [PRD-111](../prds/prd-111.md)'s stage-0 gate, and it passes.

## Why more sites could help at all, and why it is not obvious

θ is **global** — the same dye, the same chemistry, the same mobility, wherever
it is attached. The geometry is **per site**. So `N` decays are `N` datasets
sharing one parameter vector, and the joint Fisher information is the sum,
`F = Σ_s J_sᵀ J_s`.

Summing is not automatically progress. Six copies of *one* site give `6 F_s`:
identical eigenvectors, identical condition number, error bars better by √6 and
the invisible combination still invisible. The sum beats its terms **only if the
blind directions rotate with the geometry**. That is a question about protein
geometry, not about statistics, and it has to be measured.

## The sites, and what each one is sensitive to

Chosen to span quencher coverage. `dD/dS` is the sensitivity to
`free_diffusion` relative to the sensitivity to `slow_factor` — the mobility's
share of the site's information.

| site | note | voxels | AV within 6 Å of a quencher | median dye–quencher | `dD/dS` |
|---|---|---|---|---|---|
| A132.CB | PRD-110 reference | 2 180 | 22.7 % | 8.55 Å | 0.0141 |
| A124.CB | TRP126 at 3.5 Å — contact | 3 484 | 75.9 % | 3.41 Å | 0.0106 |
| A93.CB | CYS97 6.1, TRP158 7.1 — crowded | 3 494 | 70.4 % | 3.54 Å | 0.0138 |
| A19.CB | TYR18/TYR25 at ~6 Å — many, medium | 3 229 | 65.0 % | 4.24 Å | 0.0161 |
| A65.CB | HIS31 at 9.7 Å — sparse, distant | 3 087 | 33.4 % | 7.85 Å | 0.0074 |
| A53.CB | CYS54 at 6.4 Å, otherwise bare | 4 731 | 27.2 % | 9.73 Å | 0.0066 |

### The hypothesis was wrong, and usefully so

Stated in the PRD before measuring: *diffusion only matters when the dye has to
travel* — a quencher in the contact shell should quench where the dye already
sits, leaving the mobility invisible, while a distant quencher makes reaching it
the rate.

**The ordering is the opposite.** The sites with the *most* quencher coverage
(A19, A93, A132 at `dD/dS` ≈ 0.014–0.016) weight diffusion more than the sparse,
distant ones (A65, A53 at 0.0066–0.0074); the correlation with median
dye–quencher distance is **−0.56**, not positive.

The mechanism it actually follows: diffusion enters when quenching is fast
enough to **deplete**. When a large fraction of the volume quenches, the
surviving population is the part that has not yet reached a quencher, and the
mobility controls how fast that reservoir is drained. When only a small patch
quenches, the decay is nearly single-exponential at a rate set by the
equilibrium overlap, and transport never becomes rate-limiting. Distance to the
quencher is not the variable; **quenched fraction is.**

That matters for choosing sites: to measure dye mobility, label where the dye is
*strongly* quenched, not where it must travel.

## Do the blind directions rotate?

Principal angles between the two weakest directions of each pair of sites.
0° means both sites are blind the same way; 90° means one sees what the other
cannot.

| | A124 | A93 | A19 | A65 | A53 |
|---|---|---|---|---|---|
| **A132** | 0.81 / 56.3 | 2.68 / 41.9 | 4.31 / **88.4** | 14.6 / 31.3 | 7.84 / 36.6 |
| **A124** | | 0.78 / 17.3 | 7.96 / 50.7 | 2.37 / 74.4 | 2.01 / 82.6 |
| **A93** | | | 3.94 / 64.6 | 1.58 / 57.2 | 1.16 / 69.2 |
| **A19** | | | | 4.87 / 65.6 | 10.9 / 54.3 |
| **A65** | | | | | 0.88 / 37.2 |

The second angle reaches **88.4°** (A132 against A19): two sites whose second
blind directions are very nearly complementary. The first angle stays small
(0.8–14.6°) everywhere, so there is one direction *every* site is largely blind
to — identified below.

> **A rank trap, noted because the first attempt fell into it.** Two subspaces
> of dimension `r` in `Rⁿ` are forced to intersect once `2r > n`. With five
> parameters, comparing rank-**3** blind subspaces reports a 0.00° first angle
> for every pair regardless of geometry — arithmetic wearing the costume of a
> result. The default is now rank 2, and the script refuses to be read otherwise.

The angle table is a diagnostic; this is the measurement:

| | |
|---|---|
| sum of the per-site smallest eigenvalues | **0.0456** |
| smallest eigenvalue of the joint matrix | **152.7** |
| **gain from rotation** | **3 350 ×** |

If every site were blind in the same direction, the joint matrix would be about
as blind as the sum of the parts. It is 3 350× better. The rotation is real and
it is the entire effect — no ratio-of-ratios required.

## The joint Fisher information

| λ | direction | |
|---|---|---|
| 2.69 × 10⁷ | `−1.00 slow_factor` | |
| 3.50 × 10⁴ | `−0.98 rC − 0.18 kQ_scale` | |
| 1.38 × 10³ | `+0.93 free_diffusion + 0.26 contact_distance − 0.25 kQ_scale` | **now visible** |
| 5.78 × 10² | `−0.96 contact_distance + 0.27 free_diffusion` | **now visible** |
| 1.53 × 10² | `−0.95 kQ_scale − 0.23 free_diffusion + 0.20 rC` | weakest — the shared blindness |

**Condition number 1.76 × 10⁵**, against `1.96 × 10⁸` for the best single site
and `1.43 × 10¹⁰` for A132 — three orders below the best of its own terms.

Marginal uncertainties: `slow_factor` ±0.2 %, `rC` ±1.7 %, `free_diffusion`
±3.3 %, `contact_distance` ±4.1 %, `kQ_scale` ±7.7 %. All five below 50 %.

The direction that stays weakest is `kQ_scale`: the absolute PET rate scale.
Every site is relatively blind to it for the same reason at all of them — raising
`kQ` and lengthening the dwell in the contact shell do the same thing to the
decay, and no geometry undoes that. Multi-site takes it from invisible to ±7.7 %
but it remains the last one to be believed.

## The fit, which is the confirmation and also the caveat

Same standard PRD-110 was held to: eigenvalues make a claim, an optimiser
walking the real surface has to back it. Synthetic data at known θ, deliberately
displaced start, one θ across all six decays.

| | evaluations | forward solves | wall clock |
|---|---|---|---|
| joint `scipy.least_squares`, 5 parameters | **97** | **582** | **137 s** |

| parameter | true | start | one site (PRD-110) | six sites |
|---|---|---|---|---|
| `slow_factor` | 0.985 | 0.900 | −4.0 % | **−2.8 %** |
| `rC` | 1.500 | 2.400 | +5.5 % | **−0.2 %** |
| `free_diffusion` | 8.000 | 12.800 | −52.3 % | **+2.5 %** |
| `kQ_scale` | 1.000 | 0.600 | −69.3 % | **+4.3 %** |
| `contact_distance` | 6.500 | 10.400 | +15.0 % | **−33.4 %** |

The two parameters a single decay lost entirely — `kQ_scale` at −69 % and
`free_diffusion` at −52 % — come back to a few percent. Note also that the joint
fit converged in **97 evaluations against 676**: a well-conditioned problem is
not merely more correct, it is cheaper, which is the opposite of what PRD-110's
surrogate was for.

### `contact_distance` is not mis-fitted — it is degenerate, and the fit proves it

The Fisher matrix promises ±4.1 %; the fit misses by −33.4 %. The first reading
was that this is a discretisation artefact — `contact_distance` is a hard
threshold, so on a grid `F(t)` moves only when the cutoff crosses an actual
voxel–atom distance, and the derivative of a staircase describes the step it
straddles. **That reading is wrong, and two checks kill it.**

**The fit is statistically perfect.** Reduced χ² = **1.016** at 1.5 Å and
**1.012** at 2.0 Å, on 384 points and 5 parameters. The optimiser did not fail
to converge; it found a *second minimum that fits exactly as well.*

**The second minimum is reproducible, not noise.** Across resolutions it lands
in the same place, and always paired the same way:

| | `contact_distance` | `slow_factor` | reduced χ² |
|---|---|---|---|
| truth | 6.500 | 0.9850 | — |
| fit at 1.5 Å | 4.328 | 0.9574 | 1.016 |
| fit at 2.0 Å | 4.505 | 0.9604 | 1.012 |

And a coarser grid recovers `contact_distance` slightly *better* (−30.7 % against
−33.4 %), which is the opposite of what a discretisation artefact does.

**The two parameter sets are the same physical model.** Building the mobility
field both ways on A132 at 2.0 Å:

| | |
|---|---|
| correlation of `log D` over accessible voxels | **0.927** |
| median mobility (both) | 8.0 Å²/ns — the free value, so the median voxel contacts nothing |
| overlap of the equilibrium occupancy `p ∝ 1/D` | **98.7 %** |

`slow_factor` is applied **once per contacting atom**, and `contact_distance`
sets how many atoms contact. Shrinking the radius reduces the count; raising the
per-atom slowing compensates, and `slow_factor^n(contact_distance)` comes out the
same. The decay sees the field, and **the field is determined to about a
percent.** The pair that parameterises it is not.

**Consequences, and they are not the ones the artefact reading implied.**

* The two are **one parameter, not two.** No number of sites separates them —
  six did not, and the degeneracy is exact rather than merely stiff.
* **Marginal error bars from `JᵀJ` must not be quoted for either.** `slow_factor`
  reports ±0.2 %, the tightest of all five, and is *also* displaced (0.957
  against 0.985) — in exactly the compensating direction. A Fisher matrix is a
  local linearisation and cannot see a second minimum; here it is confidently
  wrong about its best-determined parameter.
* **The model should be re-parameterised**, replacing the pair with the quantity
  the data actually determines — the mobility field, or a smooth contact shell
  with one width parameter. That also makes the model differentiable in a
  parameter it currently is not.
* The three that *are* recovered — `free_diffusion` +2.5 %, `kQ_scale` +4.3 %,
  `rC` −0.2 % — are recovered honestly, and they are the ones a single decay lost
  entirely.

## Does the result depend on the grid?

Repeated end to end at 2.0 Å (1 158 accessible voxels at A132 against 2 180):

| | 1.5 Å | 2.0 Å |
|---|---|---|
| joint condition number | 1.76 × 10⁵ | 3.04 × 10⁵ |
| best single site | 1.96 × 10⁸ | 5.46 × 10⁸ |
| sum of per-site `λ_min` | 0.0456 | 0.0161 |
| joint `λ_min` | 152.7 | 90.6 |
| **gain from rotation** | **3 350 ×** | **5 623 ×** |
| largest principal angle | 88.4° | 88.4° |
| `free_diffusion` recovery | +2.5 % | +2.7 % |
| `kQ_scale` recovery | +4.3 % | +0.5 % |
| `rC` recovery | −0.2 % | +0.2 % |

Same conclusion, same order of magnitude, same failures. The result is not a
grid artefact.

> **A correction that had to be made first.** Both this page and PRD-110's
> originally recorded "2.5 Å". Every run was in fact at the builder's 1.5 Å
> default: `compute_av` *writes* `simulation_grid_resolution` into `source_info`
> from its `disc_step` argument rather than reading it back, so the benchmark's
> `--resolution` flag was inert and the resolution robustness check above had
> never actually been possible. The numbers were always internally consistent —
> only the label was wrong. The flag now reaches `disc_step`, the benchmark
> asserts the grid spacing it got, and `compute_av` raises when `source_info`
> declares a resolution that disagrees with `disc_step` instead of silently
> discarding it.

## What this means

* **Calibrating the PET model requires multi-site labelling.** That is a real
  experimental cost, and it is now a measured requirement rather than a
  preference. Four decays may do; one provably does not.
* **The right sites are the strongly quenched ones**, not the ones where the dye
  has far to travel — and they should be chosen for *complementary* blindness,
  which is computable in advance from the geometry alone, before any measurement
  (A132 with A19 was the best pair here at 88.4°).
* **`kQ_scale` remains the weakest direction** even jointly (±7.7 %), and
  **`contact_distance` and `slow_factor` are one parameter, not two** — six
  sites do not separate them, and the fit reaches reduced χ² = 1.02 at a
  `contact_distance` 33 % away from truth by building the same mobility field.
  Re-parameterise rather than collect more data.
* **None of this is yet evidence about the world.** The data are synthetic, from
  the same model being fitted, so this measures the model's internal geometry
  under a perfect model and nothing else. PRD-111 stage 1 — a per-site nuisance
  parameter, then real decays — is where that assumption gets tested, and it is
  the assumption most likely to be false.
