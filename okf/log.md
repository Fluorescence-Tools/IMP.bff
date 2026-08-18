# Update Log

## 2026-08-18 (PRD-113 stage 1: `dye/` — the species)

* **A dye had two unrelated descriptions.** The *spectral* one in
  `fret/forster.py` (names, extinction coefficients, quantum yields, curves from
  the bundled tables); the *molecular* one in `cgdye/io/template_cif.py` (atoms,
  transition-dipole atoms, formal charges). Nothing connected them, so "what is
  AlexaFluor 488" had two answers. `IMP.bff.dye.Dye` is now the one answer, and
  it is deliberately **model-independent**: `linker_length` and
  `allowed_sphere_radius` are AV *representation* parameters, and a rotamer
  library has neither.
* **R0 is derived, not supplied.** `forster_radius(donor, acceptor, kappa2, n)`
  takes two `Dye`s and the medium. The old route hard-coded the refractive index
  as the literal `1.4**4` inside the calculation, so nothing could ask what R0
  would be in a different solvent — it is now a parameter (AlexaFluor488/594:
  5.6878 nm at n = 1.4, **5.8856 nm at n = 1.33**). The name-based
  `forster_radius_from_spectra` is kept for the rotamer code and agrees exactly.
* **Names aligned to flrCIF**, checked against `python-ihm`'s FLR model rather
  than guessed. `Dye.name` → `_flr_probe_list.chromophore_name`,
  `chromophore_center_atom` → `_flr_probe_descriptor.chromophore_center_atom`,
  `lifetime` → `_flr_reference_measurement_lifetime.lifetime`, plus
  `reactive_probe_name`/`probe_origin`/`probe_link_type`. flrCIF's word for a dye
  is **probe** (the reagent) and for the fluorescent moiety **chromophore**.
  Recorded as `FLRCIF_ITEMS`, same convention as `fret/fps_schema.py`.
* **flrCIF has no item for quantum yield, extinction coefficient or a spectrum.**
  The nearest, `_flr_fret_calibration_parameters.phi_acceptor`, is an analysis
  calibration value, not a species property. Those three are bff-native and
  marked `None` deliberately rather than by omission.
* **`_flr_reference_measurement_lifetime` is `(species_fraction, lifetime,
  species_name)`** — the lifetime-spectrum output contract for stage 6 already
  exists in the dictionary, named. Worth adopting verbatim there.
* **Two defects found and fixed on the way.** `api.py` had five names appearing
  twice; one (`read_component_template_cif`) was a genuine duplicate, and the
  other four were **not duplicates at all** — they live in a second dict, the
  rename map from export name to the module's real attribute. A regex dedup
  deleted that map and broke four exports; caught by the suite, restored. The
  lesson is the obvious one about regex edits on structured data.
* **A plan assumption corrected**: PRD-113 stage 1 said to move
  `cgdye/io/template_cif.py` into `dye/`. Having read it, it is a *generic
  component* template reader used by `topology/builder`, `analysis/density` and
  scripts — only `read_dye_template_cif` is dye-specific. It stays put and moves
  to `io/` at stage 7; `dye/library.py` calls into it meanwhile.
* Suite **605 passed** (from 599), same single pre-existing `IMP.em` failure.

## 2026-08-18 (PRD-113 stage 0: the instrument layer goes)

* **`IMP.bff` is a forward-model engine** — what `tttrlib` is to ChiSurf for
  photons, `IMP.bff` is for structure. It emits **experiment-neutral** quantities
  (lifetime spectra, rate constants, κ² distributions, distances) and leaves
  convolution, IRF, pileup, counting statistics and fitting against raw data
  outside. Scoring against *structural* data (distances, as `AVNetworkRestraint`
  does) stays. Recorded as [PRD-113](prds/prd-113.md); PRD-112 is absorbed into
  it, its LabelLib stage done and standing.
* **Deleted the TCSPC instrument layer**: `src/Decay*.cpp` (10) +
  `src/PhotonStatistics.cpp`, `include/Decay*.h` (11) +
  `include/internal/PhotonStatistics.h`, `pyext/Decay*.i` (11) and their
  `%include` lines, `test/test_Decay*.py` (9), and
  `pyext/src/spectroscopy/decay.py`. 44 files. It did convolution, pileup,
  linearisation and decay *fitting* — all experiment-side — and its only consumer
  was that one Python module. **Nothing in chisurf, quest, imp-tricks or tttrlib
  referenced any of it**, and no AV/PathMap source included it, so the cut was
  clean.
* **`spectroscopy/kappa2.py` deliberately survives.** It is the opposite kind of
  thing: a forward model producing orientation-factor distributions,
  wobbling-in-a-cone averages and order parameters from structure. It moves to
  `photophysics/` in stage 4, joining `fret/kappa2.py` so that "kappa squared"
  has one home instead of two. `spectroscopy/__init__.py` re-exported the deleted
  module and is now a docstring saying which half went and why —
  `test_tttrlib_is_optional.py` imports that package, so it had to keep
  importing.
* **Verified**: cmake reconfigured, module rebuilt, `DecayCurve` and friends
  absent from `IMP.bff`, `AV`/`PathMap`/`AVNetworkRestraint` intact. Suite **599
  passed**; the single failure, `test_AccessibleVolume::test_access_av_feature`,
  is the pre-existing unrelated `IMP.em` MRC error — the gate is that it neither
  disappears nor gains company, and it did neither.

## 2026-08-18 (solver stability; the rate term)

* **The explicit-step criterion ignored the rate term, which dominates it**
  (`okf/validation/quenching_solver_stability.md`). The update coefficient is
  `1 − 6D dt/dg² − k dt`; `diffusion_stability_limit` validated only the
  diffusion half. At T4L site 19, 2.5 Å, with a 25 Å²/ns ceiling on `D` setting
  the step: diffusion 0.16, **quenching 2.01**, sum 2.17 — past *divergence* —
  and the decay reached **7 × 10³⁶** while the solver called the step safe.
* **The divergence need not look like one.** Site 124 at a sum of 2.10 returned a
  smooth, finite, monotone, plausible decay that was **2.6 % wrong** at 25 ns.
  That is why it is pinned by a test rather than left to inspection.
* **Two thresholds, and they differ by two**: the coefficient goes negative at a
  sum of 1 and the scheme diverges at 2. `diffusion_stability_limit(d_max, dg,
  k_max)` returns the positivity bound, the stricter one, because a negative
  probability density is not an acceptable answer either.
* **The rate is now integrated exactly** — `exp(-k dt)` as a factor rather than
  `1 - k dt` subtracted — so it contributes no stability constraint at all and a
  strongly quenched site costs the same as a weak one. Pinned: a uniform rate
  reproduces `exp(-k t)` to 1e-10 whatever the step.
* **A published claim retracted.** On finding this I said PRD-111 stages 0–1 were
  contaminated. Wrong twice: I read the positivity bound as the divergence bound,
  and I computed the thresholds with the *new harness's* step, 1.6× larger than
  the one those stages used. They ran at sums of 0.92–0.98, and re-measuring
  reproduces them to ~1 % (eigenvalues 8.223e5/2.992e4/736/129/2.67 against
  8.209e5/2.994e4/735/128/2.66). **The divergence belongs to the new harness**:
  tightening the `D` ceiling to save compute raised the step and pushed the
  contact sites past 2.
* **A real hazard removed**: at 2.0 Å the fits could search `kQ_scale` to 10,
  where the sum reaches 8.6 and the old scheme would have diverged. They
  converged near 1.0 and never went there, but nothing was stopping them.

## 2026-08-18 (the flux form, and the estimator)

* **The flux discretisation was wrong, and it changes every earlier quenching
  finding** (`okf/validation/quenching_flux_form.md`). PRD-110 recorded it as
  *"a convention, not a derivation"*; that was too generous. Equilibrium is
  thermodynamics and mobility is kinetics, so a dye slowed by friction with no
  attraction must still be found uniformly across its accessible volume — while
  the inherited `d[i]p[i] − d[j]p[j]` flux gives `p ∝ 1/D`. The field's own
  canonical treatment agrees: the **Haas-Steinberg** equation writes the
  diffusion operator as `D ∂/∂r [ p(r) ∂/∂r ( N/p(r) ) ]` precisely so its
  stationary state is the *given* `p(r)` for any `D`. QuEst's notebook 04
  integrates it, but at constant `D` and uniform `p(r)` — the one case where both
  discretisations coincide.
* **`flux_form="smoluchowski"` is the default; `"ito"` reproduces the inherited
  behaviour.** Both closed forms verified against the kernel at a 16× mobility
  contrast: peak/min occupancy 1.0000 against 16.0000.
* **`slow_factor` was a disguised attraction.** Reducing `D` near an atom
  concentrates the dye there and the quenchers are exactly where `D` is reduced,
  so `slow_factor^n_contacts` made the attraction exponential. Mean lifetime on
  T4L A132 at `D = 0.5` going 0.985 → 0.90: **−32.2 % under `ito`, +0.2 % under
  the default** — a factor of ~160 in how much the parameter matters.
* **Four earlier findings are artifacts and are corrected in place**: PRD-110's
  closed-form `p ∝ 1/D` (right algebra, wrong operator — the practical half
  stands), its *"`slow_factor` only meaningful within a whisker of 1.0"*, its
  *"decay becomes exactly independent of `D`"*, and PRD-111 stage 0's
  *"`contact_distance` and `slow_factor` are one parameter"*. They traded because
  `slow_factor^n` set the strength of the spurious attraction; under the default
  `slow_factor` recovers to +1.1 % and `contact_distance` is simply
  **uninformative** (±61 %, its own eigenvector at λ = 2.66).
* **Where it counts the model got better**: the joint six-site fit now recovers
  `kQ_scale` to −0.0 % and `rC` to −0.1 %, in 312 forward solves against 516. The
  leading Fisher eigenvalue falls 34×, which is the size of the information that
  was really an attraction. The multi-site argument survives (rotation gain 396×).
* **PRD-111 stage 1 was not re-measured** — the free-dye nuisance, the 9.4 σ
  detectability and the benign geometry error all ran under `ito`.
* **The model now has no stickiness at all rather than the wrong kind.** Real
  dyes stick; the correct form takes a separate equilibrium,
  `flux_ij = D_ij p_eq,ij (p_i/p_eq,i − p_j/p_eq,j)`, reducing to what is
  implemented when `p_eq` is uniform. Future work, and its own identifiability
  question.
* **The estimator was `neyman_lsq`, hand-rolled, and it is the wrong one.**
  tttrlib ships the objective registry (`fit_objectives_json()`) and my
  benchmark reproduced `statistics::neyman` down to the `max(1, ·)` clamp — which
  tttrlib itself documents as *"biased low at small counts"*, against
  `poisson_mle` (2I\*) being *"what a TCSPC decay should normally be fitted
  with"*. These decays run to nothing inside the window, so most of the axis is
  that regime. One mitigation, and it flattered the study rather than the
  reverse: σ came from the noiseless model at the true θ, an oracle weighting no
  real fit has.
* **A point estimate is the wrong output**, given a 10⁵ condition number and a
  demonstrated second minimum — the answer is a posterior. Prior information is
  being discarded too: `PET_QUENCHING_REFERENCE` ships as *"starting values meant
  to be calibrated"*, which is a prior, and a bounded uniform search over
  `kQ_scale ∈ [0.05, 10]` throws it away.
* **All of it already exists, in the right layers.** Statistics in tttrlib
  (`twoIstar`, `statistics::{neyman,poisson,pearson,gauss}`); inference in
  ChiSurf (`FitGroup.run(local_first=True)` staged group fitting, `priors.py`
  MAP, `sample.py`/`ensemble.py` MCMC, `diagnostics.py` R̂/ESS/MCSE);
  forward model here. `benchmark/kq_sensitivity_analysis.py` does all three
  inside `imp.bff` and should not grow further — **PRD-111's inference stages
  move to ChiSurf.** The Fisher analysis stands either way, being a property of
  the model and the noise.
* **χ² was never going to show the problem**: 1.01–1.02 in every run, including
  those where `contact_distance` was 30–70 % wrong or pinned at its bound. A
  scalar cannot distinguish an under-determined fit from a misfitting one; the
  residual *shape* can, and nothing here looked at it.

## 2026-08-18 (PRD-111 stage 2 scoped; the acceptor observable)

* **Stage 2 rewritten around the acceptor decay, on the owner's physics
  (2026-08-18): quenching and FRET are spatially uncorrelated and
  non-homogeneous.** Quenching is set by contact with the residues immediately
  around the dye; FRET by `1/R⁶` to an acceptor elsewhere. Both are modulated by
  the same `r(t)` through unrelated geometry. So the acceptor decay is a second,
  **independent** projection of `p(r, t)` — which is what stage 0 showed the
  parameters need — and **nothing is separable**: the decay is not a quenching
  decay times a FRET decay, and neither rate may be averaged over the volume
  before the other is applied. Both must be evaluated along the *same*
  trajectory. That makes `fret_rate_trace` (acceptor averaged over its AV) the
  fast-acceptor approximation and `fret_rate_pair_trace` the general case; the
  known ~7 % gap between them on 148l E15→E90 is this effect, not a numerical
  detail.
* **The acceptor decay carries a rise**, its rate being the donor's *total*
  depopulation — radiative plus quenching plus transfer — seen on a channel free
  of the donor's own detection artefacts, with the amplitude ratio carrying the
  transfer efficiency.
* **It costs parameters as well as adding data.** The acceptor is a different
  dye at a different site with its own local quenchers, so it does **not** share
  the donor's θ. Whether it is a net gain for identifiability is measurable by
  the same Fisher construction as stage 0, and that is the stage-2 work; the gate
  is that it must reduce the joint condition number *after* its own parameters
  are counted.
* **A quenched sub-population is not the FRET-active sub-population**, so
  stage 1's free-fraction nuisance and the transfer efficiency interact — the
  free fraction cannot be assumed FRET-inactive.
* **No experimental data exists at the moment (owner), so "fit real decays" is
  blocked rather than pending.** Until data exists every result in PRD-111 is a
  statement about the model, not about a dye on a protein, and model-perturbation
  (stage 1) is the only validation available. Recorded as a standing caveat on
  the PRD, not a to-do.
* **`fret_rate_pair_trace`'s mismatch error named the wrong lever.** It asked for
  the same `t_max` and `t_step`, which are already shared; the number of **walks**
  is what has to agree. Unequal lengths are the normal case — the two dyes are
  quenched differently because they sit at different sites — so the message now
  says so, and says why truncating is refused (a trajectory concatenates one walk
  per excitation, so an arbitrary cut splits a walk). Matching the counts is
  QuEst-side work, recorded in `../quest/okf/log.md`.

## 2026-08-18 (PRD-111, quenching identifiability across sites)
* **PRD-111 stage 0 measured; the gate passes** (`okf/prds/prd-111.md`,
  `okf/validation/quenching_multisite.md`,
  `benchmark/kq_sensitivity_analysis.py`). PRD-110 stopped because θ is not
  recoverable from one decay; the successor question is *what has to be
  measured*, and the cheapest candidate is more sites. θ is **global** (one dye,
  one chemistry, one mobility) while the geometry is **per site**, so N decays
  share one θ and the joint Fisher information is the **sum** — which beats its
  terms only if the blind directions rotate with the geometry. Six T4L chain-A
  sites spanning quencher-in-contact (A124, TRP126 at 3.5 Å) to nearly-bare
  (A53): joint condition number **1.8 × 10⁵** against 1.96 × 10⁸ for the best
  single site. The rotation is the whole effect — joint `λ_min` **152.7** against
  **0.0456** for the sum of the per-site minima, a 3 350× gain. A joint fit
  confirms it: `free_diffusion` −52 % → **+2.5 %**, `kQ_scale` −69 % → **+4.3 %**,
  in **97 evaluations against 676**.
* **`contact_distance` and `slow_factor` are one parameter, not two.** Fisher
  predicts ±4.1 % and ±0.2 %; the fit misses by −33.4 % and −2.8 % in exactly
  compensating directions — at **reduced χ² = 1.016**. Not a failed fit but a
  second, exactly equivalent minimum, reproducible across resolutions (4.33 Å /
  0.957 at 1.5 Å, 4.51 Å / 0.960 at 2.0 Å). The two build the *same mobility
  field*: log D correlates at 0.927, equilibrium occupancies overlap **98.7 %**.
  `slow_factor` applies once per contacting atom and `contact_distance` sets the
  count, so `slow_factor^n(contact_distance)` is invariant along a valley. **No
  number of sites separates them.** `JᵀJ` error bars must not be quoted for
  either — a local linearisation cannot see a second minimum, and here it is
  confidently wrong about its *best*-determined parameter. The model needs
  re-parameterising on the quantity the data determines, not more data.
* **Repeated end to end at 2.0 Å**: joint condition 3.04 × 10⁵, rotation gain
  5 623×, largest principal angle 88.4°, same three parameters recovered and the
  same pair degenerate. Not a grid artefact.
* **`compute_av` silently discarded a declared grid resolution, and it had been
  discarding one all along.** It *writes* `simulation_grid_resolution` into
  `source_info` from its `disc_step` argument rather than reading it back, so
  both quenching benchmarks passed the resolution into a field that was
  overwritten: their `--resolution` flag was inert and **every run, PRD-110's
  included, was at the 1.5 Å default while the validation pages recorded 2.5 Å.**
  The numbers were always internally consistent — only the label was wrong, and
  both pages now carry the correction. `compute_av` raises on disagreement
  instead of discarding, the benchmarks assert the spacing they got, and
  `test/fret/test_av_resolution.py` (7 tests) pins it. A resolution is the one AV
  parameter whose being wrong is invisible in the result.
* **The stated hypothesis was wrong, and recorded as such.** Predicted before
  measuring: *diffusion only matters when the dye has to travel*. The ordering is
  the opposite — sensitivity to `free_diffusion` tracks **quenched fraction**,
  and correlates with median dye–quencher distance at **−0.56**. Transport
  becomes rate-limiting when quenching is fast enough to deplete. **To measure
  dye mobility, label where the dye is strongly quenched.**
* **A diagnostic that was pure arithmetic, caught and fixed.** The first
  subspace-angle comparison used rank-3 blind subspaces in a 5-parameter space;
  two subspaces of dimension `r` in `Rⁿ` must intersect once `2r > n`, so it
  reported 0.00° for all 15 site pairs regardless of geometry. Rank is now capped
  at `n/2` and the decisive number is the eigenvalue gain, which no choice of
  rank can flatter.
* **PRD-111 stage 1 measured; that gate passes too.** Two perturbations of the
  model, both realistic. A **per-site free-dye fraction** of 2–7 % (incomplete
  labelling, a stuck rotamer; it competes with `slow_factor` because both slow
  the decay) and a **5 % linker-length error** (data at 21 Å fitted at 20 Å,
  which grows the AV from 1 158 to 1 443 voxels).
* **Unmodelled free dye wipes out the multi-site gain**: `free_diffusion`
  −60.5 %, `kQ_scale` −29.5 % — back to PRD-110's single-decay errors. Six sites
  sharpen a correct model; they do not protect a wrong one. **But it is loudly
  detectable**: reduced χ² 1.683 on 379 dof is a **9.4 σ** rejection, so the
  failure mode is a fit that visibly fails, not a wrong number quietly
  published.
* **Modelling it works, and `free_diffusion` pays.** 11 parameters (5 global +
  6 per-site) restore χ² to 1.021 (0.3 σ) and bring `kQ_scale`, `rC` and
  `slow_factor` back to a few percent, at 2.3× the forward solves — but
  `free_diffusion` settles at −13 %, five times worse than stage 0, because a
  free fraction and a fast diffusion both raise the long-time intensity.
  **Stage 0's ±3 % on `free_diffusion` was a perfect-model number.** The fitted
  fractions themselves are badly recovered (0.040 → 0.009) and must never be
  quoted as labelling efficiencies.
* **The geometry error is benign — the opposite of the expected answer.** The
  5 % linker error was included expecting it to be the damaging one, since a
  joint fit converts per-site model error into global bias. With the nuisance
  fitted, `kQ_scale` returns to **−1.6 %** and `rC` to −2.1 % (both *better* than
  under the correct geometry) and χ² rises only to 1.110 — 1.5 σ, not a
  rejection. **Worry about the sample, not the structure.** Untested and
  different in kind: a *systematic* geometry error shared by every site, which a
  per-site nuisance cannot absorb.
* **All of stage 0 is synthetic data from the model being fitted**, so it
  measures the model's internal geometry under a perfect model and nothing about
  the world. Stage 1 (a per-site nuisance parameter, then real decays) is where
  that assumption is tested, and it is the assumption most likely to be false.

## 2026-08-17 (cgdye, PRD-107)
* **PRD-107 written and stage 0 landed** (`okf/prds/prd-107.md`; plan grilled
  with the user, decisions recorded in the PRD). Hygiene found three real
  defects hidden by stale-path skips: `scripts/rotamer_hgbp1_site481.py`
  shipped with a SyntaxError; `scripts/rrt_hgbp1_site481.py` called a fictional
  `IMP.kinematics.Sampler` API (rewritten on a torsion-space RRT in
  `sampling/rrt_imp.py`); `rotamer/io.load_protein_frames` never called
  `IMP.rmf.load_frame`, so every RMF frame returned the same coordinates.
  De-skipped `test_reference_attachments` (all 33 libraries, placement equals
  the reference transform to 1e-6), `test_rotamer_sampling`, the RRT script
  tests, the integration test (builds its system.cif itself), the RMF test.
  click made optional for the library (`utils.import_click`), the fake
  "Langevin" walk renamed `dof_walk`/`sample-dof-walk` (no shims — hard-rename
  policy), `RotamerFRET` no longer writes a `log` file into cwd. Cheap suite
  `test/cgdye test/fret`: 103 passed/14 skipped → **190 passed/5 skipped**.
* **PRD-107 stage 1:** one CHARMM36 table and one LJ kernel
  (`topology/dye.py::lj_energy`) behind all three cgdye scorers, `LJ_PARAMETERS`
  gone, scalar `evaluate` routed through `evaluate_batch`, unused
  `sampling/rrt.py` deleted (`rrt_imp.py` → `rrt.py`), a latent `_site_mask`
  list bug fixed. `test_lj_single_source.py`. Cheap suite 193 passed / 5 skipped.
* **PRD-107 stage 2 / PRD-106 cgdye half:** `pyext/src/fret/strip.py` is the
  one strip engine (grammar → `StripSelection`; PDB lines, obstacle arrays,
  IMP hierarchies; non-destructive by default); `fret/av.py` imports it,
  cgdye `strip_sidechain_at_site` delegates with its own keep-set (strips CB
  — the linker replaces the side chain from CA; AV keeps CB). PRD-106 status
  → in-progress with a record of what remains (full grammar, diagnostics,
  authored-mask fixtures). Cheap suite 199 passed / 5 skipped.
* **PRD-107 stage 3 — FRETpredict parity.** Scratch env with the local
  FRETpredict checkout: Hsp90 agrees to 5e-8 (the old "1.1e-3 gap" was a
  pin recorded at T=300 vs FRETpredict's 293). pp11 (cutoff10 libraries)
  found a real bug: every library name resolved to the cutoff-30 RMF
  template, so `cutoff10/20` requests silently scored 33 rotamers instead of
  711/123 (Es off by 0.18). Fixed in `rotamer/io.py` (cutoff substituted; the
  FRETpredict DCD sets in `data/rotamer_library` are canonical, RMF templates
  a cutoff-30 fallback); `load_protein_frames` reads multi-MODEL PDBs. Fixtures
  gzip'd under `test/cgdye/rotamer/data/`, pins with provenance in
  `test/references/cgdye_fretpredict_pins.json`, `test_fretpredict_pins.py`.
  `test/cgdye/rotamer` zero skips; cheap suite 209 passed / 0 skipped.
* **PRD-107 stage 4 — physics invariants.** `test_physics_invariants.py`
  (κ², R0, single-rotamer regimes, Boltzmann/mean-field/kinetic weights,
  master-equation limits, Kronecker pair). Found and fixed two physics bugs:
  `compute_exact_efficiency` solved the untransposed system (wrong for any
  non-symmetric transition matrix; fast-exchange 0.878 vs 0.764 correct), and
  the linker Metropolis sampler scored all pairs without 1-2/1-3/1-4
  exclusions while collapsing same-named MOL2 atoms into one site (self-pairs
  at r = 0, ~1e6–1e31 per frame). `dye_internal_system` fixes both; sampler
  pins recorded. Cheap suite 226 passed / 0 skipped.
* **PRD-108 stage 3 — Langevin/Brownian dye samplers; PRD-108 implemented.**
  `LangevinDyeSampler` (md: MD + Langevin thermostat; bd: Brownian dynamics,
  dt ≤ 0.5 fs with the stiff bonds), `dye sample-langevin`, hGBP1 script;
  thermodynamic pins (equipartition, MSD = 6Dt, torsion Boltzmann histogram,
  md ⇄ bd, attached dye at 481). Found and fixed a physics bug: torsion types
  in CHARMM convention were passed to `IMP.core.Cosine` (opposite sign) — planar
  π torsions were minimal at 90°, linkers eclipsed — in `dye_restraints.py` and
  `sim/runner.py` (`topology.dye.torsion_cosine`). Cheap suite 347 passed.
* **PRD-108 stage 2 — AV↔rotamer table.** `dye rotamer compare-av` /
  `cgdye/rotamer/compare_av.py` write `okf/validation/av_vs_rotamer.md` and
  drift pins; medium test with loose sanity bounds. Result: rotamer
  ensembles differ from AVs by 3–12 Å in ⟨R_DA⟩ and fit the 99 T4L distances
  worse (Σχ² 257 vs 155); FRETpredict's soft screening leaves most weight on
  rotamers 1.3–2.5 Å from the protein. C++ `R1` in AVNetworkRestraint is not
  motivated by the data; AV-parameter calibration is the follow-on. Also:
  the `dye` CLI's `main()` never existed (now defined; `rotamer` group attached).
* **PRD-108 stage 1 — fps.json `R1`.** Schema type + fields, C++ `IMP_WARN`
  on non-AV `simulation_type` (test with an R1 position in the T4L score set),
  `fps_positions_for_docking`, `write_rotamer_fps` /
  `distances_from_ensembles` / position payloads; JSON schema regenerated.
* **PRD-108 stage 0 — `RotamerEnsemble`.** `cgdye/rotamer/ensemble.py`
  (an `AccessibleVolume` with per-rotamer centre + dipole + weight, all atoms
  kept; fps `R1`), pair kernels `fret_pair_geometry/efficiencies/distribution`
  in `fret/distance.py`; `RotamerFRET._frame_fret` refactored on top with the
  FRETpredict pins unchanged. Tests + exports. PRD-108 page written.
* **PRD-107 stage 5 — flat `IMP.bff` API, docs; PRD-107 implemented.**
  `pyext/src/api.py` + a `%pythoncode` hook in `pyext/swig.i-in` expose 82
  cgdye/fret names lazily as `IMP.bff.<Name>` (the SWIG module never included
  `pyext/src/__init__.py` — dead file, deleted). Hard renames throughout
  (`System`→`DyeForceFieldSystem`, `calculate_r0`→`forster_radius_from_spectra`
  in `fret/forster.py`, `kappa2_from_vectors`→`kappa2_from_dipoles` in
  `fret/kappa2.py`, ... — full list in the PRD), no aliases;
  `test/test_public_api_names.py` asserts docstring + naming family per
  export. Manual notebook `structure_cgdye.ipynb`, README paragraph,
  `okf/cgdye.md`; notes appended to chisurf PRD-47/93/96. Cheap suite 310
  passed / 0 skipped.
* Note for the build: `pyext/src/Files.cmake` and `test/Files.cmake` are
  regenerated by the IMP configure step (`__init__.py` excluded by design);
  new/renamed Python files need a reconfigure (`cmake -S ../imp -B $B`) to be
  symlinked into the build tree, and stale dangling symlinks there must be
  pruned by hand.

## 2026-08-17
* **PRD-105 pass 7:** structure-of-arrays cloud shared by mean and
  quadrature (serial per-AV compute 2.02 → 1.84 ms/frame, bit-identical),
  scratch/capacity reuse, pair tasks queued in start-ready order.
* **Memory pass, PRD-105 code:** `leaks --atExit` over 14 restraint
  lifecycles in all modes (incl. async, pickle, forced recompute, coarse
  mode): no leaks in libimp_bff (only SWIG's one-time registration blocks).
  Guard Malloc over the full suite, the benchmark (all modes) and the
  lifecycle script: clean. Two ownership hardenings: an unfinished
  `evaluate_async()` job is joined by its own destructor and `job_` is the
  last member of the restraint (destroyed first, before the pool and the
  AVs it uses); the registry's coordinate snapshot is shared with its maps
  (`shared_ptr`) so a map handed to Python cannot dangle if it outlives the
  registry. ASan via `DYLD_INSERT_LIBRARIES` into the conda Python spins in
  `AsanInitFromRtl` on this macOS — unusable here; Guard Malloc +
  `MallocCheckHeap` + `leaks` are the tools that work.
* **Placement (user, 2026-08-17): "split at data level — bff owns structure,
  tttrlib owns data/algos."** The decay routines in imp.bff are the
  deprecated 2.25 copies (PRD-93: they leave in the next release; tttrlib is
  canonical and already carries the same `fconv_per_cs` bounds fix in its
  wrapper). The bff clamp is a stopgap in code that is going away, not an
  investment.
* **Memory: `decay_fconv_per_cs` wrote one past the fit array** — its `stop`
  and `conv_stop` are inclusive (the sibling routines' `stop` is exclusive)
  and the Python wrapper maps `stop=-1` to `n_fit`, so `fit[n_fit]` was
  written and the tail added there: heap corruption that surfaced as random
  aborts later in the same process (Guard Malloc pinned it to
  `test_DecayRoutines.test_fconv_per_cs`). Clamped to the last channel;
  results unchanged. Full suite now clean under Guard Malloc.
* **PRD-105: compilation/locality/arrangement pass.** PGO and branch hints:
  no gain; blocked layout not pursued (working set is cache-resident).
  `AVNetworkRestraint::evaluate_async()`/`wait_score()` added: the pool run
  proceeds on its own thread after the serial Model-reading phase, so the
  caller can load the next frame / evaluate other restraints meanwhile;
  bit-identical, tested; ~1 % on T4L (frame load is 0.04 ms), meant for
  callers with heavier per-frame work.
* **PRD-105 SIMD:** NEON row-form relaxation for interior tiles (~3 %),
  run-based vectorised cloud/penalty passes (noise). Serial per-AV sum ~2.0
  ms/frame; the kernel is bookkeeping/cache-bound, not arithmetic-bound.
* **PRD-105 pass 6 (search kernel):** obstacles encoded in the cost array
  (`BLOCKED_COST`), reused bucket/queue scratch, per-shape candidate list with
  shell/interior classification (occupancy-only tiles skip the sphere test),
  inline locations, `set_origin_fast` wired. Serial per-AV compute 4.1 → 2.25
  ms/frame; threaded frame ~0.7 ms best-of-8 on a loaded box, now with the
  correct 26-neighbour metric. All results bit-identical.
* **PRD-105 correction + memory pass:** the pass-3 "26-neighbour" stencil
  was 18 (`sqrt(3.0)` radius with `d² <= r²` excludes d²=3) — that, not
  tunnelling, caused the 15–25 % smaller AVs and the search halving; fixed
  (`sqrt(3.0)+1e-6`, test asserts 26/30 offsets), pins regenerated; the honest
  26-vs-30 change is 0.42 Å rms (max 3 Å). Speed is back to ~1.0 ms/frame
  quiet. `search_mode="euclidean"` (exact DDA visibility, straight linker;
  0.82 vs 0.77 ms — a model option, not a speed one) added; the cross-check
  it enabled found the stencil bug. Memory: `set_path_map_header` now rebuilds
  the location arrays on shape change (in-place origin refresh wrote the new
  size into old arrays); the full-suite 1-in-4 abort traced to an
  out-of-bounds read inside pip LabelLib (`Grid3DExt::excludeConcentricSpheres`,
  Guard Malloc) — `medium_test_av.py` skips the LabelLib backend unless
  `IMP_BFF_TEST_LABELLIB=1`; 0/8 crashes since. AV lattice suite clean under
  Guard Malloc / MallocCheckHeap.
* **PRD-105: `search_mode="euclidean"` added** (straight-linker AV: exact
  voxel visibility by DDA, Euclidean cost, no path search; option on
  `AVNetworkRestraint`/`AV`, default stays `dijkstra`). Model, not speed:
  ~40 % of the path-search voxels, 5.6 Å rms distance change, 0.82 vs 0.77
  ms/frame. Approximate visibility chains were tried and rejected (leak or
  over-conservative).
* **PRD-105: hybrid Euclidean/geodesic search and obstacle distance
  transform evaluated and rejected** (details in the PRD): visibility-seeded
  Dijkstra was slower with safe (face-step) visibility and tunnelled with
  26-step visibility; the obstacle D-field is worth ~0.03 ms wall and breaks
  local deltas. Tree left at 97cca77.
* **PRD-105 fourth pass ("all tricks", exact):** reach-only re-centring
  shared extents, one coordinate snapshot per frame, chord-run sphere raster,
  one pipelined pool run per evaluation (rasters → searches → carves →
  pairs with dependency waits; prepare reads the pending classification),
  spinning workers, integer occupancy into penalty/carve, interior flags,
  byte visited, in-place origin. Quiet-machine best 0.69–0.71 ms/frame
  (pre-PRD 37.7, ~53×); run-to-run variance on the shared box (0.7–1.0) now
  exceeds the remaining gains. Grid-spacing sweep recorded (2.5 Å: 0.51 ms at
  1.8 Å rms change; not applied — user parameter).
* **PRD-105 third perf pass (objective < 1 ms/frame): reached — 0.81 ms
  median quiet, 0.94–0.99 loaded (pre-PRD 37.7).** Exact: SoA tile arrays for
  the lattice path (`search_lattice`/`carve_lattice`, lazy `tiles` sync),
  bucket queue without duplicates or per-bucket sort (all tiles in an active
  unit bucket are final), fused source spheres, raw-location cloud, two pool
  tasks per AV + z-slab rasters, `timing_ms_total` diagnostics. **Model
  change:** the lattice path now uses a symmetric 26-neighbour stencil
  (`search_stencil=26`, `30` = historical): the historical offset loops ran
  `-2 ≤ d < 2`, so paths could cross a one-voxel wall towards −x/−y/−z only;
  the leak made T4L AVs 15–25 % larger and shifted means by up to 2.5 Å
  (distances 1.1 Å rms, max 10 Å) — pins regenerated, 30-stencil value kept
  as a regression pin, legacy anchoring unchanged. `quad_k` default 100 → 50
  (0.025 Å max error with moments). `search_grid_factor` (coarse search,
  fine carve) implemented and measured: 2.4 Å rms error for only 20 % speed
  → kept opt-in, default 1. Hardware note: 4P+4E cores make per-task CPU sum
  ~2× the serial CPU; compute wall is ~4–5× over serial.
* **PRD-105 second perf pass**: `set_origin` moved into the threaded compute
  phase, shared rasters refreshed on threads, persistent `internal::ThreadPool`
  (longest-first dynamic scheduling), and — the big one — the lattice path now
  runs `PathMap::find_path_dijkstra_bounded` (exact lazy Dijkstra that stops
  at cost ≥ ll/h; source tile left at default cost as before) on a monotone
  unit-bucket queue (edges ≥ 1 voxel ⇒ relaxations land in later buckets;
  each bucket sorted once by (cost, idx) = a heap's pop order). All results
  bit-identical (scores unchanged to the last digit); T4L default mode
  **4.2 → 1.8 ms/frame** (pre-PRD 37.7 → 21×). Fixed a 1-in-10 flake in
  `test_AccessibleVolume.test_distance_distributions` (MC histogram tolerance
  30000 was ~1.5× the expected 2N; now 60000).
* **PRD-105 perf pass** (same day, one commit): the AVs' compute phases run
  on threads (`AV::resample_prepare/compute/finish` split; restraint
  `set_number_of_threads`, default hardware concurrency; pair sums threaded
  too), neighbours enumerated inline in the search (no per-tile edge vectors,
  index heap, penalty mirror), `fill_sphere` without per-voxel Vector3D,
  cloud from reached tiles, locality-aware skip via per-generation change
  boxes on `AVOccupancyMap` (`get_changed_since`), `PathMap::set_exact_search`
  (lazy Dijkstra; agrees with the historical search on every reached tile
  except the source tile, not faster, off by default). T4L default mode
  18.9 → **4.2 ms/frame**, repeats 0.3 ms; all modes bit-identical to before
  (legacy pins unchanged; new thread-equivalence, local-move and exact-search
  tests). One race found and fixed on the way: with legacy anchoring the quad
  caches were not pre-built, so threaded pair sums wrote them concurrently —
  every AV now refreshes its cache in the compute pass before the pair loop.
* **PRD-105 implemented** (`okf/prds/prd-105.md`, status `implemented`;
  implementation record appended to the PRD). New: `AVOccupancyMap` /
  `AVOccupancyRegistry` (integer occupancy counts on the absolute lattice,
  exact subtract/add deltas, integer-voxel rolls with slab raster, grow-on-
  demand shared extent, `read_window` fill-zero-beyond), `AV::resample_lattice`
  (window `2·floor(ll/h+½)+1` centred on `round(s/h)`, skip when nothing that
  feeds the search changed, `force_full` for the exactness proof),
  `av_distance_quadrature` (≤K lattice blocks with centroid, weight and second
  moments; second-order corrected double sum — max error < 0.005 Å at K=100 on
  T4L @2 Å where centroids alone gave 0.26 Å), `AVNetworkRestraint(...,
  space_fixed, shared_map, distance, quad_k)` + `get_diagnostics_json` +
  `get_used_av`, SWIG kwargs shadow for the overloaded restraint constructor
  (kwargs had silently stopped working when the default ctor was added for
  serialization — `test_AVNetworkRestraint.py` was failing on it). Suite:
  `test/test_av_lattice.py` (6 modes, tier-1 bit-exact vs forced full recompute
  and shared vs private, legacy pins via sha256 in
  `test/references/prd105_legacy_pins.json`, order-independence with xfail on
  mc modes, quad K-curve vs the exact double-sum oracle, diagnostics,
  occupancy-map unit tests). Benchmark `benchmark/benchmark_av_screening.py`:
  T4L default 18.9 ms/frame vs corrected legacy 42.4 (pre-PRD code measured
  37.7 today), repeats 0.8 ms. **Two latent defects found in the legacy path
  and fixed for both paths:** `PathMap::update_tiles` had
  `edge_computed.resize(false, nvox)` (arguments swapped → flags resized to
  zero → `get_edges` read stale freed bits → edges never recomputed after the
  first evaluation, so a moved structure kept its old connectivity: scores
  were evaluation-order dependent beyond MC noise, and the old baseline was
  artificially cheap); `AV::resample` carved tiles with `density *= 0` and
  never restored them, so a moving structure's AVs shrank frame after frame
  (`inf` scores late in the T4L trajectory). Legacy single-evaluation results
  are byte-identical to the old code (verified against a stash-build of
  11cb3e7); repeated-evaluation legacy pins are from the fixed code. Existing
  pins in `test_AccessibleVolume.py` / `test_AVNetworkRestraint.py`
  regenerated once to lattice values (mean position shift 0.18 Å, restraint
  score 11.92 → 13.08 deterministic). Deviations from the locked text
  (recorded in the PRD): torus phase indexing replaced by an in-place shift
  of the retained block (same work bound, plain indexing everywhere); quad
  carries second moments. Pre-existing, untouched: `test_DecayConvolution.py`
  / `test_DecayScore.py` fail with the same overloaded-ctor kwargs TypeError.
* **PRD-105 design locked via grilling (15 decisions) — scope rewritten**
  (`okf/prds/prd-105.md`, status `specified`). The original fast-update-only
  scope is superseded: measurement showed all 159 T4L beads move every frame
  (incremental raster can't pay) and rolls fire ~every frame under any window
  policy (padded windows buy nothing). Locked v1: global absolute lattice as
  exactness reference, default on, legacy `space_fixed=False` opt-out
  (deprecated); minimal linker-defined rolling windows (PBC, integer-voxel
  rolls); one shared occupancy raster per (spacing, extra-radius) class —
  17 rasters → 4 classes on T4L, ~14.7 → ~2 ms/frame; cold Dijkstra
  unchanged (warm start = v2); deterministic lattice-quadrature distances
  K=100 replace MC (~14.5 → ~0.6 ms; MC-class accuracy 0.16 Å vs oracle;
  moments rejected at 2.2 Å bias; FFT oracle test-only; numpy timings
  implementation-unfair — C++ A/B required). Distances live in this PRD,
  independent of PRD-94. Two-tier exactness contract (tier 1 bit-exact maps,
  tier 2 deterministic distances); 6 flag modes (`space_fixed`, `shared_map`
  requires it, `distance=mc|quad`), all parameterized in CI with xfail on
  MC order-dependence (pre-existing defect, fixed by replacement); one merge,
  three staged commits; done = suite + tables, no speed gate. All grilling
  A/Bs ran on the real T4L restraint (`chi2_C1_33p`) in the arm64 env.
* **PRD-105 ideas extended — "semi space-fixed" grid** (idea 3 in the
  non-committed section): PBC (torus-indexed) window whose origin is only ever
  relocated by *integer multiples of the spacing* on a fixed global lattice.
  A roll is an O(1) index remap (cached integer occupancy counts stay valid —
  voxel centres stay on the lattice); between rolls the grid is exactly
  space-fixed so the bit-exact subtract/add deltas apply; per roll only the
  newly exposed face slab is rasterised. Wrap-around phantoms are provably
  masked if the window is padded by the reach (they land ≈2·ll from the source,
  inside the `fill_sphere` blocked-beyond-linker region). Conditions recorded:
  no torus-wrapping path edges, hysteresis on the roll trigger, fallback full
  recompute must quantise the origin with the same rule (one shared lattice),
  phase-aware coordinate/export variants. User asked for this variant
  ("pbc + relocation of origin, semi space fixed"); still collect-only, not
  implemented.

## 2026-08-16
* **PRD-106 authored in this bundle** (`okf/prds/prd-106.md`, feature PRD:
  strip mech — a selection-driven mechanism to strip side chains / whole
  residues / atom masks from the obstacle set). Grounding found in code: the
  fps.json `strip_mask` field is declared in `fret/fps_schema.py:163` and
  authored in the shipped T4L + TG2/flex examples (grammar `chain X and resid N
  and not name ...`, `or`, `resname`, parentheses) but has **zero consumers** —
  the AV build substitutes `allowed_sphere_radius` source clearance for
  stripping (`fret/av.py:124-132`); a whole-residue stripper
  `_strip_residue_atoms` (`fret/av.py:557`) is dead code; cgdye's
  `strip_sidechain_at_site` (`cgdye/labeling/attachment.py:73`) hard-codes a
  keep-set (`N CA C O OXT`, strips CB) that disagrees with the AV masks' keep
  (`CA CB C N O`). Scopes one engine in `fret` (grammar-faithful PyMOL subset,
  non-destructive, two outputs — stripped hierarchy + filtered obstacle array),
  live `strip_mask` in the AV build, cgdye delegation with behaviour pinned by
  the existing hGBP1 tests, dead-code removal. Cross-repo consumer: chimol's
  `_stripped_pdb_for` seam (T-20260815-01) hand-derives the same strip today.
* **Ideas collected into PRD-105 as a non-committed design section** (space-fixed
  grid, shared/integrative AV map, better-than-dense-cube data structures). User:
  "ideally make a space fixed grid ... points appear and disappear from grid",
  "move from simple grid to a better data struct"; collect, do not implement.
  Grounding found in code: grid is re-anchored to the source on every `resample`
  (`src/AV.cpp:161-165`, `src/PathMapHeader.cpp:55-74`) so voxel indices are not
  spatially stable across frames; every AV rasterises *all* root leaves into its
  own private cube (`src/AV.cpp:149`) — 17× duplicated on T4L. Ideas: world-
  anchored grid (envelope), rebase policy with integer-spacing-only exactness
  caveat, one shared occupancy grid per frame for all AVs, blocked/chunked grid +
  per-atom voxel footprint cache + hybrid sparse-hash, octree as stretch.
* **PRD-105 authored in this bundle** (`okf/prds/prd-105.md`, feature PRD:
  AV fast update — reuse the previous computation when local geometry barely
  changed). Every `AV::resample()` rebuilds the whole path map (obstacle
  rasterise + full Dijkstra; measured 96 % of an AV build on 148L at 0.5 Å
  grid, `pyext/src/av/compute.py`). The PRD scopes a fast-update path on
  `AV`/`PathMap`: short-circuit when nothing within the label's sphere of
  influence moved, locally re-rasterise + incrementally relabel only the
  changed neighbourhood otherwise, full fallback when the change is too
  large — under an exactness contract (fast == full recompute, or the fast
  path does not run). Consumers: trajectory screening
  (`plot_AVScreening.py`, `imp_engine.screen()`), full-AV sampling,
  `dock_minimize` refinement cycles. Rule (user, 2026-08-16): bff feature
  PRDs live in `okf/prds/` of this bundle, not in the shared chisurf bundle;
  `okf/index.md` scope-rule updated to match.

## 2026-08-11
* **PRD-97 stages 0–3: `IMP.bff.fret` created.** The FRET docking engine,
  the AV backend, six algorithm modules and the single fps.json reader moved
  here from ChiSurf's `fret/core` (details in
  `../chisurf/okf/prds/prd-97.md`). Rules applied during the move (user,
  2026-08-11): **no LabelLib** (IMP.bff's AV is the only backend) and **no
  numba** (the two jitted kernels are vectorised numpy). `fps_schema.py` is
  the authored fps.json definition — flrCIF item names where flrCIF has
  them — and `data/fps_json_schema.json` is derived from it with a drift
  test. `pyext/src/fps.py` deleted (both importers now use
  `IMP.bff.fret.io`). Fixed while landing: `src/AV.cpp` wrote radius1 into
  all three radii (`set_av_parameter`); the moved AV path had two crashes the
  old LabelLib fallback had hidden (ndarray truth test, argument-less
  `DensityHeader.get_origin()`); AV source clearance now scales with linker
  width (`allowed_sphere_radius >= lw/2 + grid/2` unless the position sets
  its own). `examples/structure/TG2/flex.fps.json` `S1_val_chi2` prefixes
  fixed; `hGBP1.fps.json` `577_577` dangling reference flagged, not fixed.
  Tests: `test/fret/` (19) plus cgdye all green — 93 passed, 14 skipped.
* **Bundle created.** `imp.bff` gets its own OKF bundle for repo-local
  knowledge, scoped deliberately narrow: cross-stack knowledge and all PRDs
  stay in `../chisurf/okf/`, photon-level concerns in `../tttrlib/okf/`, and
  the shared agent board is symlinked from `tttrlib/okf/agent-board.md`.
  Initial contents: `index.md`, `overview.md` (the sibling stack — tttrlib,
  imp.bff, imp, imp-tricks, chisurf, fpsimp, quest, ucfret — and the
  where-knowledge-lives map), and `references/index.md` pointing into the
  shared bundles. `AGENTS.md`/`CLAUDE.md` updated to match.
