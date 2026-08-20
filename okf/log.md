# Update Log

## 2026-08-20 (later) — the keystone, and the three things that were waiting on it

`States`/`AccessibleVolume` (Python, `representation/distance.py`) and
`BasicAV`/`ACV` (C++, `AVModel.h`) described the same thing: a weighted point
cloud with a grid. The Python called the C++ for every computation while
holding its own idea of what the cloud *was*, which is exactly how the two came
to disagree about the axis order of `density` (PRD-113 stage 3a). There is one
now — `States` → `AccessibleVolume` → `ACV`, all C++ — and `BasicAV` is retired
rather than aliased.

What that unblocked, in the order it unblocked it:

* **`Quencher`, `PETParameters`, `ResidueQuenching` and the PET tables** →
  `PETQuenching.h`. The values were in `label.py` and the tables they read were
  in `quenching.py`, so `reference_quenchers` crossed a module boundary one way
  and every consumer crossed a language boundary the other, for a rate that
  takes two partners to define.
* **`InteractionTerm` and its three implementations** → `InteractionTerms.h`,
  as `IMP::Object`s because `total_rate` sums a list of *different* terms.
  Gated on 40 states against 60 atoms: radiative and FRET identical, PET and
  the sum to 3.6e-15.
* **`lifetime_spectrum_from_states`** → `LifetimeSpectrum.h`, which is what
  `IMP_bff.observables.i` had said it was waiting for.
* **`residue_sites` and `atomic_quenching_parameters`** → `PETQuenching.h`,
  on parallel arrays rather than a numpy structured dtype.

Shapes that changed rather than moved, and why:

* **`grid_shape` is derived from the density**, not a constructor argument that
  had to agree with it. A non-cubic grid is refused now instead of half-believed.
* **`mean_position` falls back to the attachment point** for an empty or
  zero-weight cloud. The Python `States` did; the C++ returned the origin. A
  buried site whose volume came back empty is at its attachment atom.
* **The structure door stops storing density as float32** — a storage choice of
  the deleted dataclass, costing seven digits for nothing.
* **One participant signature for every term.** Python's duck typing let each
  take what it needed; C++ takes the same two and reports through `arity` which
  it reads. `PETTerm`'s quenching atoms move into the term, because resolving
  which of a structure's atoms quench is a function of names and parameters and
  does not change per call.
* **NaN is the C++ spelling of `None`** for `quench_radius`, `attenuation_length`
  and the `Dye` numbers. `FRETTerm` tested its donor's lifetime with
  `is None or <= 0.0`, which NaN passes; `not (tau0 > 0.0)` is the test that
  catches all three.

Deleted rather than ported, all of it with no consumer anywhere:
`quencher_atom_indices`, `quencher_centers`, `rate_constants` (a named
pass-through to `total_rate`), the legacy bare-number entry in
`normalize_amino_acid_quenching`, and `BasicAV.save_xyz` — whose test had been
failing since the C++ port without anything noticing, because
`pytest test/` does not collect a `medium_test_*` file. Two more in the same
file called `ACV.from_basic_av` with a signature it never had.

**Where this leaves `pyext/src`:** 30 files / 22,736 lines → 22 / 19,285. The
remaining order is unchanged apart from the keystone being done:

1. **The fps.json layer** (`io/fps.py`) against `internal/FPSReaderWriter.h`,
   which already reads the same format in C++ for `AVNetworkRestraint`. Untyped
   dicts on the Python side; the same treatment `DyeForceFieldSystem` got.
2. `scoring.py`, `io/cif.py`, `io/structure.py`, `representation/rotamer.py`
   (its `from_site` and pair physics; the value itself is now a `States`),
   `representation/av.py`, `cgdye/*`.
3. **Programs, not ports**: `restraints/docking.py` and `cgdye/sim.py`.
   `docking.py`'s only consumer is ChiSurf, through a forwarder
   (`chisurf/plugins/modelling/fret/core/imp_engine.py`) that points at
   `IMP.bff.fret.imp_engine` — **a module path that no longer exists**, so that
   bridge is broken today and needs repointing at `IMP.bff.restraints.docking`
   in that repository.

## 2026-08-20 — emptying `pyext/src`: the rule, and nine modules gone

**The rule, stated by the user and now the repository's:** *what can be C++
must be C++ in `imp.bff`.* Tests, examples and documentation are Python;
`prototypes/` is exempt. The motivation is not performance — it is to minimise
the number of places a value has to change language, because every one of those
is a marshalling convention that can drift.

The corollary for `pyext/src` is that its target is **empty**: kernels and
value types to C++, programs to `bin/`, and only what genuinely cannot be
either into `%pythoncode` in a `.i` file. It went 30 files / 22,736 lines →
22 files / 19,936 lines in this pass.

What moved, and what each move found:

* **`restraints/greedy_olga.py` → C++.** Only the greedy loop was left; the
  kernels had moved already. `select_informative_pairs` joins them, so the
  per-step chi-squared accumulation no longer crosses SWIG `max_pairs` times.
* **`tools.py` → `DataPaths.h`.** Five path joins over `get_data_path("cgdye")`.
  Landed first as `%pythoncode` on the reasoning that `pathlib` reads better;
  the rule above says otherwise, and it was redone in C++ the same day.
* **`representation/compare.py` → `bin/imp_bff`.** A validation tool, and
  **both** of its entry points were dead: `rotamer compare-av` imported a `main`
  that was never written, and `av-vs-rotamer` called nine names it never
  imported. Neither has a `test_*.py`, and the module's own test is a
  `medium_test_*`, which plain `pytest test/` does not collect.
* **`analysis.py` → `bin/imp_bff`.** A program end to end; its two path
  defaults pointed at the *package* directory, which an installed module may
  not write to.
* **`observables.py` → C++ and `.i`.** Two names were already C++ and the module
  re-exported them. `compute_distance_distributions` was filed under the
  output contract but wrote a CSV, which that contract explicitly excludes; its
  only caller anywhere is `restraints/docking.py`, where it now lives.
  `test_no_convolution_anywhere_in_the_observables_package` globbed
  `pyext/src/observables/*.py`, a directory that stopped existing when the
  package became one module — it had been scanning nothing and passing.
* **`dye.py` → `DyeLibrary.h`.** `Dye` and `Spectrum` as C++ values, the
  `_bff_dye` CIF read through the vendored `ihm` C reader, the mtime cache with
  it. Gated against the Python over all 38 bundled dyes: identical spectra,
  worst relative overlap difference 4.3e-15, worst R0 difference 4.4e-15 nm.
  An unknown number is NaN now, and `FRETTerm`'s `is None or <= 0.0` test did
  not catch it — `not (tau0 > 0.0)` does.
* **`restraints/simple_av_network.py`, `restraints/direct_labeling.py` →
  `LabelingRestraints.h`.** `BasicAV` was already C++, so the network restraint
  holds the volumes rather than proxies for them. `chi2_score`,
  `fret_efficiency` and `distance_from_fret_efficiency` went to `AVDistance.h`
  with them.

Three guards in the test suite had to change, because each turned a module
*leaving* `pyext/src` into a test edit: `len(DOMAINS) >= 12`, an exact set of
subpackage directories, and `len(modules) > 25`. They now check that what is
there is well-formed and that nothing new appears, so the list is free to
shrink to nothing.

**Still to do**, in dependency order — the keystone first:

1. **`States`/`AccessibleVolume` (`representation/distance.py`) against
   `BasicAV` (`AVModel.h`).** These are two descriptions of the same thing: a
   weighted point cloud with an attachment point, one in Python and one in C++.
   Nothing above them can move until they are one. `States` adds
   `orientations` and a free-form `params` dict; `BasicAV` adds the density
   grid and the three distances.
2. **The interaction terms** (`photophysics.py`) — they consume `States`, so
   they follow it. `rate_constants` and `lifetime_spectrum_from_states` are
   parked in `IMP_bff.observables.i` waiting on exactly this.
3. **The fps.json layer** (`io/fps.py`) against `internal/FPSReaderWriter.h`,
   which already reads the same format in C++ for `AVNetworkRestraint`. The
   Python side is untyped dicts; the same treatment `DyeForceFieldSystem` got.
4. **The quenching tables** (`quenching.py`) — `QUENCHER_ATOMS`,
   `PET_QUENCHING_REFERENCE` and the per-residue normalisation, then
   `Quencher`/`PETParameters` in `label.py` that read them.
5. `scoring.py`, `io/cif.py`, `io/structure.py`, `representation/rotamer.py`,
   `representation/av.py`, `cgdye/*`.
6. **Programs, not ports**: `restraints/docking.py` and `cgdye/sim.py` are
   driver code. `docking.py`'s only consumer is ChiSurf, through a forwarder
   (`chisurf/plugins/modelling/fret/core/imp_engine.py`) that points at
   `IMP.bff.fret.imp_engine` — **a module path that no longer exists**, so that
   bridge is broken today and needs repointing at
   `IMP.bff.restraints.docking` in that repository.

## 2026-08-19 — AV3 was AV1 with two ignored numbers; now it is Olga's AV3

`get_radii()` dropping `radius2` (below) turned out to be the small half of the
problem. **The geometry never read `radius2` or `radius3` at all.**
`src/AV.cpp` carved by `get_radius1()` in every path (`:355`, `:424`, `:527`), so
an `AV3` position produced the AV1 volume of whichever radius happened to be
first — and was therefore order-dependent, where the reference is not.

Measured on T4L 172L site 22, before:

| radii | bff voxels | densities |
|---|---|---|
| `(3.5, 0, 0)` | 14252 | {1} |
| `(3.5, 2.5, 1.0)` | **14252** | {1} |
| `(1.0, 2.5, 3.5)` | **18890** | {1} |

AV3 identical to AV1, and the answer changing when the same three radii are
listed in a different order.

### What the reference actually does

Read from LabelLib's source, not inferred from its behaviour —
`Grid3DExt::excludeConcentricSpheres`, `FlexLabel/src/FlexLabel.cxx:235`
(LabelLib 2af43ac, vendored at `../chisurf/junk/LabelLib`):

```cpp
std::sort(effR.data(), effR.data() + numClashes);
const VectorXf rhos = VectorXf::LinSpaced(numClashes + 1, 0.0f, maxRho);
...
for (iClash...) for (; neighbours[iNei].r <= curR; ++iNei) ref = min(ref, rhos[iClash]);
```

The radii are **sorted**, `rhos` is `{0, 1/3, 2/3, 1}`, and each atom writes
`min(ref, rhos[i])` over the shell between consecutive `atom_vdW + radius[i]`.
So the density is **the fraction of probe radii that fit**, and the sort is why
the result cannot depend on their order. `dyeDensity` builds the AV1 case by
passing a one-element vector to the same routine, so AV1 and AV3 are one
algorithm there, not two.

A consequence worth stating because bff's convention differs: **a zero radius is
a real probe**, one that fits everywhere. `dyeDensityAV3(3.5, 0, 0)` is *not*
`dyeDensityAV1(3.5)` — measured, 21625 voxels against 16370, densities {2/3, 1}.
bff writes AV1 as `(r, 0, 0)` and the fps schema requires `radius2`/`radius3`
positive for an `AV3` position (`io/fps.py:888`), so bff selects the **positive**
radii: one gives the AV1 carve, three give the AV3 carve. That reproduces both
conventions without changing any existing AV1 result.

### The change

* `AV::get_active_radii()` — the positive radii, with the zero-sentinel
  reasoning recorded on it.
* `PathMap::carve_lattice_fractional()` — density `k/n` from one occupancy-count
  array per radius. `n == 1` reproduces `carve_lattice()` exactly, so AV1 is
  untouched.
* `AVLatticeState` carries one occupancy source per dye radius (registry-backed
  and private branches both), and the carve reads them all.
* `resample_legacy` gets the same treatment, so both dispatch paths agree.
* **`compute_av` stopped binarising the density.** It did `np.where(d > 0, 1, 0)`,
  which is a no-op for AV1 and destroys AV3 — it returns the AV1 volume of the
  *smallest* radius and silently discards the other two. Same defect as
  `get_radii()`: an AV3 that looks like it worked. Every existing consumer
  (`compare.py`, `distance.py`, `simple_av_network.py`) passes AV1 radii, where
  the carve already emits exactly 0 or 1, so nothing else moves.

After: `(3.5, 2.5, 1.0)` → 18890 voxels, densities **{1/3, 2/3, 1}**, and the
three permutations of those radii are identical.

`test_av3_matches_labellib_rule` pins the four properties the rule implies, each
of which failed before: density is the mean of the per-radius indicators; the
level set at `k/n` equals the AV1 volume of the k-th largest radius;
`AV3(r,r,r) == AV1(r)`; and the order of the radii does not matter.

### The search metric: the default is now LabelLib's, with a speed option

bff's AV and LabelLib's differed by ~15 % **before any dye radius is applied**, so
this was the core search, not the carve. With **no obstacles at all** and a 20 Å
linker, where the answer must be a sphere and needs no reference to check:

| | voxels | vs LabelLib |
|---|---|---|
| analytic sphere | 33510 | — |
| LabelLib | 30688 | 1.0000 |
| bff, stencil 26 | 26146 | 0.8520 |
| **bff, stencil 74 (now default)** | **30682** | **0.9998** |

**Cause.** LabelLib's `essentialNeighbours()` (`FlexLabel.cxx:149`) keeps the
shells with squared offset length {1,2,3,5,6} — 74 edges — noting that a shorter
list gives "isopath surfaces that are cubic instead of spherical". bff had
{1,2,3}. It did not lose volume at the rim; it lost it in *every* shell,
under-reaching along the diagonals. `set_search_stencil(74)` sets the neighbour
radius to √6, which also admits the six squared-length-4 axis jumps — exactly
redundant, since their cost equals two unit steps — so the path lengths are
LabelLib's.

**Three things now exist that did not.**

1. **The stencil is an option in the Python API.** It was reachable only on the
   decorator, so `compute_av` / `compute_av_from_structure` could not select it
   at all; both now take `search_stencil`, and it is recorded in the returned
   `params`.
2. **74 is the default.** The reference metric is what you get unless you ask
   otherwise.
3. **26 is the speed option, and it is a real one**: on T4L site 22 at 1.0 Å a
   re-resample is 1.7 ms against 3.5 ms.

**Compensation, so the speed option is not also a volume option.** A coarse
stencil overestimates path length, so its volume is that of a *shorter* linker —
and the bias is a property of the stencil, not the structure. Measured
obstacle-free over L = 12–25 Å and dye radii 1.0–3.5 Å, the linker-length scale
that makes 26 reproduce 74 is **1.0551 ± 0.0021** (range 1.0526–1.0592); stencil
30 measures the same, being 26's asymmetric variant. Applying it:

| config | voxels (T4L 22) | vs reference | time |
|---|---|---|---|
| 74, reference | 19155 | 1.000 | 3.5 ms |
| 26, raw | 15996 | 0.835 | 1.7 ms |
| **26, compensated** | **19068** | **0.995** | **2.1 ms** |

So the speed option costs ~0.5 % of the volume rather than ~16 %, at 60 % of the
runtime. `set_compensate_stencil()` is **opt-in**: asking for a stencil gives
that stencil's own answer unless compensation is requested. That is not
fastidiousness — the first version defaulted it on and broke
`test_search_stencil_and_grid_factor_options`, which pins the *historical*
stencil-30 metric and is entitled to get it.

**What moved, and it is not small.** Changing the default changes every AV.
Pinned values updated with the reason recorded beside each: mean AV position by
up to **2.6 Å**, ⟨R_DA⟩ from 54.56 to **55.86 Å**, R_E from 53.78 to **54.82 Å**,
FRET efficiency from 0.450 to **0.421**, and the R_DA distribution regenerated
(as the mean of 15 × 10k runs, so the reference is not itself one noisy draw).
Anything that compared bff numbers to previously recorded ones must be re-read.

**The cost of the default, stated plainly.** 74's longest jump is √6 ≈ 2.45
voxels, so it **crosses gaps thinner than that**, which 26 (√3 ≈ 1.73) cannot.
`test_empty_av_gives_nan` was a concrete instance: a site sealed at 1.5 Å spacing
with `linker_width` 0.5 returned 27 voxels 16–20 Å away. LabelLib has the same
property — it is what the reference metric *is* — so a grid coarse enough to
leave an obstacle layer under ~2.5 voxels is unsafe with either. That test now
makes its empty AV with a linker too short to reach anywhere, which is
stencil-independent, and the tunnelling is documented on
`AV::get_search_stencil` rather than hidden behind a passing test.

### Still open on the AV, and not chased

A residual gap remains **with obstacles**: at T4L 22 the default now gives 22374
voxels against LabelLib's 21348 — 4.8 % *larger*, where obstacle-free the two
agree to 0.02 %. The sign and the dependence on obstacles point at the linker
contour-length cutoff (bff blocks beyond a Euclidean sphere of the linker length;
LabelLib thresholds on the Dijkstra **path length**, `setAboveThreshold(linkerLength, -4)`),
which would keep voxels LabelLib drops. Not investigated further.

**Consequence to know about:** `prototypes/quench_pinn` builds its AVs with
LabelLib directly, so it is unaffected by any of this; but its volumes and bff's
are now within ~5 % rather than ~15 %.

### Two things this leaves open

**Duplication.** LabelLib already implements this, and bff now implements it
again — PRD-112 stage 1 removed LabelLib as a backend on 2026-08-11, so the
duplication is that decision's cost, not something introduced here. It is worth
revisiting deliberately rather than by drift: the alternative is that AV1/AV3
geometry is *one* implementation somewhere and bff consumes it.

**AV3 is still untested in practice.** Nothing in either repo *uses* AV3 — no
`fps.json` in the tree sets `simulation_type: "AV3"`, and no LabelLib AV3 call
existed anywhere. That is why both defects survived. One fix plus one unit test
does not make the path exercised; a real AV3 position in an example or a
regression fixture is what would.

Verification: `ninja IMP.bff` clean; full suite **902 passed, 3 xfailed, 30
subtests** (`medium_test_av.py` excluded — it has 3 failures that pre-date this,
`BasicAV.save_xyz` missing from the source entirely and `ACV.from_basic_av`'s
signature drifted, neither reachable from this change).

## 2026-08-19 — `AV::get_radii()` never returned `radius2`

`include/AV.h:316` built its result from `radius3` twice:

```cpp
return IMP::algebra::Vector3D({get_radius1(), get_radius3(), get_radius3()});
```

so **`radius2` was never returned**, contradicting the method's own docstring two
lines above it. `get_parameter()` at `:145` reads all three correctly, which is
what made the discrepancy visible at all.

**Why it survived.** The only consumer is
`pyext/src/restraints/network.py:118`, `r_mean = max(av.get_radii())`, so the bug
is invisible unless `radius2` is the strict maximum — which needs an AV3 dye, and
**nothing in the repo uses AV3**. Every current path is AV1, where
`radius2 == radius3 == 0` and the duplication cannot be seen. The existing
assertion in `test_AccessibleVolume.py` pinned exactly that case,
`get_radii() == (3.5, 0, 0)`, so it passed either way.

That is the general shape of it: AV3 is implemented end to end —
`compute_av_from_structure(radii=(r1, r2, r3))`
(`pyext/src/representation/av.py:441`), the fps schema (`pyext/src/io/fps.py:537`,
`radius2`/`radius3` validated positive at `:888`), and the decorator
(`include/AV.h:194`, `src/AV.cpp:804`) — but it is **untested in practice**, so a
defect on that path had nothing to trip over.

**Fixed**, and pinned with a test that fails on the old code rather than one that
merely passes on the new: `test_get_radii_returns_all_three` sets three
**distinct** radii (5.0 / 4.5 / 1.5, and a case with `radius2` largest so the
`max()` consumer is exercised). Verified in both directions by reverting the
header, rebuilding, and watching it fail with
`ACTUAL [5.0, 1.5, 1.5]` against `DESIRED [5.0, 4.5, 1.5]`.

Found while reviewing whether AV3 could give the prototype's flat dye a
three-radius representation — see
`prototypes/quench_pinn/okf/validation/s20-dye-model-review.md`. **AV3 should not
be trusted until this is in**, and it is worth treating the absence of any AV3
test as the real gap: one fix does not make an untested path safe.

Verification: `ninja IMP.bff` clean; `test_AccessibleVolume.py`,
`test_av_lattice.py`, `test_label.py`, `io/test_formats.py` — 58 passed,
3 xfailed, 30 subtests passed.

## 2026-08-19 (the objects move to C++: bff starts becoming C++-carried)

pmi's shape was a waypoint, not the target. `atom` is 314 lines of Python over
13,137 of C++ and `core` is 495 over 10,305, and they are shims **not because
their kernels are in C++ but because their domain objects are**. bff had 16,404
executable Python lines against 37,000 of C++, and only **1 % of its calls
reached a C++ kernel** — the kernels were done. No amount of further kernel
porting changes the shape while an object is a Python class holding numpy
arrays.

Four objects and three kernels moved, each gated against the Python it replaced
*before* that Python was deleted:

| | removed | speed | gate |
|---|---|---|---|
| `LifetimeSpectrum` | 161 | — | its 22 tests, untouched |
| `BasicAV` + `ACV` | 465 | — | cloud and all three distances bit-identical |
| `DyeDiffusionSimulation` | 171 | `sample_grid` **11.6×** | **bit-identical**, incl. 4 concatenated walks |
| `GridDiffusionSolver` | 292 | — | adjoint gradients to 2e-14 relative |
| κ² sampler + 3 kernels | 267 | `kappasq_dwt` **43×** | limits exact; 5e-4 on 200k draws |

Ratio C++ : executable Python went **2.26 → 2.48 : 1**. For calibration: `isd`
is 1.37, `core` 21, `atom` 42.

**The bit-exact gate was only possible because the seed derivation moved
unchanged.** `(base + i * 104729) % (2**31 - 1)` — a large prime stride, so
walks from one base seed do not share a low-order pattern. A different
derivation in C++ would have forced a distributional comparison and a far
weaker gate.

**Two mistakes, both caught by the suite rather than by my own gate**, because
the gate called everything by keyword and the suite does not:

* **the `GridDiffusionSolver` constructor order.** It is `(diffusion_map,
  bounds, density, ...)` and callers pass positionally. Swapping `D` and `k`
  surfaces as a stability error a long way from the cause.
* **`k_max` does not belong in the stability bound.** The rate is carried as
  `exp(-k dt)`, the exact solution over the step, so it constrains nothing; it
  is the `1 - k dt` form that diverges. Including it cut the allowed step by a
  third on a site with `k_max = 96 1/ns`. `diffusion_stability_limit`'s own
  docstring describes the older form, which is what misled me.

And one caught by the gate, which is what gates are for: the distance-ratio
transform fills **zero** outside the sampled range (`np.interp(left=0, right=0)`)
and renormalises *after* interpolation. Clamping instead was a 0.12 error in
the weights while the axis and ⟨κ²⟩ still matched to 2e-16 — the shape of error
that looks like a working port if only the easy outputs are checked.

`test/test_cpp_objects.py` writes all of that down as invariants, because a
gate against deleted code cannot be re-run.

**Two corrections to earlier claims in this log's neighbourhood**, both from
measuring the wrong thing:

* **cgdye is pinned by one import, not seven.** Of the files binding
  `IMP.bff.cgdye`, 31 are in imp-tricks' `build/lib/` — a stale copy of imp.bff
  itself — 2 are examples that already fail (`IMP.bff.cgdye.rotamer` was
  removed by an earlier migration), and the single live file binds
  `cgdye.cli.dye` plus `cgdye.rotamer.cli`, the second also already broken. The
  earlier count measured *paths that resolve*, not *live consumers*.
* **counting numpy calls is a poor way to pick a port target.** Of the four
  functions it ranked highest, `compute_rotamer_score` is string-keyed atom
  selection around an inner loop that is already C++, `kappa_distance` was
  already a wrapper, and `convolve_distance_with_k2_ratio` already called
  `outer_product_histogram`. Counting *arithmetic* rather than array operations
  says so immediately.

**What is left, honestly.** The decision-free numerics are exhausted. Of the
remaining ~15,700 executable lines, the arithmetic-heavy remainder is
concentrated in `cgdye` (3,542 lines, 141 `IMP.atom`/`core` calls and zero
numpy — C++-native work) and `io/cif` (2,855, and `ihm` has no C++ equivalent
in the stack, so that is *writing* a CIF parser rather than moving one). Both
need an owner decision. `_chisq_rt_cdf_python` must **stay** Python: it is the
reference the C++ `chi2_right_tail` is gated against.

## 2026-08-19 (pyext/src consolidated: 113 files to 50)

`pmi`, the most Python-heavy module in IMP, is 23,394 lines in **26 files**.
`bff` was 26,900 in **113**, averaging 236 lines against pmi's 900. The tree was
deep because the files were fragmented, so the files were merged rather than
re-filed — no directory scheme turns 236-line modules into a shape IMP has
anywhere. It is now **50 files averaging 539**.

**What decided each case, and it was not taste.** A domain became one flat
module unless it had a reason not to:

* **Something binds its submodules by path.** Of the 99 dotted `IMP.bff.*`
  names referenced across chisurf, imp-tricks, quest and ucfret, only **19
  still resolve** — the rest were broken by the earlier migrations and nobody
  noticed. The 19 concentrate in `quenching` (five), `restraints` (two) and
  `cgdye` (most of it). Those stayed packages.
* **Size.** One `representation.py` would be 4,900 lines and one `io.py` 3,800,
  both larger than anything in IMP (`pmi/macros.py`, the biggest, is 2,803).
  Those stayed packages holding a handful of substantial modules.

**Merging is what finds the collisions.** Three were live, and invisible only
because the two definitions sat in different modules:

* **four `FLRCIF_ITEMS`**, one per dataclass, in `dye/species`, `dye/spectra`,
  `label/site`, `label/quencher`;
* **two `compute_av`** — arrays vs. a PDB plus an fps position.
  `IMP.bff.compute_av` resolved to the second while
  `IMP.bff.representation.av.compute_av` resolved to the first, so *the route
  decided which function you got*. Now `compute_av_from_structure`, with both
  flat names keeping their meaning;
* **`validate` shadowed by a `validate` flag** — `io/fps_schema.validate(payload)`
  merged into a module whose `read_fps_json` takes `validate: bool`, making the
  call `validate(payload)` on a boolean. Now `validate_fps`.

**One rule, learned by breaking it.** A click command is a decorated function,
so `import click` runs at module scope; merging `representation/rotamer/cli.py`
into `rotamer.py` made `import IMP.bff.representation.rotamer` require click.
**A CLI does not live beside the code it drives** — every entry point outside
cgdye is in `cli.py`.

**Three structural tests had quietly stopped covering anything**, and the
consolidation is what exposed it:

* `test_import_discipline` derived its domain list from *directories*, so nine
  domains fell out of the acyclic-graph check the moment they were merged;
* its `len(modules) > 50` tripwire was set when there were 114 modules, and
  fired on the success it was meant to be blind to;
* two "no decorated function anywhere" assertions stood in for "`@njit` is
  gone", and stopped meaning that when `@abc.abstractmethod` and `@property`
  were merged in beside the kernels.

All three now say what they mean. `test_a_directory_only_exists_where_something_binds_its_submodules`
is new and holds the rule above.

## 2026-08-19 (PRD-115 proposed: a differentiable lattice on a shared NN core)

* **PRD-115 written** ([prds/prd-115.md](prds/prd-115.md)): make the field
  solver differentiable — hand adjoint of `diffusion_propagate` (the sweep is
  linear in `cur`, so the reverse pass is the transposed 7-point stencil;
  Smoluchowski self-adjoint up to the `bounds` mask, Itô not; √n
  checkpointing on the `n_out` states), per-voxel features over reachable
  voxels only, and a learned `D(r)`/`k(r)` (later a drift potential) on the
  vendored `MlpCore.h`. Built on PRD-111's finding that `slow_factor` and
  `contact_distance` are one parameter — the data determine the *field* — and
  that every gradient in the stack is a finite difference through a 1–10 s
  solve. Not a PINN: nothing replaces the solver; no torch/jax.
* **The NN side is done in tttrlib** (T-20260819-01, commits `933a4cc7a`,
  `8ac9b0a52`, `cd0cdb000`): header-only, std-only `MlpCore.h` with
  `backward(dL/dy)`, Taylor-augmented passes for losses on `dy/dx` and
  `d²y/dx²`, smooth activations, flat parameters, whole-model `MlpModel` +
  scalers + JSON templated on the JSON type. Vendored verbatim at
  `include/internal/MlpCore.h`; `test/test_vendored_mlpcore.py` compares
  sha256 with `../tttrlib` and compiles `test/cpp_snippets/mlpcore_eval.cpp`
  against bff's own `internal/json.h` under `IMP::bff::internal`, loading a
  model trained by the conda-packaged tttrlib 0.27.0 — predictions to 1e-12,
  `dL/dparams` and `dL/dx` FD-checked. Perf was asserted properly (thread
  CPU time, pre-refactor code from git, interleaved) after wall-clock had
  misled twice: parity.
* Found by the survey, fixed: `benchmark/quenching_identifiability.py`
  imported `IMP.bff.quenching.solver`, which PRD-113 moved to
  `IMP.bff.sampling.smoluchowski`; both quenching benchmarks were
  un-runnable. Noted for PRD-115: `#pragma omp` is inert in bff builds; no
  orientational degree of freedom exists anywhere in bff, so anisotropy is a
  separate PRD, not something `MlpCore.h` covers.

## 2026-08-19 (overnight: seven kernels to C++, and the marshalling that was the real cost)

Seven ports, each gated as equality against the Python it replaced. **The suite
went from 374 s to 91 s.**

| what | before | after | |
|---|---|---|---|
| rotamer interaction energies | 621 ms | 78 ms | 7.9× |
| mean-field pair energies | 267 ms | 6 ms | 47× |
| greedy Olga candidate scoring | 4111 ms | 1311 ms | 3.1× |
| explicit-dye LJ energy | 323 ms | 65 ms | 5.0× |
| hierarchy frame read (rotamer suite) | ~44 s | 6 s | 7× |
| dye library parse (repeat) | 58.7 ms | 0.21 ms | 281× |
| walk on a 101³ grid, 60 steps | 33.8 ms | 0.34 ms | 99× |

**Profiling beat guessing, twice.** I was about to port the rotamer pair
matrices and profiled first: the cost was CIF parsing (`read_dye_library` called
21 times for one file) and SWIG traffic (`_collect_frame`, a dozen crossings per
atom), neither of which the intended port would have touched.

**Marshalling turned out to be the dominant cost**, not arithmetic. Converting a
numpy array into a `std::vector` runs ~34 ns per element; a *returned* vector
becomes a Python tuple, while an **out-parameter** stays a wrapper walked one
`__getitem__` at a time. Three consequences, now applied consistently: flags
ride in the returned array rather than in out-parameters; large inputs take
`(pointer, length)` through numpy.i's `IN_ARRAY1`; small ones stay
`std::vector`.

> **Corrected 2026-08-19, twice.** "Converts at C speed" was wrong: a returned
> `std::vector` is not free. Measured properly — work held constant while the
> array size varies, on two kernels — the return path costs **35–40 ns per
> element** (8–16 ns for SWIG to build the tuple, ~27 ns for `np.asarray` to
> walk it back), not the 66 ns first recorded. An out-parameter walked as a
> proxy is ~340 ns. The right answer for a large array is a numpy **view** over
> the kernel's own buffer (`ARGOUTVIEWM`), which is free.
>
> The second correction is the more useful one. **The input side was the worse
> half, and it was invisible because every adapter did the natural thing.**
> Handing `np.ascontiguousarray(x).ravel()` to a `const std::vector<double>&`
> parameter costs **32–37 ns per element** — *more* than passing the same
> numbers as a Python list (8–13 ns), because SWIG walks the sequence and each
> ndarray element access mints a fresh Python float. So the package was paying
> the most expensive of the three ways in, everywhere, by writing the code that
> looked fastest.
>
> Fixed in one place rather than forty: `IMP_bff.types.i` now carries an `in`
> typemap for `const std::vector<double>&` that bulk-copies a 1-D contiguous
> float64 array (**4.4 ns/element**) and falls back to SWIG's own converter for
> lists, tuples and `VectorDouble`. No kernel signature changed.
> `diffusion_propagate` — four `ng³` grids in, one out — went from **9.54 ms to
> 0.296 ms** of boundary crossing at `ng = 41`: the cost of 53 solver steps, and
> `equilibrium_occupancy` makes up to 200 such calls per solve.
>
> The full table now lives in `include/internal/OutputView.h`; the guards are
> `test/test_zero_copy_views.py` (ownership, aliasing, leaks) and
> `test/test_vector_input_typemap.py` (the fallbacks, and that a strided array
> is not read as contiguous).

**Defects the ports found**, all pre-existing and all invisible to a gate that
compares against the previous implementation:

* **`_fast_convolve_loop` was called but defined nowhere.** It arrived in the
  κ² migration as a call to something that stayed in ChiSurf, so `use_fast=True`
  — the default — raised `NameError` for every input over 1000 products.
* **The mean-field cross term used transposed parameters** whenever `d1 > d2`:
  an `(n₁, n₂)` distance matrix against `(n₂, n₁)` parameters. Flattened lengths
  match, so numpy broadcast it silently. **98.9 % wrong** on dyes whose elements
  differ in radius.
* **The dye-library cache exposed a test writing through a frozen dataclass**
  with `object.__setattr__`, which had been harmless only because every call
  re-parsed.
* **`compute_rotamer_score` validated its `potential` argument after an early
  return**, so a misspelling was accepted whenever the atom selection was empty.
* **No `#pragma omp` in this module has ever run.** IMP's CMake leaves
  `OpenMP_CXX_FLAGS` empty, so every figure above is single-threaded against
  single-threaded numpy — and a lower bound.
  (`okf/validation/openmp_is_not_enabled.md`.)

**And `pyext/` itself**: 21 entries to 6. Fourteen single-line `.i` files folded
into `swig.i-in`, which is IMP's own convention (`core` `%include`s 101 headers
directly); the three with real content renamed `IMP_bff.*.i` so they are
distinguishable from the vendored `numpy.i`. Two name collisions removed —
`fret/kappa2.py`, a 17-line re-export shim, and `distribution.py` beside
`distributions.py` in one directory meaning unrelated things.

## 2026-08-18 (PRD-113 stages 6-7: the output contract, io/, and the public surface)

* **`observables/`** — `LifetimeSpectrum`, `(amplitude, rate)` pairs. The
  contract is stated once and *tested*: no convolution, IRF, pileup, counting
  noise or binning anywhere in the package. C++ for the two real loops —
  evaluating the unconvolved decay, and coarse-graining 30 000 species to 128
  while preserving `sum(a)` and `sum(a·k)` exactly (population and initial
  slope) to 1e-12, max decay deviation 5.3e-5. The claimed error scaling is
  measured, not asserted in prose: halving the bin width quarters the error.
* The `exact` flag is a **checked** claim. A constant quenching rate — nothing
  to average over — must reproduce the spectrum's lifetime, and an alternating
  rate must visibly not. Without the second test the first proves nothing.
* **`io/`** — `fret/io.py` was 682 lines of three unrelated formats. Split on
  AST boundaries into `fps_schema`/`fps`/`fps_legacy`/`structure`. A first pass
  extracting only *function* spans silently dropped `AV_SIMULATION_TYPES`, a
  module-level constant in the gap between two functions; caught by diffing the
  set of top-level names before against after, which is now how the split is
  verified. A second check for unresolved names found the one real coupling —
  `read_fps_json` dispatches into the legacy reader by file extension.
* **`api.py` inverted** — `BY_DOMAIN` authored, `EXPORTS` derived, and the
  naming-families regex deleted. It matched export *names* against a token list
  that grew with every feature and could not tell a misfiled name from a
  well-filed one. Replaced by checks that do not grow: every export's module
  must lie inside its declared domain, no name in two domains, every domain
  importable.
* **The import-discipline test paid for itself on its first run**: a
  cross-domain relative import written to look local, an invalid escape
  sequence nothing surfaced, and a test-module **basename collision**
  (`test/io/test_io.py` vs `test/cgdye/test_io.py`) — pytest imports by
  basename, so one file was reported as a single ERROR line at the bottom of a
  780-test run and its assertions simply did not run. Basename uniqueness is
  now a test, because the failure mode is silent.
* **`dynamics/`** — the Brownian walk, the Smoluchowski solver and the
  excited-state Monte Carlo out of `quenching/`. They were filed under the first
  physics that used them, not under what they are. `pyext/src/` now matches
  PRD-113's target layout. The `api.py` domain check failed the instant they
  moved (eight exports still filed under `quenching`), which is the drift the
  old regex could not see.
* The runtime domain graph is acyclic and reads as a layering:
  `photophysics ← representation ← av ← restraints`, `dynamics` a leaf,
  `quenching → {dynamics, fret}`, `fret → {io, photophysics, representation,
  restraints}`, `cgdye → {dye, fret, io, representation}`.
* Suite **799 passing**, one pre-existing `IMP.em` MRC failure.

## 2026-08-18 (PRD-113 stages 3-5, tranches 7-10: **numba reaches zero**)

* `AVDistance.h`, `DistanceCalibration.h`, `OrientationFactor.h`,
  `BrownianWalk.h`, `PhotonSimulation.h` — the last **23** kernels.
  `_jit.py` deleted; `test/test_no_numba.py` parses every source file so it
  stays deleted. Suite **741 passed**, and it now runs in **134 s** against 614 s.
* **The gate had to be split three ways**, because "identical output to numba"
  is achievable for some kernels and impossible for others:
  * integer/index arithmetic — **bit-exact** (`split_contact_volume`, voxel for
    voxel over 8 random grids).
  * float reductions — **≤ 7.1e-15**, and the residual is *FMA contraction*:
    clang fuses `dg * i + r0`, numba's LLVM does not. Chasing it would be
    chasing codegen, so the tolerance says so.
  * RNG-carrying — **closed forms and distributions only**. No C++ generator
    reproduces numba's stream. Two `test_frozen_reference` tests asserted that
    stream (`n_accepted == 94567`; 6297 photons to 1e-9) and were replaced by
    the distributions they came from, with both spreads recorded in the
    docstrings.
* **Gating on closed forms found two defects that gating on the old numbers
  had preserved for years.**
  * **`kappasq_all` returned half the orientation factor.** `np.random.random(3)`
    normalised fills the cube's positive octant, not the sphere: `⟨κ²⟩ = 0.333`
    where the rigid isotropic limit is `2/3`. C++ gives 0.6663. **Fixed**;
    no callers, so nothing downstream carried it.
  * **The particle walk diffused at 3D.** Per-component step variance was
    `6 D dt` (the total 3-D MSD used as one component's width) where the
    convention is `2 D dt`. `GridDiffusionSolver` gives `⟨x²⟩ = 2Dt` exactly.
    **Fixed**, after tracing it to QuEst's 2019 docstring — a Berkeley teaching
    page's step *magnitude* transcribed as a per-component σ, unchanged through
    Cython → numba → C++ — and after correcting my own **10× unit error**: I
    wrote that `D = 40 Å²/ns` was implausibly high because "free Alexa488 is
    around 4 Å²/ns", and used that to argue the default had absorbed the factor
    of three. `1 Å²/ns = 10 µm²/s`, so 400 µm²/s **is** 40 Å²/ns — the default
    is the free-solution value entered correctly, and there is no calibration
    of `D` anywhere in the stack to have absorbed anything. Defaults unchanged.
    The fix does expose a real inconsistency the wrong width was masking: the
    particle model defaults `D = 40` (free) and the field model
    `free_diffusion = 8` (tethered), 5× apart under one name.
    `okf/validation/particle_vs_field_diffusion.md`.
  * **`gaussian_rmp_to_rda_mean` settled**: σ is the per-component width of the
    *separation vector*, so the correction is `σ²/Rmp` — the two transverse
    components give `⟨|ε⊥|²⟩ = 2σ²`, halved by the expansion. The canonical
    version had `σ²/(2Rmp)`, i.e. half; one function now, and `rmp ≤ 0` returns
    0 rather than the 3.6e11 Å the clamped denominator gave. Anything that
    fitted `sigma_rda` through the old canonical version has a σ that is √2 too
    large.
* **`import IMP.bff.av` as a first import had been raising ImportError** at
  every prior commit — `representation/__init__` → `distribution` → `av` →
  `representation`. Nothing in 676 tests caught it because every test imports
  `IMP.bff` first. Fixed by deferring the `av` import in `distribution.py`,
  which sits above the builder; all eight subpackages now import standalone.
* **De-duplication the port forced**: one `distance_sample_statistics` behind
  three reductions that disagreed at the limits; one Horner evaluator behind two
  `polynomial_transfer`s that read coefficients in opposite orders (the
  ascending caller now reverses and delegates); one `brownian_walk_in_volume`
  behind two near-identical walks, with mobility as a *field* rather than a
  scalar-and-mask pair. `chi2_score` is now literally one object.
* **Still an owner decision**: `gaussian_rmp_to_rda_mean` is `Rmp + s²/Rmp` in
  one module and `Rmp + s²/(2Rmp)` in the other, 45.8 against 45.4. A test
  asserts they still differ, so it cannot be closed by accident.
* **Two improvements the C++ made possible**: `photon_trace` is parallel *and*
  reproducible (each photon's generator seeded from `(seed, index)`, so the
  result does not depend on the scheduler — numba had to choose one or the
  other), and `quenched_decay` takes a seed at all, which the numba never did.

## 2026-08-18 (PRD-113 stage 5, tranches 5-6: maps, FRET traces, the solver — and a 40× regression, found and fixed)

* `QuenchingMap.h`, `FRETRateTrace.h`, `DiffusionSolver.h` — 7 more kernels.
* **Gates**: `slow_near_atoms`, `fret_map`, both grid stampers and the
  Smoluchowski step are **bit-exact**. `quenching_map` differs by 2.6× machine
  epsilon (`np.exp` vs `std::exp`), the trace kernels by 3.2×, the Itô step by
  4e-19. Each was isolated relatively rather than judged on an absolute number:
  3.6e-12 looked alarming until measured against values of order 141.
* **The port made the quenching suite 40× slower — 14 s to 558 s — and that is
  the interesting part.** numba's kernels were `parallel=True` over `prange`;
  the C++ was serial, allocated a grid per call, and crossed the SWIG boundary
  **once per step**. A solve is 1 000–10 000 steps.
* **The fix was to move the loop, not to tune the kernel.** `diffusion_propagate`
  runs the whole ping-pong loop in C++ with OpenMP over x-slabs: 1.33 ms/step →
  **0.023 ms**, 58×. Wiring `run()` alone changed nothing measurable — the hot
  path was `equilibrium()`, which takes tens of thousands of steps. With both
  inside: **488 s → 18.3 s**, now slightly better than the numba it replaced.
* The lesson generalises to the rest of the port: **a kernel called in a tight
  Python loop must take the loop with it**, or the binding cost dominates
  whatever the kernel gains.
* numba: **25 → 23**. Suite **676 passed**.

## 2026-08-18 (PRD-113 stage 5, tranche 3: ASA and the PET rate to C++)

* `SolventAccessibleSurface.h` / `.cpp` — `sphere_points`,
  `solvent_accessible_surface_area` (Shrake–Rupley), `quenching_rate_per_frame`.
* **Gate**: `quenching_rate_per_frame` bit-exact. The other two differ by up to
  2.4e-6 — **and that is the Python being wrong, not the port.** It built the
  golden spiral in `float32`; the C++ agrees with an independent `float64`
  construction to **1e-15**, where the old one was off by 2.4e-6. So the port is
  the more accurate of the two, and the difference is not a tolerance to relax
  later.
* **Another uncovered function, found the same way as the last one.** `pet.py`
  was missing `import IMP` after the delegation, and **139 quenching tests
  passed anyway** — nothing calls `quenching_rate_per_frame`. That is the second
  time in two tranches that a green suite has vouched for a function it never
  runs (`worm_like_chain_linker` was the first).
  `test/quenching/test_ported_kernels.py` now covers all three, including the
  neighbour-cutoff bug the ASA kernel was ported with: an occluder at 5 Å must
  be seen, which the old 2.45 Å cutoff would have missed.
* **The SWIG name collision recurs by construction.** Any C++ function SWIG
  binds into `IMP.bff` takes that flat name, so `api.py` must not also claim it
  — `sphere_points` hit this exactly as `normal_distribution` did. Worth
  expecting for every remaining tranche rather than rediscovering.
* numba: **35 → 33**. Suite **676 passed**.

## 2026-08-18 (PRD-113 stage 5, tranche 2: polymer chains to C++ — and a live breakage found)

* `PolymerChain.h` / `.cpp` — `gaussian_chain_ree`, `gaussian_chain`,
  `worm_like_chain` (the Becker–Rosa–Everaers multi-piece solution, both sides
  of the κ = 0.125 branch), `worm_like_chain_linker`.
* **Gate: agreement to ≤ 3 × 10⁻¹⁷** on every case. One line reported a relative
  failure at κ = 0.05 — on entries of order 1e-18, where a relative tolerance is
  the wrong instrument; the absolute difference is 6e-18.
* **The port found a function that had been raising for a whole commit.**
  `worm_like_chain_linker` was numba-jitted and called `normal_distribution`,
  which tranche 1 had just turned into a C++ delegation numba cannot type — so
  every call raised `TypingError`. **The suite did not notice: 656 tests passed
  with it broken**, because that function had no test at all.
* So its equality gate had no working reference. It was checked against an
  **independent numpy convolution** instead — agreement to 7e-18 — and
  `test/test_polymer_chain.py` now covers it, along with both branches of the κ
  expression and the contour-length cutoff.
* The lesson is about the gate, not the port: an equality gate only proves
  something when the thing being compared against still runs. Where it does not,
  say so and find another reference.
* **IMP concatenates a module's sources**, so a helper in an anonymous namespace
  in two `.cpp` files is a redefinition rather than two private copies. Shared
  as `include/internal/Normalize.h`.
* numba: **40 → 35**. Suite **665 passed**.

## 2026-08-18 (PRD-113 stage 5, tranche 1: distributions to C++)

* **numba must go entirely** (owner). 44 jitted functions across 13 files; this
  is the first four, and it establishes the pattern for the rest: C++ header +
  source, a SWIG `.i`, `Files.cmake`, a Python wrapper keeping the old
  signature, and **an equality gate against the numba before anything is
  deleted**.
* `Distributions.h` / `Distributions.cpp` — `poisson_0toN`,
  `normal_distribution`, `generalized_normal_distribution`,
  `distance_between_gaussian`.
* **Gate: 7 of 9 cases bit-for-bit identical**, the other two within
  **7 × 10⁻¹⁸** — the last bit, from a different summation order in the
  normalisation.
* **Two things the port caught by forcing the formulas to be read.**
  `generalized_normal_distribution` is *not* the exponential-power family its
  name suggests: it skews by transforming the axis,
  `z = −log(1 − κ(x−μ)/σ)/κ`, and evaluates the **standard** normal at `z` —
  `loc` and `scale` are already folded in. My first header documented the wrong
  distribution and the wrong default (`shape=2.0`, when the Python default is
  `0.0`, i.e. no skew). Both corrected against the body rather than the name.
* **A name collision, resolved deliberately.** SWIG binds C++ functions into
  `IMP.bff` directly, so `IMP.bff.normal_distribution` was already taken and
  `api.py`'s lazy hook could never fire — the test caught it. The C++ functions
  **are** the public surface; `IMP.bff.distributions` keeps numpy-returning
  wrappers for use inside the package, and `api.py` records why those names are
  deliberately absent.
* numba: **44 → 40**. Suite **656 passed**.

## 2026-08-18 (PRD-113 stage 4c: interaction terms)

* **The photophysics abstraction exists now**, shaped like a force field's
  terms: a functional form, an **arity**, and parameters looked up by type.
  1-body `RadiativeTerm`, 2-body `PETTerm` (dye × quencher atoms) and `FRETTerm`
  (dye × dye), with the N-body case shaped for but not implemented.
* **Rates add, and that is a first-class property.** `total_rate()` sums the
  channels because parallel deactivation channels add — which is exactly what
  `GridDiffusionSolver` already relies on when it sums a quenching map and a
  FRET map. Burying that inside each observable is how `fret_rate_trace` and
  `fret_rate_map` became two incompatible calls.
* **They are a consolidation, not a seventh layer, and the tests say so.**
  `PETTerm` reproduces `quenching_rate_map`'s law at the voxel centres to
  1e-9 (minus the `1/tau0` floor, which is the radiative term's job — that
  separation is the point of having terms). `FRETTerm` reproduces
  `fret_rate_trace` exactly. The kernels can move underneath the terms later
  without changing an answer.
* **`R0` is derived inside the term** from the two dyes, the medium's refractive
  index and κ² — not passed in. And the term reports
  `used_isotropic_kappa2`, so a representation that cannot resolve orientations
  (an accessible volume) is *told* it fell back to 2/3 rather than silently
  averaged.
* **Terms are representation-agnostic**: they consume
  `IMP.bff.representation.States`, so one implementation serves an AV, a rotamer
  library, a coarse-grained model and an MD trajectory. `needs_orientations`
  says which terms care.
* Suite **656 passed**.

## 2026-08-18 (PRD-113 stage 4b: the distance layers, and two that disagree)

* **`fret/distance.py` → `representation/distance.py`.** A distance between two
  labels is a property of what represents them, not of the FRET engine. It was
  the most complete of the implementations — pair statistics, the RDA histogram,
  the transfer polynomial, orientation-resolved pair geometry — so it is the
  canonical one. `fret/distance.py` is a re-export while the rest of `fret/`
  migrates.
* **I was wrong that there were six duplicates.** That claim came from matching
  *function names* (`comm -12` over `def` lines) without ever comparing
  signatures. Deleting on that basis broke 8 tests, which is how it was caught.
  The truth, checked properly:

  | | verdict |
  |---|---|
  | `fret_efficiency`, `distance_from_fret_efficiency` | identical — removed, re-exported |
  | `chi2_score` | identical but for a parameter name — kept, it is the only symbol `restraints/` imports |
  | `av_pair_statistics` | **different jobs**: array reduction here, sample-and-reduce there |
  | `gaussian_rmp_to_rda_mean`, `polynomial_transfer` | **different answers** |

* **Two same-named, same-signature functions give different numbers.**
  `gaussian_rmp_to_rda_mean(45, 6)` is `Rmp + σ²/Rmp` = **45.8** here and
  `Rmp + σ²/(2 Rmp)` = **45.4** there — a factor of two in the correction term.
  `polynomial_transfer(45, [0, 1, 0.02])` is **85.5** here (ascending
  coefficients) and **45.02** there (descending, `np.polyfit` order).
* **Which is live matters and is now known.** `fret/engine.py` calls the
  `representation.distance` versions, and its `fit_transfer_polynomial` produces
  `np.polyfit` coefficients — which only the descending evaluator reads
  correctly, so that pairing is consistent. The `distance_metrics` versions have
  **no consumers**. This is a trap that had not yet sprung, not a live defect.
* **Which `gaussian_rmp_to_rda_mean` is correct is a physics question** — whether
  `sigma` is the per-component width of an isotropic 3-D cloud or the width of
  the distance distribution — and is left for the owner rather than guessed.
  Both are documented side by side at the top of `distance_metrics.py`.
* Suite **641 passed**.

## 2026-08-18 (PRD-113 stage 4a: `photophysics/`, κ² unified, `spectroscopy/` gone)

* **933 lines of anisotropy modelling were unreachable.**
  `spectroscopy/kappa2.py` — κ² distributions, wobbling in a cone, order
  parameters, the κ²→distance-ratio conversion — had **no consumers at all**: it
  was never exported through `api.py`, and it imported `numba` *directly* rather
  than through `IMP.bff._jit`, so an installation without numba could not even
  load it. Both fixed; it is `photophysics/orientation.py` and is now exported.
* **κ² has one home.** The geometry (`kappa2_from_dipoles`, 44 lines, real
  consumers) and the distributions (933 lines, none) were in unrelated packages
  named the same thing. `photophysics/kappa2.py` and
  `photophysics/orientation.py` now sit side by side. `fret/kappa2.py` is a
  re-export while the rest of `fret/` migrates.
* **`spectroscopy/` is gone** — its two halves were an experiment-side decay
  wrapper (deleted in stage 0) and a forward model (moved here).
* Three docstrings became raw strings: LaTeX in them raised a `SyntaxWarning` on
  every import.
* Suite **641 passed**.

## 2026-08-18 (the sampler: the table transfers, the short test did not)

* PRD-111's posterior harness used the **stretch** move on a posterior whose
  condition number is ~10⁵. ChiSurf benchmarks exactly that case
  (`docs/development/benchmarks.md`) and reports 0.0036 → 0.0109 ESS per
  evaluation from stretch to an adaptive-covariance slice at κ=100, noting the
  gap grows with correlation. The machinery had been surveyed and then not used.
* **An 80-step comparison said the opposite** — slice 6× *worse* per evaluation —
  because a slice step costs 5.7 log-probability evaluations and neither chain
  had mixed (R̂ 2.3–2.5), so the ESS being compared was meaningless.
* **600 steps settles it**, one site, 12 walkers, 3.0 Å:

  | | evaluations | wall | ESS/eval | ESS/s | worst R̂ |
  |---|---|---|---|---|---|
  | stretch | 5 986 | 103 s | 0.00113 | 0.066 | **3.040** |
  | slice (adaptive covariance) | 35 428 | 753 s | **0.00223** | **0.105** | **1.137** |

  2× per evaluation and 1.6× per second — but the number that decides it is R̂:
  **stretch does not mix at all**, and got *worse* from 80 steps (2.333) to 600
  (3.040), so its ESS never meant anything. The slice sampler is the default in
  the harness now.
* The forward model still dominates: one six-site evaluation is ~0.36 s. A
  cheaper solve is a bigger lever than the proposal — the explicit scheme's
  `dt ≤ dg²/(6D)` sets the step count, and an implicit or exponential treatment
  of the *diffusion* term would lift it the way the exponential rate integration
  already lifted the reaction constraint.

## 2026-08-18 (PRD-113: stage 3e was mis-scoped — six distance layers, and they agree)

* **The C++ already has the distance API.** `AV.h` exports `av_distance`,
  `av_distance_quadrature`, `av_random_points`, `av_random_distances` and
  `av_distance_distribution`, operating on `AV` decorators. Stage 3e was written
  as "port `av/_kernels.py` to C++"; the six numba kernels are a **duplicate**
  of that, taking point arrays instead of decorators. It is de-duplication, not
  a port, and it belongs with the other distance layers in stage 4.
* **There are six implementations of dye-pair distances**: `av/_kernels.py`
  (numba), `av/basic.py` (wrapping it), `fret/distance.py`,
  `distance_metrics.py`, `representation/distribution.py` (once on each concrete
  class), and the C++ `av_distance` family.
* **They agree.** T4L A132×A65, 2180×3087 points, 2×10⁵ samples: `<R_DA>` 51.945
  against 51.940 Å, `<R_DA>_E` 51.696 against 51.692 Å, `Rmp` identical to all
  printed digits. The spread is Monte-Carlo sampling noise. **Unlike the
  transposed density, no defect is hiding in this one** — the merge is safe and
  the C++ is the implementation to keep.
* Remaining work for it: C++ overloads taking point arrays, plus SWIG typemaps
  for a *second* 2-D input array — `BFF.types.i` defines one
  (`double *input, int n_input1, int n_input2`) and a two-cloud signature needs
  two. Recorded rather than started, so the change lands in one piece.

## 2026-08-18 (PRD-113 stage 3d: the label/representation collision resolved)

* **`LabelDistribution*` moved to `representation/`.** `label/` means the
  *system* — which dye is attached where — and a label *distribution* is a
  representation of where it can be. One word, two questions; the collision was
  recorded in stage 2 and is now removed. `label/` exports only `Label`,
  `Quencher`, `PETParameters` and the reference tables.
* **They are a third implementation of the `States` surface** (`points`,
  `mean_position`, `n_points`) and a **fourth** copy of the distance layer —
  `dRmp`/`dRDA`/`dRDAE`/`pRDA` duplicated across both concrete classes, and
  again in `BasicAV`, and again in `fret/distance.py`. A `states` view bridges
  them for now; folding four distance layers into one is stage 4, and is a
  behaviour change that does not belong in a move.
* `representation/` now holds every way of describing where a dye is:
  `AccessibleVolume`, `LabelDistributionAV`, `DyeDistributionNormal`, and the
  rotamer library as a sibling — all over `States`.
* Gate 12/12. Suite **636 passed**.

## 2026-08-18 (PRD-113 stage 3c: States, and rotamers stop being an AV)

* **`States` is the abstraction that was missing.** An AV grid point, a rotamer,
  a coarse-grained conformer and an MD frame are the same kind of thing — a state
  the dye can occupy, with a weight, a position and possibly an orientation.
  Everything downstream consumes states, so distances and κ² are written once.
* **`RotamerEnsemble` no longer inherits `AccessibleVolume`.** That inheritance
  was *why* distance code worked for rotamers — by inheritance, not by design —
  and it forced a rotamer library to carry a grid it does not have:
  `density=zeros((0,0,0))`, `grid_step=0.0`, `grid_shape=(0,0,0)`. **No consumer
  ever read them**: every one used `points`, `mean_position`, `n_points` or
  `has_volume`, which is exactly the `States` surface. Both are now siblings.
* **One `AccessibleVolume`.** There were two, in `av/compute.py` and
  `fret/av.py`, with identical seven-field definitions in identical order —
  which is how they came to disagree about the axis order of `density` without
  anything noticing. Now `representation/types.py`, `AccessibleVolume(States)`.
* **`mu` and `orientations` are one array.** The rotamer code keeps its name;
  representation-agnostic code asks for `orientations`. Being the same object
  means no caller can set one and read a stale other — and it is what lets a
  consumer see that a rotamer library resolves dipoles where an AV does not.
* **A regex over-matched and took `_av_imp_bff` with the dataclass.** Caught by
  the suite; redone by locating the class and the next top-level definition
  rather than by pattern. Same lesson as the `api.py` rename map: structured
  code is not a place for a regex.
* Gate still 12/12 byte-identical. Suite **633 passed**.

## 2026-08-18 (PRD-113 stage 3b: one path-map core, two front doors)

* **`representation/pathmap.py` is now the one place `IMP.bff.AV` is driven.**
  Both builders set up the same decorator, call the same `resample()` and read
  the same `PathMap`; both had to learn the same lessons, and each had learned
  only some. The array door knew about the x-fastest voxel order; the structure
  door did not, and shipped a mirrored density (stage 3a). Sharing the core is
  what stops that happening a third time.
* **Gate: 12/12 byte-identical** — 3 sites × 2 resolutions × 2 doors, comparing
  `points`, `density`, `grid_origin`, `attachment_point`, `grid_step` and
  `grid_shape` by SHA-256. Held after each door was repointed and again after
  the dead code was removed. `av/compute.py` 405 → 350 lines.
* **What is deliberately *not* shared, and why.** The two doors disagree on
  conventions: the array door binarises the density to 0/1 float64 and keeps
  IMP's point weights; the structure door keeps raw float32 values and forces
  the weights to one. Those are **behaviour**, and this stage forbids behaviour
  changes — so `PathMapReading` returns the raw readings and each door applies
  its own conventions, with the disagreement recorded in the dataclass docstring
  rather than silently resolved. Choosing between them is a later, deliberate
  decision.
* The core keeps both hard-won details, one from each door: the AV is decorated
  onto its **own** particle with the source passed separately (otherwise the
  resampled map sits at the coordinate origin), and the build is serialised
  while `resample()` is not (which is ~96 % of the wall clock, so the split is
  what makes a threaded caller worth having).
* Suite **631 passed**, same single pre-existing `IMP.em` failure.

## 2026-08-18 (PRD-113 stage 3a: the AV density was transposed)

* **`IMP.bff.compute_av` returned `density` mirrored relative to its own
  `points`.** IMP numbers voxels with *x* fastest (`i = x + nx*y + nx*ny*z`), so
  `fret/av.py`'s C-order reshape into `(nx, ny, nz)` — which makes the *last*
  axis fastest — transposed the volume. Measured on T4L A132 at 1.5 Å: of 2180
  cloud points, **1548 (71 %)** landed on a voxel the density called occupied,
  against **100 %** after transposing. Confirmed at 1.5, 2.0 and 2.5 Å.
* **It survived because a mirrored volume looks right.** Same voxel count, same
  bounding box, same total volume — only the per-voxel relationship to the
  protein is wrong. Nothing but a voxel-by-voxel comparison catches it, and
  nothing was doing one.
* **The fix was already written, in the other builder.** `IMP.bff.av.compute`
  (the array path, from imp-tricks) does
  `reshape((nz, ny, nx), order="C").transpose(2, 1, 0)` and documents the hazard
  in a docstring — *"a transpose still overlaps the truth by ~83 % of its voxels,
  so the point cloud has to be compared against the grid"*. It also reads the
  point cloud from IMP rather than deriving it from the grid, **deliberately, so
  that two independent readings can disagree**. The structure path had neither
  the fix nor the test. This is the case for the merge in one example.
* **Downstream effect measured, not assumed**: on T4L A132 at 2.0 Å the mean
  donor lifetime is **3.5626 ns** correct against **3.3974 ns** as computed —
  4.9 % — with `k_max` 2.54 against 2.81.
* **PRD-110 and PRD-111 read `av.density` and built their maps from it**, so
  every number in them was computed on an accessible region mirrored relative to
  the protein. Their *conclusions* concern the model's sensitivity structure —
  θ unidentifiable from one decay, the multi-site rotation, `contact_distance`
  uninformative, Smoluchowski against Itô — and a mirrored AV is still a
  plausible AV-shaped region, so those are unlikely to hinge on the reflection.
  **That is an expectation, not a measurement**: re-running them is owed.
* **`test/fret/test_av_grid_registration.py`** pins it at three resolutions, plus
  a guard on the guard — asserting that a *transposed* density actually fails the
  test, so it cannot quietly become vacuous.
* No quenching test broke, because they pin against analytic results rather than
  recorded numbers. Suite **631 passed**.

## 2026-08-18 (PRD-113 stage 2: `label/`, PET as a pair property, CIF as the data format)

* **`Label` and `Quencher` replace two untyped dicts.** A labelling site was a
  `source_info` dict threaded through the AV builder, the quenching model and
  every benchmark — and it mixed *where the dye is attached* with *how its
  accessible volume is computed*. `Label` carries only the first half;
  `linker_length`, `radius1..3`, `allowed_sphere_radius` and
  `simulation_grid_resolution` are AV *representation* parameters, and a rotamer
  library has none of them. Pinned by a test.
* **PET is a pair property, not a quencher property** (owner, 2026-08-18). `kQ`
  depends on the redox potentials of *both* partners — a rhodamine, an oxazine
  and a cyanine see the same tryptophan differently — so a `Quencher` carries
  identity only (`comp_id`, the redox-active `atom_ids`, optional location) and
  `PETParameters` is keyed by `(dye, comp_id)`. My first cut put `rate_constant`
  on the `Quencher`, which would have hardened exactly that error.
* **The bundled table says so itself.** Its docstring reads *"reference PET
  parameters for a xanthene dye (Alexa488-like)"* — and the package applies it to
  every dye. `reference_pet_parameters(dye)` keeps that transfer but makes it
  visible: each entry records `measured_for`, so `is_transferred` reports when an
  assumption was made rather than a measurement used.
* **CIF is the data format** (owner). The dye library was 41 CSV files — one
  extinction/QY table and one curve file per dye — with no category names, no
  units and no provenance. Now one `dye_library.cif` (38 dyes, 26 169 spectrum
  points) read through `ihm.format`, the CIF layer the package already uses. The
  CSVs are gone.
* **Two CIF traps, both caught by checking rather than by assuming.**
  `ihm.format.CifReader` reads its keywords from the **handler's `__call__`
  signature**, so a `**kwargs` handler is silently handed nothing — the first
  reader returned zero dyes. And `CifWriter._repr` formats floats with `"%.3f"`,
  which rounded the spectra to three decimals and moved R0 by 0.0065 %
  (5.687415 vs 5.687781 nm); its docstring names the escape hatch — hand it a
  string. After that the CIF route reproduces the CSV R0 **bit for bit**
  (5.687781310385014).
* **`_bff_dye` and `_bff_dye_spectrum` are bff-native, and that is recorded.**
  No dictionary in the stack defines an item for a quantum yield, an extinction
  coefficient or a spectrum — checked across all ten `.dic` files in
  `../mmfdb/src/mmfdb/data`; the only matches are
  `_em_detector.detective_quantum_efficiency` and the NMR spectral categories.
  They should be proposed for `mmfdb_flr_ext.dic`. Identifiers still follow
  flrCIF (`chromophore_name` is `_flr_probe_list.chromophore_name`).
* **`_flr_fret_forster_radius.index_of_refraction` and `.kappa_squared` exist** —
  added by `mmfdb_flr_ext.dic`, not upstream IHM-FLR — so stage 1's new
  `refractive_index` parameter has a dictionary item after all, and a stored R0
  can be reproducible rather than a bare number.
* **A name collision found**: `label/` already meant label *distribution* (a
  representation concept) and now also means label = dye at a site (a system
  concept). `LabelDistribution*` is re-exported meanwhile and moves to
  `representation/` in stage 3.
* Suite **626 passed**, same single pre-existing `IMP.em` failure.

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

## 2026-08-19 — programs into `bin/`, and what the improper gap was hiding

- `scripts/` and `pyext/src/cgdye/scripts/` folded into `bin/`, IMP's place for
  installed programs: `imp_bff_traj2bcif`, `imp_bff_dye_pdb2cif`. The
  conversion half of the dye converter joined `io.structure`; the one-system
  Langevin driver became an example. Adding a program to `bin/` **disables the
  module** until it has a documented section in `README.md` — `setup_module.py`
  enforces it, and the symptom is a missing `_IMP_bff.so`, not a doc warning.
- A package-wide dead-code audit found essentially nothing (2 candidates in 570
  definitions, one a false positive). Its first three revisions each measured
  themselves rather than the code — see 1d1dfb5.
- Chasing a duplicated exclusion derivation surfaced
  [`impropers_are_dropped.md`](impropers_are_dropped.md): the combined-system
  builder never fills `impropers`, which left 81 of them unbuilt across the two
  shipped components, hid a `TypeError` in `scoring.build_dye_restraints`
  (fixed, 192a763), and is the only reason the two exclusion derivations agree.
