# Update Log

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
