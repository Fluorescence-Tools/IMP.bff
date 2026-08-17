# Update Log

## 2026-08-17
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
