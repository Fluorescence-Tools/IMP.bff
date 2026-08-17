---
okf_version: "0.2"
---

# PRDs driving imp.bff

Feature PRDs that drive imp.bff work are authored in this bundle. Cross-stack
PRDs (repository split, build discipline, ecosystem-wide decisions) stay in
[`../../../chisurf/okf/prds/`](../../../chisurf/okf/prds/index.md), per the
scope rule in [`../index.md`](../index.md).

## bff-local PRDs

* ✅ [PRD-105: AV evaluation on a semi space-fixed lattice — rolling windows, one shared obstacle map, deterministic distances](prd-105.md) — **status: implemented (2026-08-17; three staged commits on `dev`).** Global absolute lattice is the default (`space_fixed=False` = deprecated legacy opt-out, byte-identical single-evaluation pins), minimal linker-defined windows rolling by whole voxels, one shared occupancy raster per (spacing, extra-radius) class (`AVOccupancyMap`/`AVOccupancyRegistry`, 4 classes on T4L), cold Dijkstra unchanged, deterministic lattice-quadrature distances (K=100, block centroids + second moments: max error < 0.005 Å vs exact on T4L @2 Å; MC 50k ~0.25 Å range). Measured T4L screening after the same-day perf pass (threaded compute phases, inline neighbour enumeration, locality-aware skip): **~0.8–1.0 ms/frame** vs 36.6 corrected legacy / 37.7 pre-PRD (~40×; three perf passes: threads, bounded exact bucket-queue search, SoA tiles, symmetric 26-neighbour stencil — a physics fix, the historical stencil leaked through one-voxel walls in −x/−y/−z, AVs 15–25 % smaller now — and K=50); same-frame repeats 0.18 ms. Two latent legacy defects found and fixed on the way (stale Dijkstra edges after the first evaluation from a swapped `resize` argument; carved tiles never restored). 6-mode suite in `test/test_av_lattice.py`; benchmark in `benchmark/benchmark_av_screening.py`.
* ✏️ [PRD-106: strip mech — one selection-driven mechanism to strip side chains (and whole residues / atom masks) from the obstacle set](prd-106.md) — The fps.json `strip_mask` field is authored in the shipped examples (T4L, TG2/flex) with a PyMOL-shaped grammar yet **read by no code**; the AV build substitutes clearance for stripping (`allowed_sphere_radius`, `fret/av.py:124-132`), `_strip_residue_atoms` (`fret/av.py:557`) is dead code, and cgdye's `strip_sidechain_at_site` uses a hard-coded keep-set that disagrees with the AV convention. Scopes one strip engine in `fret` (grammar-faithful, non-destructive, two outputs — stripped hierarchy + filtered obstacle array), wiring `strip_mask` into the AV build and cgdye labelling, retiring the dead/duplicated strippers. Consumers: AV build, cgdye `attach_dyes(strip_site_sidechain=True)`, the chimol `_stripped_pdb_for` seam (T-20260815-01).