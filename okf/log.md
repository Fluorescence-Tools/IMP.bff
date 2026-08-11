# Update Log

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
