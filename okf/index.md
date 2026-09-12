---
okf_version: "0.2"
---

# IMP.bff knowledge bundle

An [Open Knowledge Format](https://github.com/GoogleCloudPlatform/knowledge-catalog)
bundle for **repo-local** knowledge about `imp.bff` — the coordinate-level
fluorescence-modelling IMP module (accessible volumes, path maps, κ², cgprobe).

**Scope rule:** repo-local knowledge plus **the PRDs that drive imp.bff work**
live here ([`prds/`](prds/index.md)). Cross-stack knowledge — the ecosystem
map, the build workflows, the repository-split rules — stays in
[`../../chisurf/okf/`](../../chisurf/okf/index.md) and is *linked from* here,
never copied. A duplicated page forks the knowledge and both copies rot.
Photon-level concerns live in
[`../../tttrlib/okf/`](../../tttrlib/okf/index.md). A bff feature PRD is
authored in `okf/prds/` of this repo; it is linked (not copied) into the
shared index.

# Agent message board

* [agent-board.md](agent-board.md) - **read this before starting work.**
  Shared coordination channel for agents across tttrlib, chisurf, and imp.bff.
  The board lives in `tttrlib/okf/agent-board.md` and is symlinked here so all
  projects see one board. Claim work, post blockers, hand off.

# Concepts

* [Overview](overview.md) - What imp.bff is, the sibling-repository stack it
  sits in, and where each kind of knowledge lives.
* [`.drot` v10](drot-format.md) - The rotamer-library container, byte for
  byte: what each member holds, how the grids are laid out, which versions
  read, and why the default rung is lossless.
* [The fit graph](log.md) - not a page of its own yet, but the four entries
  that describe it are `log.md` **2026-09-01 (9)** through **(12)**: the whole
  optimisation in C++ (`Minimizer`), the grouping (`JointChiSquared`), the
  group builder in ChiSurf, and the TCSPC decay (`TcspcDecay`) whose kernels
  are tttrlib's. Together they are the arrangement the three repositories now
  fit with -- **bff builds the network, tttrlib computes the curves, ChiSurf
  is not between them.**
* [cgprobe](cgprobe.md) - The explicit-dye layer: what it does, where the physics
  lives (R0/κ² in `fret`), the flat `IMP.bff.*` surface, conventions (CB
  strip), tests and pins.
* [Labelizer correspondence](labelizer-correspondence.md) - What became what in
  the 1:1 port of the Labelizer label-site score (PRD-120): every reference
  module and function mapped to its `labelizer_*`/`Labelizer*` counterpart, where the
  constants went, what was deliberately not ported, and which reference
  defects are reproduced on purpose. Its measured agreement is in
  [validation/labelizer_ab.md](validation/labelizer_ab.md).

# Handovers

Written to be picked up cold. Each states what exists, what the next task is,
and the traps that cost the previous session time.

* [handover-fit-graph-2026-09-01.md](handover-fit-graph-2026-09-01.md) -
  **`GlobalFitModel` onto `JointChiSquared`** (ticket `T-20260901-02`,
  **done** 2026-09-01, see `log.md` 2026-09-01 (11)). A plain ChiSurf fit
  optimises entirely in C++ at 2.8x; a `FitGroup` -- which is what the GUI
  builds -- now does too, at 1.74x, its members' `ChiSquared` nodes under one
  `JointChiSquared` and its shared parameters as port links. Kept because its
  eight traps still apply to anything touching this code.
* [handover-expression-engine.md](handover-expression-engine.md) - the
  string-over-arrays evaluator shared with tttrlib.
* [handover-phase1-2026-08-21.md](handover-phase1-2026-08-21.md) - the
  `%pythoncode`-to-C++ pass (PRD-117).
* [handoff-pythoncode-to-cpp.md](handoff-pythoncode-to-cpp.md) - the method
  that pass followed.

# Subdirectories

* [references](references/index.md) - Pointers into the shared bundles: the
  ecosystem map, build workflow, module conventions, and the PRDs that drive
  imp.bff work.

# Log

* [log.md](log.md) - Chronological update history for this bundle.
