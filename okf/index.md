---
okf_version: "0.2"
---

# IMP.bff knowledge bundle

An [Open Knowledge Format](https://github.com/GoogleCloudPlatform/knowledge-catalog)
bundle for **repo-local** knowledge about `imp.bff` — the coordinate-level
fluorescence-modelling IMP module (accessible volumes, path maps, κ², cgdye).

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
* [cgdye](cgdye.md) - The explicit-dye layer: what it does, where the physics
  lives (R0/κ² in `fret`), the flat `IMP.bff.*` surface, conventions (CB
  strip), tests and pins.

# Subdirectories

* [references](references/index.md) - Pointers into the shared bundles: the
  ecosystem map, build workflow, module conventions, and the PRDs that drive
  imp.bff work.

# Log

* [log.md](log.md) - Chronological update history for this bundle.
