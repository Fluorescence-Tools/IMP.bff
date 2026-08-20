# CLAUDE.md

Read [`AGENTS.md`](AGENTS.md) first — it is the authoritative note for agents
in this repo. Key facts:

- **Sibling repos**: `tttrlib`, `imp.bff` (this repo, "bff"), `imp`,
  `imp-tricks`, `chisurf`, `fpsimp` (fpsim), `quest`, `ucfret` are all
  checkouts under the same parent directory (`../`). Layering:
  tttrlib → imp.bff → imp-tricks → chisurf.
- **Knowledge base / PRDs**: this repo has a small OKF bundle at `okf/`
  (start: `okf/index.md`; shared agent board: `okf/agent-board.md`) holding
  **repo-local** knowledge only. Everything cross-stack — PRDs, ecosystem
  map, build workflows — lives in `../chisurf/okf/` (PRD index:
  `../chisurf/okf/prds/index.md`); write cross-stack findings back there,
  never copy pages between bundles. `tttrlib` keeps its own bundle at
  `../tttrlib/okf/` for photon-level concerns.
- **Never commit to the `../imp` checkout** (enforced by a pre-commit hook
  there; see PRD-93).
- Placement rule for code: photons/curves → tttrlib, coordinates → imp.bff,
  neither → chisurf. When input and consumer disagree, the consumer wins.
- **Language rule for this repo: what can be C++ must be C++.** Tests, examples
  and docs are Python; `prototypes/` is exempt; programs go in `bin/`; the rest
  of the Python surface belongs in `pyext/*.i`. `pyext/src` is meant to end up
  empty — see the AGENTS.md section "What language a thing is written in".
