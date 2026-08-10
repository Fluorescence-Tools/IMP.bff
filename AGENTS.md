# AGENTS.md

Notes for anyone — human or agent — working in this repository.

## The knowledge base lives in ChiSurf

This repository has **no OKF bundle of its own**. The durable, agent-readable
knowledge for the whole fluorescence-modelling stack — this repo, `imp.bff`,
`imp-tricks`, `chisurf` and how they fit together — lives in one place:

    ../chisurf/okf/

Read it before starting work, and **write findings back into it** rather than
into a new bundle here. A second bundle would fork the knowledge and both
copies would rot.

Start with:

* [`../chisurf/okf/index.md`](../chisurf/okf/index.md) — the bundle root.
* [`../chisurf/okf/references/imp-ecosystem.md`](../chisurf/okf/references/imp-ecosystem.md) —
  which of the four checkouts owns which symbol, and the scope boundaries
  between them. **A name missing from `IMP.bff` is more likely shadowed than
  absent**; that page says how to tell.
* [`../chisurf/okf/workflows/imp-local-build.md`](../chisurf/okf/workflows/imp-local-build.md) —
  building IMP from source against the `arm64` conda env, and wiring the build
  tree in so rebuilds land with no re-install.
* [`../chisurf/okf/subsystems/imp-module-conventions.md`](../chisurf/okf/subsystems/imp-module-conventions.md) —
  the SWIG, cereal, generated-file and deprecation conventions an out-of-tree
  IMP module has to follow.
* [`../chisurf/okf/workflows/change-tracking.md`](../chisurf/okf/workflows/change-tracking.md) —
  the loop every material change follows: update the matching concept, append a
  dated bullet to `okf/log.md`, leave a resume point, commit locally.

`tttrlib` keeps its own bundle at `../tttrlib/okf/` for photon-level concerns,
and the two share one agent message board.
