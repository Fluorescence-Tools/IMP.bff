# AGENTS.md

Notes for anyone — human or agent — working in this repository.

## Sibling repositories — one stack, one parent directory

This repo is one checkout in a multi-repository stack. All of these live as
**siblings under the same parent directory** (`../` from here):

| Checkout | Also called | Role |
|---|---|---|
| `../tttrlib` | — | Photon-level data (TTTR streams, histograms, correlation). Bottom of the stack. |
| `../imp.bff` | **bff** (this repo) | Coordinate-level fluorescence modelling as an IMP module (AVs, path maps, κ², cgprobe). |
| `../imp` | — | The IMP checkout imp.bff builds against. **Never commit there.** |
| `../imp-tricks` | — | Downstream IMP experiments; shadows `IMP/bff/` — see PRD-93. |
| `../chisurf` | — | Fitting-model glue, GUI, data-IO. Top of the stack. Hosts the shared OKF knowledge bundle. |
| `../fpsimp` | fpsim | Fluorescent-protein sampling; a *consumer* of imp.bff (PRD-96), not a second implementation. |
| `../quest` | — | Sibling project in the same fluorescence ecosystem. |
| `../ucfret` | — | Sibling project in the same fluorescence ecosystem. |

The layering is **tttrlib → imp.bff → imp-tricks → chisurf**; the placement
test is *what is the input* (photons/curves → tttrlib, coordinates → imp.bff,
neither → chisurf). See `../chisurf/okf/prds/prd-93.md`.

## The compute/display line — the rule above the language rule

> **Keep it all in bff and tttrlib. Only the things that get displayed move
> into chisurf/Python — so that the data stay local.** (owner, 2026-09-01)

Two halves, and the second is the one that gets forgotten. Computation is
C++; **and the data stay where the computation is.** It is not enough for the
arithmetic to be in C++ if the arrays are marshalled back and forth to drive
it. The test is *who consumes the value*: a human looking at it (a plot, a
table, a fitted number) crosses; another computation does not.

A value crossing once per `run()` is fine. The same value crossing once per
*iteration* is what the rule exists to prevent, and the reason is measured
rather than assumed — replacing an optimiser but keeping the callback bought
1.02x, and wrapping the callback in a C++ loop was a *regression*
(`okf/log.md` 2026-09-01 (9)).

The full statement — the decision procedure, what legitimately stays in
chisurf, where it is enforced, and the measured gaps still open against it —
is cross-stack and lives in the shared bundle:
[`../chisurf/okf/architecture/compute-display-line.md`](../chisurf/okf/architecture/compute-display-line.md).

## What language a thing is written in

**What can be C++ must be C++.** Not for speed — to minimise the number of
places a value changes language, because every one of those is a marshalling
convention that can drift. This repository has paid for that twice already: a
density that a C-order reshape returned transposed, and a `radius1` the C++
wrote into all three radii while the Python believed otherwise. Both survived
because the two sides each held their own idea of the same object.

The exceptions are exact:

* **tests, examples and documentation** are Python;
* **`prototypes/`** is exempt entirely — it is where an idea is tried, not
  where it ships;
* **programs** live in `bin/`, and a program is Python because a `click`
  command is;
* what is genuinely neither — a `%pythoncode` shim restoring a caller's array
  shape, a table of dictionary item names — goes in a `pyext/*.i` file.

The consequence is that **`pyext/src` is meant to end up empty**: kernels and
value types to C++, programs to `bin/`, the remainder into the `.i` files. See
[`okf/log.md`](okf/log.md) for what has moved and what is next, and
[`pyext/src/README.md`](pyext/src/README.md) for the shape that is left.

## Two bundles, one scope rule

This repository has a small OKF bundle of its own at [`okf/`](okf/index.md) —
start at [`okf/index.md`](okf/index.md), read the shared
[`okf/agent-board.md`](okf/agent-board.md) before starting work, and log
changes in [`okf/log.md`](okf/log.md). It holds **repo-local knowledge**:
build quirks, design decisions specific to this tree, and — as of
2026-08-16 — **the feature PRDs that drive imp.bff work**, authored in
[`okf/prds/`](okf/prds/index.md).

Everything cross-stack — the ecosystem map, cross-repo PRDs (the repository
split, build discipline, ecosystem-wide decisions), the build workflows, the
repository-split rules — stays in the shared bundle:

    ../chisurf/okf/

The split is *what the PRD decides about*: a decision that is about
`IMP.bff` alone lives here in `okf/prds/`; a decision that spans repositories
lives in `../chisurf/okf/prds/`. Link across bundles rather than copying. A
duplicated page forks the knowledge and both copies rot.

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
* [`okf/prds/index.md`](okf/prds/index.md) — the PRDs driving imp.bff work:

`tttrlib` keeps its own bundle at `../tttrlib/okf/` for photon-level concerns,
and the two share one agent message board.
