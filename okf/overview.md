---
type: Codebase
title: IMP.bff
description: Out-of-tree IMP module for coordinate-level fluorescence modelling — accessible volumes (AV), path maps, kappa-squared, and coarse-grained dye (cgprobe) sampling.
resource: https://github.com/fluorescence-tools/imp.bff
tags: [imp, fluorescence, fret, accessible-volume, kappa2, cgprobe, cpp, swig]
timestamp: '2026-08-11T00:00:00Z'
---

# What it is

`imp.bff` (the *Bayesian Framework for Fluorescence*) is an out-of-tree
[IMP](https://integrativemodeling.org) module. Its input is **coordinates**:
given a structure or ensemble, it models what fluorescence experiments would
see — accessible-volume dye clouds, path-map densities, orientation factors
(κ²), and coarse-grained explicit-dye sampling.

# The sibling stack

All checkouts live under one parent directory. The layering is
**tttrlib → imp.bff → imp-tricks → chisurf**; the placement test is *what is
the input* (photons/curves → tttrlib, coordinates → imp.bff, neither →
chisurf), and when input and consumer disagree, the consumer wins. See
[PRD-93](../../chisurf/okf/prds/prd-93.md).

| Checkout | Role |
|---|---|
| `../tttrlib` | Photon-level data. Own OKF bundle at `../tttrlib/okf/`. |
| `../imp.bff` | **This repo** ("bff"). Coordinate-level modelling. |
| `../imp` | The IMP checkout we build against. **Never commit there** (pre-commit hook enforces it). |
| `../imp-tricks` | Sibling project: downstream IMP experiments. Shadows `IMP/bff/` — a subpackage resolves to the *first* matching directory (PRD-93). |
| `../chisurf` | Fitting glue, GUI, data-IO. Hosts the shared bundle and all PRDs. |
| `../fpsimp` | Sibling project ("fpsim"): fluorescent-protein sampling, a consumer of imp.bff (PRD-96). |
| `../quest` | Sibling project in the fluorescence ecosystem. |
| `../ucfret` | Sibling project in the fluorescence ecosystem. |

# Where knowledge lives

* **This bundle** (`okf/`): repo-local facts only — build quirks specific to
  this module, design decisions taken inside this tree.
* **`../chisurf/okf/`**: everything cross-stack — the
  [ecosystem map](../../chisurf/okf/references/imp-ecosystem.md),
  [local IMP build workflow](../../chisurf/okf/workflows/imp-local-build.md),
  [module conventions](../../chisurf/okf/subsystems/imp-module-conventions.md),
  and the [PRDs](../../chisurf/okf/prds/index.md).
* **`../tttrlib/okf/`**: photon-level concerns and the shared
  [agent board](agent-board.md).

Every material change follows the
[change-tracking loop](../../chisurf/okf/workflows/change-tracking.md):
update the matching concept, append a dated bullet to the bundle's `log.md`,
leave a resume point, commit locally.
