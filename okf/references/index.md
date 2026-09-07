---
okf_version: "0.2"
---

# References

Pointers into the shared bundles. These pages are the authority; do not copy
their content here.

# The shared ChiSurf bundle (`../../chisurf/okf/`)

* [imp-ecosystem](../../../chisurf/okf/references/imp-ecosystem.md) - which of
  the four IMP-stack checkouts owns which symbol, and the scope boundaries
  between them. **A name missing from `IMP.bff` is more likely shadowed than
  absent**; that page says how to tell.
* [imp-local-build](../../../chisurf/okf/workflows/imp-local-build.md) -
  building IMP from source against the `arm64` conda env, and wiring the build
  tree in so rebuilds land with no re-install.
* [imp-module-conventions](../../../chisurf/okf/subsystems/imp-module-conventions.md) -
  the SWIG, cereal, generated-file and deprecation conventions an out-of-tree
  IMP module has to follow. Includes "Wrapping numpy arrays": use the stock
  numpy.i suites via `%apply`, never a hand-written typemap.
* [change-tracking](../../../chisurf/okf/workflows/change-tracking.md) - the
  loop every material change follows.
* [PRD index](../../../chisurf/okf/prds/index.md) - all PRDs. The ones that
  drive imp.bff work directly: [PRD-93](../../../chisurf/okf/prds/prd-93.md)
  (the four-repository split), [PRD-94](../../../chisurf/okf/prds/prd-94.md)
  (AV Gaussian moments), [PRD-95](../../../chisurf/okf/prds/prd-95.md)
  (cereal/JSON serialization), [PRD-96](../../../chisurf/okf/prds/prd-96.md)
  (fpsimp as a consumer), [PRD-97](../../../chisurf/okf/prds/prd-97.md)
  (FRET docking moves here). bff-local PRDs live in
  [../prds/index.md](../prds/index.md): [PRD-105](../prds/prd-105.md) (AV
  evaluation on the semi space-fixed lattice — implemented 2026-08-17),
  [PRD-106](../prds/prd-106.md) (strip mech).

# The FPS toolkit (`../../../chisurf/junk/fps`), specified

Two pages that let the C# FPS toolkit be ported without opening the C# again.
Both cite `File.cs:line` for every non-obvious claim and both carry a defect
index — reading them is cheaper than rediscovering what they found. See
[PRD-121](../prds/prd-121.md).

* [FPS export formats](fps-export-formats.md) — the six things FPS writes
  after a run (per-result PyMOL script, Overlay, OverlayStates, R table, chi2
  table, `SimulationResults.bin`), byte by byte, with the arithmetic behind
  each number and what is lossy to regenerate. **`SimulationResult.RMSD` has a
  sign error**, and `BestFitRotation` is never computed by any `Save*` method.
* [FPS sampling and protocol](fps-sampling-and-protocol.md) — error estimation
  (a parametric bootstrap on a **sign-split** normal), Metropolis sampling
  (`rkT` is `1/kT`, so larger is colder), and every protocol knob that changes
  a docked answer, with all five modes' shipped defaults. **Refinement
  recomputes the accessible volumes in the docked complex**, which is how FPS
  reconciles per-subunit AVs with inter-subunit occlusion.

# The tttrlib bundle (`../../tttrlib/okf/`)

* [tttrlib bundle root](../../../tttrlib/okf/index.md) - photon-level
  concerns: bindings, PTO containers, HMM, testing.
