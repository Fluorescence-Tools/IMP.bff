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
  IMP module has to follow.
* [change-tracking](../../../chisurf/okf/workflows/change-tracking.md) - the
  loop every material change follows.
* [PRD index](../../../chisurf/okf/prds/index.md) - all PRDs. The ones that
  drive imp.bff work directly: [PRD-93](../../../chisurf/okf/prds/prd-93.md)
  (the four-repository split), [PRD-94](../../../chisurf/okf/prds/prd-94.md)
  (AV Gaussian moments), [PRD-95](../../../chisurf/okf/prds/prd-95.md)
  (cereal/JSON serialization), [PRD-96](../../../chisurf/okf/prds/prd-96.md)
  (fpsimp as a consumer), [PRD-97](../../../chisurf/okf/prds/prd-97.md)
  (FRET docking moves here).

# The tttrlib bundle (`../../tttrlib/okf/`)

* [tttrlib bundle root](../../../tttrlib/okf/index.md) - photon-level
  concerns: bindings, PTO containers, HMM, testing.
