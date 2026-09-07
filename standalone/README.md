# The standalone (IMP-free) build

`standalone/` is what makes the core of IMP.bff build without IMP (PRD-137):

- `include/` -- the `IMP/` tree the standalone build puts on the include
  path *instead of* IMP's: `IMP/bff/bff_config.h` (what IMP's tooling
  generates for the module build, hand-written), `IMP/Object.h`,
  `IMP/Pointer.h`, `IMP/constants.h`, `IMP/algebra/VectorD.h` and
  `Vector3D.h`. They implement the subset of IMP's vocabulary the core uses,
  under IMP's names, so a core source compiles against either.
  `include/Base.h` (the module's one door for IMP's macros) includes them in
  its `IMPBFF_STANDALONE` branch.
- `CMakeLists.txt` -- the build of `libimp_bff` from every source outside
  the connection layer (`src/imp/`, reached in the IMP build by
  `src/ImpLayer.cpp`), plus `src/standalone/` and the vendored parser below.
  It needs Eigen3, cereal, Boost headers, RMF and (optionally) OpenMP, and
  nothing of IMP. `cmake -S standalone -B build -DCMAKE_PREFIX_PATH=$CONDA_PREFIX`.
- `thirdparty/ihm/` -- python-ihm's C mmCIF/BinaryCIF parser, which the
  core's CIF reader uses; the IMP build takes IMP's copy.

The IMP-module build is the top-level `CMakeLists.txt` and is untouched by
any of this. `test/expensive_test_standalone_core_compiles.py` is the gate
that the core needs nothing of IMP; `test/test_base_header.py` that Base.h's
standalone branch reaches only the shims.
