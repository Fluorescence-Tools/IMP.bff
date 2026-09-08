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
- `pyext/` -- the Python module without IMP (`-DIMPBFF_PYTHON=ON`, the
  default; needs SWIG 4 and numpy). `IMP_bff_standalone.i` wraps the very
  same `pyext/IMP_bff.core.i` the IMP build wraps (the IMP build adds
  `IMP_bff.layer.i` after it); `IMP_bff_standalone.macros.i` supplies what
  IMP's kernel interface supplies there: the `IMP_SWIG_VALUE/OBJECT/
  DIRECTOR/VALUE_SERIALIZE_IMPL` equivalents, `Vector3D`/`Vector4D` as
  tuples, IMP's exception family (`IMP.bff.ValueException` is a
  `ValueError`, `IOException` an `IOError`, ... -- the same classes and the
  same `std::` mapping as IMP's kernel), IMP's director registry, and
  `%implicitconv`. The build writes `python/IMP/bff/{__init__.py,
  _IMP_bff.so}`; `python/IMP/__init__.py` is a build-tree-only marker so
  that this `IMP/` wins over an installed IMP on `PYTHONPATH`.

Running it from the build tree:

    cmake -S standalone -B build -DCMAKE_PREFIX_PATH=$CONDA_PREFIX
    cmake --build build -j
    IMP_BFF_DATA=$PWD/data PYTHONPATH=build/python python -c "import IMP.bff"

`IMP_BFF_DATA` / `IMP_BFF_EXAMPLES` point `get_data_path()` at a checkout
(a PATH-like list of directories, first hit wins); an installed core finds
`share/IMP/bff/` on its own. The test suite serves both builds: with no
`IMP.atom` importable, `test/conftest.py` leaves out every file that names
another IMP module (or a connection-layer name) and says so in the report
header -- `pytest test` then runs the core's lane.

## Packaging (PRD-137 step 6d)

Two packages, one import name:

- **`bff`** -- this core. `pip install bff` (wheels built by cibuildwheel
  from the root `pyproject.toml` over this directory, scikit-build-core;
  `utility/wheel_deps.sh` gives a wheel image cereal and an RMF 1.7.1 built
  without its HDF5 backend, since neither is on PyPI) and `conda install bff`
  (`conda-recipe/recipe-core.yaml`, `build_core.sh`). No `imp` anywhere.
- **`imp.bff`** -- the IMP module build, conda only (`conda-recipe/recipe.yaml`,
  today's package): the same core plus the connection layer.

Both install `site-packages/IMP/bff/` and are alternatives, not companions
(the recipes constrain each other; `IMP.bff.get_build()` says which one is
present, and the core warns once at import when an IMP kernel sits beside it).

The wheel carries the small part of `data/` in `IMP/bff/data` and lists
`rotamer_library/` and `cgprobe/` (62 MB) in `data/registry.json` with their
sha256 sums (`utility/data_registry.py`, regenerate after touching either
directory). Those are fetched on first use from `IMP.bff.DATA_URL`
(`https://www.peulen.xyz/downloads/imp-bff-data/`, the two directories
uploaded with their paths kept; `IMP_BFF_DATA_URL` overrides) into
`pooch`'s cache (`IMP_BFF_CACHE` overrides), through `get_data_path()`,
`fetch_data()` or the `imp_bff_fetch_data` script. Conda packages ship all
of `data/` under `share/IMP/bff/data` and never fetch.

The IMP-module build is the top-level `CMakeLists.txt` and is untouched by
any of this. `test/expensive_test_standalone_core_compiles.py` is the gate
that the core needs nothing of IMP; `test/test_base_header.py` that Base.h's
standalone branch reaches only the shims.
