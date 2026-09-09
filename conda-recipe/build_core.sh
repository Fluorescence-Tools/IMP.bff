#!/bin/bash
# The core (package `bff`) from standalone/CMakeLists.txt: libimp_bff and the
# headers under $PREFIX, the Python module into site-packages/IMP/bff, all of
# data/ and examples/ under share/IMP/bff (the compiled-in data directory).
set -euo pipefail

# The large data (data/registry.json) is not in git: fetch it into the source
# copy so the package ships all of data/. IMP_BFF_DATA_CACHE (CI: a directory
# actions/cache keeps between runs) makes this a copy instead of a download.
python "$SRC_DIR/utility/data_registry.py" --fetch "$SRC_DIR/data" --quiet \
  ${IMP_BFF_DATA_CACHE:+--cache "$IMP_BFF_DATA_CACHE"}
mkdir -p build-core
cd build-core
SCCACHE_ARGS=()
if command -v sccache >/dev/null 2>&1; then
  SCCACHE_ARGS=(-DCMAKE_C_COMPILER_LAUNCHER=sccache -DCMAKE_CXX_COMPILER_LAUNCHER=sccache)
fi
# The compute backend is on by default and costs nothing -- it is a plugin
# that links only libc and is never loaded unless a wgpu library is present.
# The switch exists so a build that does not want it cannot be broken by it,
# and is spelt the same way as in build.sh (the IMP module build), which takes
# the same variable.
GPU_ARGS=()
if [ "${IMPBFF_WITH_GPU:-1}" = "0" ]; then
  GPU_ARGS=(-DIMPBFF_WITH_GPU=OFF)
fi
cmake ../standalone -G Ninja -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_PREFIX_PATH="$PREFIX" -DCMAKE_INSTALL_PREFIX="$PREFIX" -DCMAKE_INSTALL_LIBDIR=lib \
  -DPython3_EXECUTABLE="$PYTHON" -DIMPBFF_PYTHON_INSTALL_DIR="$SP_DIR" \
  "${SCCACHE_ARGS[@]}" "${GPU_ARGS[@]}"
ninja -j "${CPU_COUNT:-4}"
ninja install
