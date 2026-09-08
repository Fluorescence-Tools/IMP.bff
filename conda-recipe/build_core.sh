#!/bin/bash
# The core (package `bff`) from standalone/CMakeLists.txt: libimp_bff and the
# headers under $PREFIX, the Python module into site-packages/IMP/bff, all of
# data/ and examples/ under share/IMP/bff (the compiled-in data directory).
set -euo pipefail
mkdir -p build-core
cd build-core
SCCACHE_ARGS=()
if command -v sccache >/dev/null 2>&1; then
  SCCACHE_ARGS=(-DCMAKE_C_COMPILER_LAUNCHER=sccache -DCMAKE_CXX_COMPILER_LAUNCHER=sccache)
fi
cmake ../standalone -G Ninja -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_PREFIX_PATH="$PREFIX" -DCMAKE_INSTALL_PREFIX="$PREFIX" -DCMAKE_INSTALL_LIBDIR=lib \
  -DPython3_EXECUTABLE="$PYTHON" -DIMPBFF_PYTHON_INSTALL_DIR="$SP_DIR" \
  "${SCCACHE_ARGS[@]}"
ninja -j "${CPU_COUNT:-4}"
ninja install
