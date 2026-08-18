#!/bin/bash

mkdir build
cd build

# sccache, when the CI runner puts one on PATH, cuts a rebuilt-from-scratch
# recipe (rattler-build always configures into a fresh build/ dir) down to
# whatever actually changed. A local build without sccache installed just
# skips the launcher flags.
SCCACHE_ARGS=()
if command -v sccache >/dev/null 2>&1; then
  SCCACHE_ARGS=(-DCMAKE_C_COMPILER_LAUNCHER=sccache -DCMAKE_CXX_COMPILER_LAUNCHER=sccache)
fi
cmake .. -DCMAKE_BUILD_TYPE=Release -G Ninja -DCMAKE_INSTALL_LIBDIR=lib -DCMAKE_PREFIX_PATH=$PREFIX -DCMAKE_INSTALL_PREFIX=$PREFIX "${SCCACHE_ARGS[@]}"
ninja install -k 0 -j 4

# Copy examples. -p: ninja install's own doc/examples install rules
# already create this directory for some imp versions.
mkdir -p $PREFIX/share/doc/IMP/examples/bff
cp -r ../examples/* $PREFIX/share/doc/IMP/examples/bff
