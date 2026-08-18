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

# -fpermissive: the installed imp package's own SWIG glue
# ($PREFIX/include/IMP/internal/swig_helpers_base.h) has a Python-2 fallback
# branch (PyString_Check/PyInt_Check/...) that is dead code under Python 3
# but is still name-resolved at template-definition time under strict
# two-phase lookup, which newer GCC (13+) now enforces as a hard error --
# GCC's own diagnostic names this exact flag as the fix. Upstream IMP's bug,
# not this tree's; drop this once IMP's headers stop referencing removed
# Python 2 C API symbols.
cmake .. -DCMAKE_BUILD_TYPE=Release -G Ninja -DCMAKE_INSTALL_LIBDIR=lib -DCMAKE_PREFIX_PATH=$PREFIX -DCMAKE_INSTALL_PREFIX=$PREFIX -DCMAKE_CXX_FLAGS=-fpermissive "${SCCACHE_ARGS[@]}"
ninja install -k 0 -j 4

# Copy examples
mkdir $PREFIX/share/doc/IMP/examples/bff
cp -r ../examples/* $PREFIX/share/doc/IMP/examples/bff
