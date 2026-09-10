#!/bin/bash

# The large data (data/registry.json) is not in git: fetch it into the source
# copy so the package ships all of data/. IMP_BFF_DATA_CACHE (CI: a directory
# actions/cache keeps between runs) makes this a copy instead of a download.
python "$SRC_DIR/utility/data_registry.py" --fetch "$SRC_DIR/data" --quiet \
  ${IMP_BFF_DATA_CACHE:+--cache "$IMP_BFF_DATA_CACHE"}

# The module's swig run %includes RMF.i (through IMP's rmf fragments), and
# IMP's swig share dir is on the include list while RMF's own is not -- the
# built IMP's build_info/RMF descriptor carries the path, but it does not
# reach the flags. Place the fragments where the flags already point.
if [ -d "$PREFIX/share/RMF/swig" ]; then
  cp -R "$PREFIX/share/RMF/swig/." "$PREFIX/share/IMP/swig/"
fi

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
# IMP's module tooling has no per-module options, so the compute door is
# switched with an ordinary compiler flag. IMPBFF_WITH_GPU=0 in the
# environment compiles the loader out and skips the plugin below.
GPU_ARGS=()
if [ "${IMPBFF_WITH_GPU:-1}" = "0" ]; then
  GPU_ARGS=(-DCMAKE_CXX_FLAGS=-DIMPBFF_WITH_GPU=0)
fi
cmake .. -DCMAKE_BUILD_TYPE=Release -G Ninja -DCMAKE_INSTALL_LIBDIR=lib -DCMAKE_PREFIX_PATH=$PREFIX -DCMAKE_INSTALL_PREFIX=$PREFIX "${SCCACHE_ARGS[@]}" "${GPU_ARGS[@]}"
# Memory-aware job count, as conda-recipe/imp/build.sh has it: the unity
# build file bff_all.cpp is one large translation unit; one job per 6 GB of
# available memory, capped by the core count, IMPBFF_BUILD_JOBS overriding.
JOBS=${IMPBFF_BUILD_JOBS:-$(python -c "import os;c=os.cpu_count() or 1;m=sum(int(l.split()[1]) for l in open('/proc/meminfo') if l.startswith('MemAvailable'))//1024 if os.path.exists('/proc/meminfo') else 0;print(max(1,min(c,m//6)) if m else max(1,c//2))")}

echo "build: -j ${JOBS}"
ninja install -k 0 -j ${JOBS}

# The compute backend, for the module package as well as the wheel.
#
# It is a run-time-loaded plugin: one C file that links nothing but libc and
# resolves every wgpu entry point by dlopen from a path the Python side finds.
# IMP's module tooling cannot express a second shared library --
# `modules/*/CMakeLists.txt` is generated from `tools/build/cmake_templates`
# and says "any changes will be lost", and no IMP module declares an
# `add_library` of its own -- but it does not need to. Nothing links this and
# it is located by path at run time, so compiling it here in one command and
# dropping it beside the module is enough: `IMP.bff.enable_gpu()` then finds
# it exactly as it does in the wheel, and still only uses it if a wgpu library
# (wgpu-py, or IMP_BFF_WGPU_LIBRARY) is present.
#
# A failure here is not fatal. The CPU path is the default, and the whole
# point of a plugin is that its absence costs nothing.
if [ "$(uname)" = "Darwin" ]; then
  _gpu_lib=libimp_bff_wgpu.dylib
  _gpu_extra=""
else
  _gpu_lib=libimp_bff_wgpu.so
  _gpu_extra="-ldl"
fi
_impbff_py=$(ls -d "$PREFIX"/lib/python*/site-packages/IMP/bff 2>/dev/null | head -1)
if [ "${IMPBFF_WITH_GPU:-1}" = "0" ]; then
  echo "IMP.bff: IMPBFF_WITH_GPU=0; no compute backend built"
elif [ -n "$_impbff_py" ] && [ -f "$SRC_DIR/gpu/imp_bff_wgpu.c" ]; then
  if "${CC:-cc}" -std=c11 -O2 -fPIC -shared -fvisibility=hidden \
       -I"$SRC_DIR/gpu" "$SRC_DIR/gpu/imp_bff_wgpu.c" \
       -o "$_impbff_py/$_gpu_lib" $_gpu_extra; then
    echo "IMP.bff: compute backend built into $_impbff_py"
  else
    echo "IMP.bff: no compute backend built; the kernels stay on the CPU"
    rm -f "$_impbff_py/$_gpu_lib"
  fi
else
  echo "IMP.bff: no module directory found; no compute backend built"
fi

# Copy examples. -p: ninja install's own doc/examples install rules
# already create this directory for some imp versions.
mkdir -p $PREFIX/share/doc/IMP/examples/bff
cp -r ../examples/* $PREFIX/share/doc/IMP/examples/bff
