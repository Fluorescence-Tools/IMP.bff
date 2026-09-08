#!/usr/bin/env bash
# The two dependencies a wheel image does not have (PRD-137 step 6d):
# cereal (headers) and RMF 1.7.1 built without its deprecated HDF5 backend.
# Installs both under the prefix given ($1, default /usr/local), where
# standalone/CMakeLists.txt finds them through CMAKE_PREFIX_PATH. Boost
# (RMF links filesystem/thread/program_options/iostreams) and Eigen come from
# the image's package manager -- see [tool.cibuildwheel] in pyproject.toml.
set -euo pipefail
PREFIX="${1:-/usr/local}"
RMF_VERSION="${RMF_VERSION:-1.7.1}"
# RMF is off by default here, as it is in the wheel (PRD-139): it links
# Boost.Iostreams, and that drags Boost.Regex and ICU -- 40 MB of bundled
# libraries for four I/O functions. WITH_RMF=1 builds it.
WITH_RMF="${WITH_RMF:-0}"
CEREAL_VERSION="${CEREAL_VERSION:-1.3.2}"
# IMP, when the wheel is to carry the connection layer (IMPBFF_WITH_IMP,
# PRD-139). Off by default: the plain wheel is the IMP-free core.
IMP_VERSION="${IMP_VERSION:-2.25.0}"
WITH_IMP="${WITH_IMP:-0}"
# What the layer actually names, measured: kernel, algebra, display,
# score_functor, core, container, atom (10.1 MB), plus em and rotamer for
# one function each (1.4 MB more). Everything else of IMP is left out.
IMP_MODULES="${IMP_MODULES:-kernel algebra display score_functor core container atom em statistics rotamer}"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT
mkdir -p "$PREFIX/include" "$PREFIX/lib"

# cereal: headers only
curl -sSL "https://github.com/USCiLab/cereal/archive/refs/tags/v${CEREAL_VERSION}.tar.gz" | tar -xz -C "$WORK"
cp -R "$WORK/cereal-${CEREAL_VERSION}/include/cereal" "$PREFIX/include/"

if [ "$WITH_RMF" = "1" ]; then
# RMF: only the library target; its CMake also builds bin/, swig/, tests,
# which a wheel does not need and whose install rules would drag them in, so
# the library and the headers (the source tree's and the generated ones) are
# copied by hand.
curl -sSL "https://github.com/salilab/rmf/archive/refs/tags/${RMF_VERSION}.tar.gz" | tar -xz -C "$WORK"
RMF_SRC="$WORK/rmf-${RMF_VERSION}"
# RMF's CMakeLists runs FindHDF5 directly (include(FindHDF5)), which no
# CMake switch can turn off, and links HDF5 whenever it is on the machine.
# The wheel wants the avro backend only: make the lookup come back empty.
sed -i.bak -e 's/^include(FindHDF5)$//' -e 's/^find_package(HDF5)$/set(HDF5_INCLUDE_DIRS "HDF5_INCLUDE_DIRS-NOTFOUND")/' "$RMF_SRC/CMakeLists.txt"
grep -q 'HDF5_INCLUDE_DIRS-NOTFOUND' "$RMF_SRC/CMakeLists.txt"
cmake -S "$RMF_SRC" -B "$WORK/rmf-build" -DCMAKE_BUILD_TYPE=Release -DRMF_DEPRECATED_BACKENDS=0 \
      -DCMAKE_INSTALL_PREFIX="$PREFIX" -DCMAKE_POSITION_INDEPENDENT_CODE=ON
cmake --build "$WORK/rmf-build" --target RMF-lib -j "$(nproc 2>/dev/null || sysctl -n hw.ncpu)"
cp -R "$RMF_SRC/include/RMF" "$PREFIX/include/"
cp -R "$WORK/rmf-build/include/RMF" "$PREFIX/include/"
cp "$WORK/rmf-build/include/RMF.h" "$PREFIX/include/"
cp -R "$WORK/rmf-build/lib/"libRMF* "$PREFIX/lib/"
# RMF's build leaves each dylib's install name pointing into the build
# directory, which is gone by the time anything loads it. Anything linked
# against it then fails at import with "Library not loaded". delocate hides
# this in a wheel; a plain build has to be told.
if [ "$(uname)" = "Darwin" ]; then
  for lib in "$PREFIX/lib/"libRMF*.dylib; do
    [ -f "$lib" ] || continue
    install_name_tool -id "$lib" "$lib" 2>/dev/null || true
  done
fi
fi

if [ "$WITH_IMP" = "1" ]; then
  # A minimal IMP: only the modules the connection layer names, only their
  # libraries. IMP's CMake guards every pyext directory with
  # `if(NOT IMP_STATIC)`, and building the -lib targets alone skips them
  # anyway, so no SWIG wrapper of IMP's is built and nothing links Python.
  # Measured on an M-series Mac: 57 s wall for the ten modules, 11.5 MB.
  curl -sSL "https://github.com/salilab/imp/archive/refs/tags/${IMP_VERSION}.tar.gz" | tar -xz -C "$WORK"
  IMP_SRC="$WORK/imp-${IMP_VERSION}"
  # IMP's build tooling lives in RMF's submodule, which a tarball does not
  # carry; fetch RMF's source for it if the RMF step above did not.
  if [ ! -d "${RMF_SRC:-}/tools/dev_tools" ]; then
    curl -sSL "https://github.com/salilab/rmf/archive/refs/tags/${RMF_VERSION}.tar.gz" | tar -xz -C "$WORK"
    RMF_SRC="$WORK/rmf-${RMF_VERSION}"
  fi
  rm -f "$IMP_SRC/tools/dev_tools"
  cp -R "$RMF_SRC/tools/dev_tools" "$IMP_SRC/tools/dev_tools"
  # every module that is not wanted, by name -- the tooling resolves
  # dependencies by name, so the list has to be explicit
  disabled=""
  for m in "$IMP_SRC"/modules/*/; do
    name="$(basename "$m")"
    case " $IMP_MODULES cgal " in *" $name "*) ;; *) disabled="${disabled}${name}:" ;; esac
  done
  cmake -S "$IMP_SRC" -B "$WORK/imp-build" -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_INSTALL_PREFIX="$PREFIX" -DIMP_MAX_CHECKS=NONE \
        -DIMP_DISABLED_MODULES="${disabled%:}" -DCMAKE_POSITION_INDEPENDENT_CODE=ON
  targets=""
  for m in $IMP_MODULES; do targets="$targets IMP.$m-lib"; done
  cmake --build "$WORK/imp-build" -j "$(nproc 2>/dev/null || sysctl -n hw.ncpu)" --target $targets
  cp -R "$WORK/imp-build/include/IMP" "$PREFIX/include/"
  cp -R "$WORK/imp-build/lib/"libimp_* "$PREFIX/lib/"
  # IMP reads its own data at run time (IMP::atom's top.lib and par.lib for
  # the CHARMM topology); without it every read_pdb fails with "IMP is not
  # installed or set up correctly".
  mkdir -p "$PREFIX/share/IMP"
  for m in $IMP_MODULES; do
    [ -d "$IMP_SRC/modules/$m/data" ] && cp -R "$IMP_SRC/modules/$m/data" "$PREFIX/share/IMP/$m"
  done
  echo "IMP ${IMP_VERSION} (${IMP_MODULES}) installed under $PREFIX"
fi

echo "cereal ${CEREAL_VERSION} installed under $PREFIX${WITH_RMF:+ (RMF: $WITH_RMF)}"
