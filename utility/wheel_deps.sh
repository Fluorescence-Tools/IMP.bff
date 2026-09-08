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
CEREAL_VERSION="${CEREAL_VERSION:-1.3.2}"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT
mkdir -p "$PREFIX/include" "$PREFIX/lib"

# cereal: headers only
curl -sSL "https://github.com/USCiLab/cereal/archive/refs/tags/v${CEREAL_VERSION}.tar.gz" | tar -xz -C "$WORK"
cp -R "$WORK/cereal-${CEREAL_VERSION}/include/cereal" "$PREFIX/include/"

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
echo "cereal ${CEREAL_VERSION} and RMF ${RMF_VERSION} installed under $PREFIX"
