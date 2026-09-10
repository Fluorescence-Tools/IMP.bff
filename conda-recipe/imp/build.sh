#!/bin/bash

# Help CMake to find CGAL

# Make sure the default encoding for files opened by Python 3 is UTF8
export LANG=en_US.UTF-8

# Don't build the scratch or cnmultifit modules
DISABLED=EMageFit:bayesianem:bff:cgal:cnmultifit:domino:em2d:emseqfinder:example:foxs:gsl:integrative_docking:kmeans:misc:modeller:mpi:multi_state:multifit:nestor:npc:npctransport:parallel:pepdock:pmi1:sampcon:saxs_merge:scratch:spatiotemporal:spb:symmetry:test

# Avoid running out of memory on by splitting up IMP.cgal and IMP.spb

# Force C++17 compilation to build successfully with newer protobuf
if [ `uname -s` = "Darwin" ]; then
  CXX_FLAGS="-std=c++17"
fi

mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release -DIMP_DISABLED_MODULES=${DISABLED} \
      -G Ninja \
      -DIMP_USE_SYSTEM_RMF=on \
      -DIMP_USE_SYSTEM_IHM=on \
      ${CMAKE_ARGS} \
      -DCMAKE_CXX_FLAGS="${CXX_FLAGS}" \
      -DPython3_FIND_FRAMEWORK=NEVER \
      ..

# Make sure all modules we asked for were found (this is tested for
# in the final package, but quicker to abort here if they're missing)
python "${RECIPE_DIR}/check_disabled_modules.py" ${DISABLED} || exit 1

if [ `uname -s` = "Darwin" ]; then
  ninja install
else
  ninja install -j 4
fi

# Activation scripts: our kernel's compiled-in data path does not survive
# rattler-build's prefix relocation (the search list arrives empty at run
# time; conda-build builds don't hit this), so IMP_DATA -- the first thing
# IMP::internal::get_data_prefixes checks -- is exported on activation.
mkdir -p ${PREFIX}/etc/conda/activate.d ${PREFIX}/etc/conda/deactivate.d
echo "export IMP_DATA=\$CONDA_PREFIX/share/IMP" > ${PREFIX}/etc/conda/activate.d/imp.sh
echo "unset IMP_DATA" > ${PREFIX}/etc/conda/deactivate.d/imp.sh

# Don't distribute example application
rm -f ${PREFIX}/bin/imp_example_app
