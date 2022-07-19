#!/bin/bash
git submodule update --recursive --init --remote
rm -r -f build
cd doc
doxygen
$PYTHON ../tools/doxy2swig.py _build/xml/index.xml ../pyext/documentation.i
cd ..

mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release -G Ninja -DCMAKE_INSTALL_LIBDIR=lib -DCMAKE_PREFIX_PATH=$PREFIX -DCMAKE_INSTALL_PREFIX=$PREFIX ..
ninja install -k 0

# Copy examples in place
cp -r ../examples/* $PREFIX/share/doc/IMP/examples/
