#!/bin/bash
git submodule update --recursive --init --remote
cd doc
mkdir -p _build/html/stable/api
doxygen
$PYTHON ../tools/doxy2swig.py _build/xml/index.xml ../pyext/documentation.i
cd ..

mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release -G Ninja ..
ninja install -k 0

# Copy examples
mkdir $PREFIX/share/doc/IMP/examples/bff
cp -r ../examples/* $PREFIX/share/doc/IMP/examples/bff

