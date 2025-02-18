#!/bin/bash
set -eu

#Some cleanup
rm -rf install

#Download development version
git clone https://github.com/SciCompKL/adjoint_petsc.git src
mkdir install
mkdir build

cmake src \
    -B build \
    -DBUILD_TESTING=OFF \
    -DCMAKE_INSTALL_PREFIX=$ISSM_DIR/externalpackages/adjointpetsc/install \
    -DCoDiPack_DIR=$ISSM_DIR/externalpackages/codipack/install/cmake \
    -DPESTc_DIR=$ISSM_DIR/externalpackages/petsc/install
if [ $# -eq 0 ]; then
  cmake --build build
else
  cmake --build build -j $1
fi
cmake --install build
