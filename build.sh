#!/usr/bin/env bash

mkdir -p build
cd build

# Pass the Conda version variable ($PKG_VERSION) to CMake
cmake -DPYTHON_EXECUTABLE=$PYTHON \
      -DCMAKE_BUILD_TYPE=RELEASE \
      -DFASTPPM_VERSION=$PKG_VERSION \
      ..

make -j${CPU_COUNT}

# Make sure target directories exist before copying
mkdir -p $PREFIX/bin
mkdir -p $SP_DIR

# Install the binary and the python module
cp src/fastppm-cli $PREFIX/bin
cp src/fastppm*.so $SP_DIR