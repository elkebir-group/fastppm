#!/usr/bin/env bash

mkdir -p build
cd build

cmake -DPYTHON_EXECUTABLE=$PYTHON -DCMAKE_BUILD_TYPE=RELEASE ..

make -j${CPU_COUNT}

cp src/fastppm-cli $PREFIX/bin
cp src/fastppm*.so $SP_DIR

