#!/bin/sh
set -eu
autoreconf -fi
mkdir -p build
cd build
../configure --prefix="$PREFIX"
make -j"${CPU_COUNT:-1}"
make install
