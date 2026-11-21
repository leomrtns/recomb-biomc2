#!/bin/sh
set -eu
export TMPDIR="$SRC_DIR/build-tmp"
mkdir -p "$TMPDIR"
autoreconf -fi
biomc2_build_dir=$(mktemp -d "$SRC_DIR/conda-build.XXXXXX")
cd "$biomc2_build_dir"
"$SRC_DIR/configure" --prefix="$PREFIX"
make -j"${CPU_COUNT:-1}"
make check
make install
