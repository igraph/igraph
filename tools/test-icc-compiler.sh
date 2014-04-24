#!/bin/sh
# Test igraph compilation with Intel's C compiler

set -e

ICC_DIR=/opt/intel
source ${ICC_DIR}/bin/compilervars.sh intel64

CC=icc
CXX=icpc
LD=xild
AR=xiar
LANG=en
LANGUAGE=en
LC_ALL=C

export CC CXX LD AR LANG LANGUAGE LC_ALL

IGRAPH_ROOT=`dirname $0`/..
cd ${IGRAPH_ROOT}

rm -rf build-icc
mkdir build-icc
cd build-icc
../configure
make 2>stderr.log
cd ..