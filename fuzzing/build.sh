#!/bin/bash -eu

mkdir build && cd build
cmake .. -DIGRAPH_WARNINGS_AS_ERRORS=OFF
make -j$(nproc)

# Create seed corpus
zip $OUT/read_gml_fuzzer_seed_corpus.zip \
        $SRC/igraph/examples/simple/karate.gml

cd $SRC/igraph

$CXX $CXXFLAGS -I$SRC/igraph/build/include -I$SRC/igraph/include -o read_gml_fuzzer.o -c ./fuzzing/read_gml_fuzzer.cpp
$CXX $CXXFLAGS $LIB_FUZZING_ENGINE read_gml_fuzzer.o \
        -o $OUT/read_gml_fuzzer ./build/src/libigraph.a
