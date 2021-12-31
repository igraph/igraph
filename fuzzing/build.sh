#!/bin/bash -eu

mkdir build && cd build
cmake .. -DIGRAPH_WARNINGS_AS_ERRORS=OFF
make -j$(nproc)

# Create seed corpus
zip $OUT/read_gml_fuzzer_seed_corpus.zip \
        $SRC/igraph/examples/simple/karate.gml

cd $SRC/igraph

for TARGET in read_gml_fuzzer bliss_fuzzer vertex_connectivity_fuzzer edge_connectivity_fuzzer vertex_separators_fuzzer
do
  $CXX $CXXFLAGS -I$SRC/igraph/build/include -I$SRC/igraph/include -o $TARGET.o -c ./fuzzing/$TARGET.cpp
  $CXX $CXXFLAGS $LIB_FUZZING_ENGINE $TARGET.o -o $OUT/$TARGET ./build/src/libigraph.a
done
