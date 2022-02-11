#!/bin/bash -eu

mkdir build && cd build
# CMAKE_BUILD_TYPE=None is an arbitrary value that prevents the automatic Release
# build type setting, allowing OSS-Fuzz to pass on its own optimization flags.
cmake .. -DIGRAPH_WARNINGS_AS_ERRORS=OFF -DCMAKE_BUILD_TYPE=None
make -j$(nproc)

# Create seed corpus
zip $OUT/read_gml_fuzzer_seed_corpus.zip \
        $SRC/igraph/examples/simple/*.gml \
        $SRC/igraph/tests/regression/*.gml \
        $SRC/igraph/fuzzing/test_inputs/*.gml

zip $OUT/read_pajek_fuzzer_seed_corpus.zip \
        $SRC/igraph/examples/simple/links.net \
        $SRC/igraph/tests/unit/bipartite.net \
        $SRC/igraph/tests/unit/pajek*.net \
        $SRC/igraph/tests/regression/*.net \
        $SRC/igraph/fuzzing/test_inputs/*.net

zip $OUT/read_dl_fuzzer_seed_corpus.zip \
        $SRC/igraph/examples/simple/*.dl \
        $SRC/igraph/tests/unit/*.dl \
        $SRC/igraph/fuzzing/test_inputs/*.dl

zip $OUT/read_lgl_fuzzer_seed_corpus.zip \
        $SRC/igraph/examples/simple/*.lgl \
        $SRC/igraph/tests/unit/*.lgl \
        $SRC/igraph/fuzzing/test_inputs/*.lgl

zip $OUT/read_ncol_fuzzer_seed_corpus.zip \
        $SRC/igraph/examples/simple/*.ncol \
        $SRC/igraph/tests/unit/*.ncol \
        $SRC/igraph/fuzzing/test_inputs/*.ncol

zip $OUT/read_graphml_fuzzer_seed_corpus.zip \
        $SRC/igraph/examples/simple/*.graphml \
        $SRC/igraph/tests/unit/*.graphml \
        $SRC/igraph/tests/regression/*.graphml \
        $SRC/igraph/fuzzing/test_inputs/*.graphml

cd $SRC/igraph

for TARGET in read_gml_fuzzer read_pajek_fuzzer read_dl_fuzzer read_lgl_fuzzer read_ncol_fuzzer bliss_fuzzer vertex_connectivity_fuzzer edge_connectivity_fuzzer vertex_separators_fuzzer
do
  $CXX $CXXFLAGS -I$SRC/igraph/build/include -I$SRC/igraph/include -o $TARGET.o -c ./fuzzing/$TARGET.cpp
  $CXX $CXXFLAGS $LIB_FUZZING_ENGINE $TARGET.o -o $OUT/$TARGET ./build/src/libigraph.a
done
