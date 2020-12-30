!/bin/bash -eu

sed 's/use_all_warnings/\#use_all_warnings/g' -i $SRC/igraph/src/CMakeLists.txt
mkdir build && cd build
cmake ..
make -j$(nproc)

# Create seed corpus
zip $OUT/read_lgl_fuzzer_seed_corpus.zip \
        $SRC/igraph/examples/simple/karate.gml

cd $SRC/igraph

$CC $CFLAGS -I$SRC/igraph/build/include -I$SRC/igraph/include -o read_lgl_fuzzer.o -c ./fuzzing/read_lgl_fuzzer.c
$CC $CFLAGS $LIB_FUZZING_ENGINE read_lgl_fuzzer.o \
        -o $OUT/read_lgl_fuzzer ./build/src/libigraph.a
