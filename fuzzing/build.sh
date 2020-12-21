#!/bin/bash -eu

./bootstrap.sh
./configure
make -j$(nproc)

# Create seed corpus
zip $OUT/read_lgl_fuzzer_seed_corpus.zip \
	$SRC/igraph/examples/simple/karate.gml


mv $SRC/igraph/fuzzing/read_lgl_fuzzer.c .
$CC $CFLAGS -I./include -o read_lgl_fuzzer.o -c read_lgl_fuzzer.c
$CC $CFLAGS $LIB_FUZZING_ENGINE read_lgl_fuzzer.o \
	-o $OUT/read_lgl_fuzzer /src/igraph/src/.libs/libigraph.a
