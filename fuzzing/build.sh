#!/bin/bash -eu

export DEPS_PATH=/src/deps
mkdir $DEPS_PATH

# Build libxml2 without ICU support, https://github.com/igraph/igraph/issues/1992
# Building libxml2 from scratch is also required for MemorySanitizer.
# It may be necessary to leave CMAKE_BUILD_TYPE empty and specify LIBXML2_WITH_MODULES=OFF
# in order for fuzz introspector builds to succeed (details unverified).
cd $SRC
wget https://download.gnome.org/sources/libxml2/2.12/libxml2-2.12.6.tar.xz
tar xf libxml2-2.12.6.tar.xz
cd libxml2-2.12.6
mkdir build && cd build
CFLAGS="$CFLAGS -O2" cmake .. -DCMAKE_INSTALL_PREFIX=$DEPS_PATH -DBUILD_SHARED_LIBS=OFF -DLIBXML2_WITH_ICU=OFF -DLIBXML2_WITH_PYTHON=OFF -DLIBXML2_WITH_TESTS=OFF -DLIBXML2_WITH_ZLIB=OFF -DLIBXML2_WITH_LZMA=OFF -DLIBXML2_WITH_PROGRAMS=OFF -DLIBXML2_WITH_MODULES=OFF
make install -j$(nproc)

# Build igraph
cd $SRC/igraph
mkdir build && cd build
# CMAKE_BUILD_TYPE=None is an arbitrary value that prevents the automatic Release
# build type setting, allowing OSS-Fuzz to pass on its own optimization flags.
cmake .. -DIGRAPH_WARNINGS_AS_ERRORS=OFF -DCMAKE_BUILD_TYPE=None -DCMAKE_PREFIX_PATH=$DEPS_PATH -DFLEX_KEEP_LINE_NUMBERS=ON
make -j$(nproc)


# Create seed corpus
zip $OUT/read_edgelist_seed_corpus.zip \
        $SRC/igraph/fuzzing/test_inputs/*.el

zip $OUT/read_dimacs_flow_seed_corpus.zip \
        $SRC/igraph/examples/simple/*.max \
        $SRC/igraph/tests/unit/*.max \
        $SRC/igraph/fuzzing/test_inputs/*.max

zip $OUT/read_dl_seed_corpus.zip \
        $SRC/igraph/examples/simple/*.dl \
        $SRC/igraph/tests/unit/*.dl \
        $SRC/igraph/fuzzing/test_inputs/*.dl

zip $OUT/read_gml_seed_corpus.zip \
        $SRC/igraph/examples/simple/*.gml \
        $SRC/igraph/tests/regression/*.gml \
        $SRC/igraph/tests/unit/*.gml \
        $SRC/igraph/fuzzing/test_inputs/*.gml

zip $OUT/read_graphdb_seed_corpus.zip \
        $SRC/igraph/fuzzing/test_inputs/*.A?? \
        $SRC/igraph/fuzzing/test_inputs/*.B?? \
        $SRC/igraph/examples/simple/*.A??

zip $OUT/read_graphml_seed_corpus.zip \
        $SRC/igraph/examples/simple/*.graphml \
        $SRC/igraph/tests/unit/*.graphml \
        $SRC/igraph/tests/regression/*.graphml \
        $SRC/igraph/fuzzing/test_inputs/*.graphml

zip $OUT/read_lgl_seed_corpus.zip \
        $SRC/igraph/examples/simple/*.lgl \
        $SRC/igraph/tests/unit/*.lgl \
        $SRC/igraph/fuzzing/test_inputs/*.lgl

zip $OUT/read_ncol_seed_corpus.zip \
        $SRC/igraph/examples/simple/*.ncol \
        $SRC/igraph/tests/unit/*.ncol \
        $SRC/igraph/fuzzing/test_inputs/*.ncol

zip $OUT/read_pajek_seed_corpus.zip \
        $SRC/igraph/examples/simple/links.net \
        $SRC/igraph/tests/unit/bipartite.net \
        $SRC/igraph/tests/unit/pajek*.net \
        $SRC/igraph/tests/regression/*.net \
        $SRC/igraph/fuzzing/test_inputs/*.net

zip $OUT/write_all_gml_seed_corpus.zip \
        $SRC/igraph/examples/simple/*.gml \
        $SRC/igraph/tests/regression/*.gml \
        $SRC/igraph/tests/unit/*.gml \
        $SRC/igraph/fuzzing/test_inputs/*.gml

zip $OUT/write_all_graphml_seed_corpus.zip \
        $SRC/igraph/examples/simple/*.graphml \
        $SRC/igraph/tests/unit/*.graphml \
        $SRC/igraph/tests/regression/*.graphml \
        $SRC/igraph/fuzzing/test_inputs/*.graphml

cd $SRC/igraph

XML2_FLAGS=`$DEPS_PATH/bin/xml2-config --cflags --libs`

# disabled:
#  - nothing at the moment
# disabled for UBSan:
#  - read_dimacs_flow, needs a complete rewrite, see https://github.com/igraph/igraph/issues/1924
#  - write_all_gml|graphml, uses `(igraph_integer_t) x == x` to check representability as integer; this triggers UBSan
TARGETS="read_edgelist read_dl read_gml read_graphdb read_graphml read_lgl read_ncol read_pajek bliss edge_connectivity vertex_connectivity vertex_separators basic_properties_directed basic_properties_undirected linear_algos_directed linear_algos_undirected centrality community weighted_centrality weighted_community misc_algos"
if [ "$SANITIZER" != "undefined" ]
then
  TARGETS="$TARGETS read_dimacs_flow write_all_gml write_all_graphml"
fi

for TARGET in $TARGETS
do
  if [ -e ./fuzzing/$TARGET.dict ]
  then
    cp ./fuzzing/$TARGET.dict $OUT
  fi
  $CXX $CXXFLAGS -I$SRC/igraph/build/include -I$SRC/igraph/include -o $TARGET.o -c ./fuzzing/$TARGET.cpp
  $CXX $CXXFLAGS $LIB_FUZZING_ENGINE $TARGET.o -o $OUT/$TARGET ./build/src/libigraph.a $XML2_FLAGS
done
