FROM gcr.io/oss-fuzz-base/base-builder:v1
RUN apt-get update && apt-get install -y cmake bison flex
COPY . $SRC/igraph
WORKDIR igraph
COPY fuzzing/build.sh $SRC/
