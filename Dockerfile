# Container image with the igraph C library installed under /usr/local,
# ready for compiling and linking programs against igraph.
#
# Build the image from the root of the source tree:
#
#   docker build -t igraph .
#
# The version number is taken from Git, so the build context includes .git.
# In a shallow clone without tags, pass it explicitly instead:
#
#   docker build --build-arg IGRAPH_VERSION=1.0.1 -t igraph .
#
# Run the test suite by building the "test" stage:
#
#   docker build --target test .
#
# Compile a program inside the image, e.g. with pkg-config:
#
#   docker run --rm -v "$PWD":/work -w /work igraph \
#     sh -c 'cc prog.c $(pkg-config --cflags --libs igraph) -o prog && ./prog'

ARG UBUNTU_VERSION=24.04

# Toolchain and the external libraries igraph links to.
# The -dev packages stay in the final image so that both shared and static
# builds of igraph can be linked against.
FROM ubuntu:${UBUNTU_VERSION} AS base
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        build-essential \
        cmake \
        ninja-build \
        pkg-config \
        libarpack2-dev \
        libblas-dev \
        liblapack-dev \
        libglpk-dev \
        libgmp-dev \
        libplfit-dev \
        libxml2-dev \
        # Needed by `pkg-config --static` for libxml2, but not pulled in by libxml2-dev
        liblzma-dev \
    && rm -rf /var/lib/apt/lists/*

FROM base AS build
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        bison \
        flex \
        git \
    && rm -rf /var/lib/apt/lists/*

ARG IGRAPH_VERSION=
ARG CMAKE_BUILD_TYPE=Release
ARG BUILD_SHARED_LIBS=ON

COPY . /src
WORKDIR /src
RUN if [ -n "$IGRAPH_VERSION" ]; then echo "$IGRAPH_VERSION" > IGRAPH_VERSION; fi \
    && cmake -S . -B build -G Ninja \
        -DCMAKE_BUILD_TYPE="$CMAKE_BUILD_TYPE" \
        -DBUILD_SHARED_LIBS="$BUILD_SHARED_LIBS" \
        -DCMAKE_INSTALL_PREFIX=/usr/local \
    && cmake --build build \
    && DESTDIR=/staging cmake --install build

FROM build AS test
RUN cmake --build build --target check

FROM base
COPY --from=build /staging/ /
RUN ldconfig
WORKDIR /work
