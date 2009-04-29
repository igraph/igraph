#!/bin/bash

IGRAPH_VERSION=0.5.2
DEST_DIR=$HOME/packages

function install_build_dependencies {
  apt-get -y --no-install-recommends install \
          debhelper automake autoconf gcc g++ pkg-config \
          build-essential libxml2-dev python-all-dev \
          python-central python-epydoc texlive-latex-base
}

function make_destdir {
  mkdir -p ${DEST_DIR}
}

function create_igraph_debian_pkg {
  wget http://switch.dl.sourceforge.net/sourceforge/igraph/igraph-$1.tar.gz
  tar -xvvzf igraph-$1.tar.gz
  cd igraph-$1
  # Patching igraph.pc.in if needed
  patch -p1 << EOF
--- old/igraph.pc.in    2009-04-29 11:36:21.000000000 +0100
+++ new/igraph.pc.in    2009-04-29 11:36:21.000000000 +0100
@@ -9,2 +9,2 @@
 Libs: -L\${libdir} -ligraph
-Cflags: -I\${includedir}/igraph
+Cflags: -I\${includedir}
EOF
  debian/prepare
  dpkg-buildpackage -tc -rfakeroot
  cd ..
}

function install_igraph_debian_pkg {
  dpkg -i libigraph_$1_*.deb libigraph-dev_$1_*.deb
}

function remove_igraph_debian_pkg {
  dpkg -r libigraph libigraph-dev
}

function compile_python_interface {
  wget http://pypi.python.org/packages/source/p/python-igraph/python-igraph-$1.tar.gz
  tar -xvvzf python-igraph-$1.tar.gz
  cd python-igraph-$1
  debian/prepare
  dpkg-buildpackage -tc -rfakeroot
  cd ..
}

function move_packages {
  mv igraph_* libigraph* python-igraph_* python-igraph-doc_* ${DEST_DIR}
}

function remove_build_dir {
  cd /
  rm -rf ${BUILD_DIR}
}

if [ `id -u` != 0 ]; then
  echo This script must be run as root.
  exit 1
fi

BUILD_DIR=`mktemp -d`
cd ${BUILD_DIR}

install_build_dependencies
make_destdir
create_igraph_debian_pkg ${IGRAPH_VERSION}
install_igraph_debian_pkg ${IGRAPH_VERSION}
compile_python_interface ${IGRAPH_VERSION}
remove_igraph_debian_pkg
move_packages
# remove_build_dir
