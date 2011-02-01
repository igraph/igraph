#!/bin/bash

if [ `id -u` != 0 ]; then
	echo "This script must be run as root."
	exit 2
fi

if [ $# -ne 3 ]; then
    echo "Usage: $0 igraph_version debian_revision series"
	echo "where igraph_version is the version number of igraph,"
	echo "debian_revision is the Debian package revision number"
	echo "and series is the distribution series (e.g., unstable,"
	echo "karmic, jaunty etc)."
	echo ""
	echo "For the Launchpad PPA, you have to prepare a package"
	echo "for the two most recent Ubuntu series"
	exit 1
fi

IGRAPH_VERSION=$1
DEBIAN_REVISION=$2
SERIES=$3

DEST_DIR=$HOME/packages/$SERIES
CURR_DIR=`pwd`
BZR_IGRAPH_ROOT=$HOME/bzr/igraph

function prechecks {
  if [ ! -d ${BZR_IGRAPH_ROOT} ]; then
    echo ${BZR_IGRAPH_ROOT} does not exist or is not a directory!
    exit 1
  fi
}

function install_build_dependencies {
  apt-get update
  apt-get -y --no-install-recommends install \
          debhelper devscripts fakeroot automake autoconf gcc g++ \
          pkg-config flex bison build-essential libxml2-dev libglpk-dev \
          libarpack2-dev libgmp3-dev libxml2-dev libblas-dev liblapack-dev \
          python-all-dev python-central python-epydoc texlive-latex-base 
}

function make_destdir {
  mkdir -p ${DEST_DIR}
}

function bazaar_update {
  pushd ${BZR_IGRAPH_ROOT}/0.6-main
  bzr update
  popd
}

function create_igraph_debian_pkg {
  if [ -f ${CURR_DIR}/igraph_$1.orig.tar.gz ]; then
    cp ${CURR_DIR}/igraph_$1.orig.tar.gz .
  else
    wget -O igraph_$1.orig.tar.gz http://switch.dl.sourceforge.net/sourceforge/igraph/igraph-$1.tar.gz
  fi
  tar -xvvzf igraph_$1.orig.tar.gz
  cd igraph-$1
  cp -r ${BZR_IGRAPH_ROOT}/0.6-main/debian .
  cat debian/changelog.in | sed -e "s/@VERSION@/@VERSION@-${DEBIAN_REVISION}/g" >debian/changelog.in.new
  mv debian/changelog.in.new debian/changelog.in
  debian/prepare
  cat debian/changelog | sed -e "s/unstable/${SERIES}/g" >debian/changelog.new
  mv debian/changelog.new debian/changelog
  # dpkg-buildpackage -tc -rfakeroot
  debuild -b -us -uc
  debuild -S -sa
  cd ..
}

function install_igraph_debian_pkg {
  dpkg -i libigraph0_$1-$2_*.deb libigraph-dev_$1-$2_*.deb || exit 3
}

function remove_igraph_debian_pkg {
  dpkg -r libigraph0 libigraph-dev
}

function compile_python_interface {
  wget -O python-igraph_$1.orig.tar.gz http://pypi.python.org/packages/source/p/python-igraph/python-igraph-$1.tar.gz
  tar -xvvzf python-igraph_$1.orig.tar.gz
  cd python-igraph-$1
  cat debian/changelog.in | sed -e "s/@VERSION@/@VERSION@-${DEBIAN_REVISION}/g" >debian/changelog.in.new
  mv debian/changelog.in.new debian/changelog.in
  debian/prepare
  cat debian/changelog | sed -e "s/unstable/${SERIES}/g" >debian/changelog.new
  mv debian/changelog.new debian/changelog
  debuild -b -us -uc
  # debuild -S -sa
  cd ..
}

function move_packages {
  mv igraph_* libigraph* python-igraph_* python-igraph-doc_* ${DEST_DIR}
}

function remove_build_dir {
  cd /
  rm -rf ${BUILD_DIR}
}

prechecks

BUILD_DIR=`mktemp -d`
cd ${BUILD_DIR}
# trap "{ cd /; rm -rf \"${BUILD_DIR}\"; exit 255; }" INT EXIT TERM

install_build_dependencies
make_destdir
bazaar_update
create_igraph_debian_pkg ${IGRAPH_VERSION}
install_igraph_debian_pkg ${IGRAPH_VERSION} ${DEBIAN_REVISION}
compile_python_interface ${IGRAPH_VERSION}
remove_igraph_debian_pkg
move_packages
remove_build_dir

echo ""
echo "Commands to upload packages to the Launchpad PPA:"
echo "dput ppa:igraph/ppa <source.changes>"
