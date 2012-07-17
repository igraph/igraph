#!/bin/bash
# Creates the OS X installer package and puts it in a disk image

FATLIB=../../fatbuild/libigraph.dylib
PYTHON_VERSIONS="2.5 2.6 2.7"

function check_universal {
	if [ `file $1 | grep -c "binary with 2 architectures"` -lt 1 ]; then
		echo "$1 is not a universal binary"
		exit 2
	fi
}

function get_dependent_libraries {
	local LIBS=`otool -L $1 | awk 'NR >= 2 { print }' | cut -f 2 | cut -d ' ' -f 1`
	echo "$LIBS"
}

function check_library_paths {
	local LIB
	for LIB in $2; do
		DIR=`dirname $LIB`
		if [ x$DIR != x/usr/lib -a x$DIR != x/usr/local/lib ]; then
			echo "$1 links to disallowed library: $LIB"
			exit 3
		fi
	done
}

function check_mandatory_library_linkage {
	local LIB
	for LIB in $2; do
		if [ x$LIB = x$3 ]; then
			return
		fi
	done
	echo "$1 does not link to required library: $3"
	exit 4
}

# Check whether we are running the script on Mac OS X
which hdiutil >/dev/null || ( echo "This script must be run on OS X"; exit 1 )

# Find the directory with setup.py
CWD=`pwd`
while [ ! -f setup.py ]; do cd ..; done

# Extract the version number from setup.py
VERSION=`cat setup.py | grep "version =" | cut -d '=' -f 2 | tr -d "', "`

# Ensure that the igraph library we are linking to is a fat binary
if [ ! -f ${FATLIB} ]; then
  pushd ../.. && tools/fatbuild.sh && popd
  if [ ! -f ${FATLIB} ]; then
    echo "Failed to build fat igraph library: ${FATLIB}"
    exit 1
  fi
fi

# Ensure that we are really linking to a fat binary and that it refers
# to libxml2 and libz
check_universal ${FATLIB}
LIBS=$(get_dependent_libraries ${FATLIB})
check_library_paths ${FATLIB} "${LIBS}"
check_mandatory_library_linkage ${FATLIB} "${LIBS}" /usr/lib/libxml2.2.dylib
check_mandatory_library_linkage ${FATLIB} "${LIBS}" /usr/lib/libz.1.dylib

# Clean up the previous build directory
rm -rf build/

# Set up ARCHFLAGS to ensure that we build a multi-arch Python extension
export ARCHFLAGS="-arch i386 -arch x86_64"

# For each Python version, build the .mpkg and the .dmg
for PYVER in $PYTHON_VERSIONS; do
  python$PYVER setup.py build_ext --no-pkg-config -I ../../include:../../build/include -L `dirname $FATLIB` || exit 3
  python$PYVER setup.py bdist_mpkg --no-pkg-config || exit 4

  # Ensure that the built library is really universal
  LIB=build/lib.macosx-*-${PYVER}/igraph/_igraph.so
  check_universal ${LIB}
  DEPS=$(get_dependent_libraries ${LIB})
  check_library_paths ${LIB} "${DEPS}"
  check_mandatory_library_linkage ${LIB} "${DEPS}" /usr/local/lib/libigraph.0.dylib

  MPKG="dist/python_igraph-${VERSION}-py${PYVER}-macosx10.5.mpkg"
  if [ ! -f $MPKG ]; then
    MPKG="dist/python_igraph-${VERSION}-py${PYVER}-macosx10.6.mpkg"
    if [ ! -f $MPKG ]; then
      MPKG="dist/python_igraph-${VERSION}-py${PYVER}-macosx10.7.mpkg"
    fi
  fi
  DMG=dist/`basename $MPKG .mpkg`.dmg
  rm -f ${DMG}
  hdiutil create -volname "python-igraph ${VERSION}" -layout NONE -srcfolder $MPKG $DMG
  rm -rf ${MPKG}
done

cd $CWD
