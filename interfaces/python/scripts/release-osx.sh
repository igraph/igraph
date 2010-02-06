#!/bin/sh
# Creates the OS X installer package and puts it in a disk image

FATLIB=../../fatbuild/.libs/libigraph.dylib
PYTHON_VERSIONS="2.5 2.6"

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
if [ `file ${FATLIB} | grep -c "binary with 4 architectures"` -lt 1 ]; then
	echo "${FATLIB} is not a 4-arch universal binary"
	exit 2
fi

# Clean up the previous build directory
rm -rf build/

# For each Python version, build the .mpkg and the .dmg
for PYVER in $PYTHON_VERSIONS; do
	python$PYVER setup.py build_ext -I ../../include -L `dirname $FATLIB` || exit 3
	python$PYVER setup.py bdist_mpkg || exit 4
	MPKG="dist/python_igraph-${VERSION}-py${PYVER}-macosx10.5.mpkg"
	if [ ! -f $MPKG ]; then
	  MPKG="dist/python_igraph-${VERSION}-py${PYVER}-macosx10.6.mpkg"
	fi
	DMG=dist/`basename $MPKG .mpkg`.dmg
	rm -f ${DMG}
	hdiutil create -volname "python-igraph ${VERSION}" -layout NONE -srcfolder $MPKG $DMG
	rm -rf ${MPKG}
done

cd $CWD
