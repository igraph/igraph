#!/bin/sh
# Creates the OS X installer package and puts it in a disk image
which hdiutil >/dev/null || ( echo "This script must be run on OS X"; exit 1 )
CWD=`pwd`
while [ ! -f setup.py ]; do cd ..; done
VERSION=`cat setup.py | grep "version =" | cut -d '=' -f 2 | tr -d "', "`
python setup.py bdist_mpkg || exit 2
MPKG="dist/python_igraph-${VERSION}-py2.5-macosx10.5.mpkg"
DMG=dist/`basename $MPKG .mpkg`.dmg
rm -f ${DMG}
hdiutil create -volname "python-igraph ${VERSION}" -layout NONE -srcfolder $MPKG $DMG
rm -rf ${MPKG}
cd $CWD
