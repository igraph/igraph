#!/bin/sh
#
# Script that should be run whenever we bump the version number of
# igraph.
#
# This script adjusts the version numbers in the following files:
#
#   - configure.in
#   - interfaces/java/build.xml
#   - interfaces/R/configure.in
#   - examples/simple/gml.out
#   - examples/simple/cattributes2.out
#   - msvc/igraphtest/igraphtest.vcproj
#   - tools/launchpad_nightly.recipe
#   - debian/changelog

set -e
set -u

if [ $# -lt 1 ]; then
	echo "Usage: $0 version"
	exit 1
fi

VERSION="$1"

# Step to the root of the source tree
cd `dirname $0`/..

# Adjust configure.in
sed -e "s/AC_INIT(igraph, [^,]*,/AC_INIT(igraph, ${VERSION},/" \
	-e "s/AM_INIT_AUTOMAKE(igraph, [^)]*)/AM_INIT_AUTOMAKE(igraph, ${VERSION})/" \
	-i configure.in

# Adjust interfaces/java/build.xml
sed -e "s/property name=\"package\.version\" value=\"[^\"]*\"/property name=\"package.version\" value=\"${VERSION}\"/" \
	-i interfaces/java/build.xml

# Adjust interfaces/R/configure.in
sed -e "s/AC_INIT(igraph, [^,]*,/AC_INIT(igraph, ${VERSION},/" \
	-i configure.in

# Adjust examples/simple/gml.out
sed -e "s/igraph version [^ ]*/igraph version ${VERSION}/" \
	-i examples/simple/gml.out

# Adjust examples/simple/cattributes2.out
sed -e "s/igraph version [^ ]*/igraph version ${VERSION}/" \
	-i examples/simple/cattributes2.out

# Adjust msvc/igraphtest/igraphtest.vcproj
sed -e "s/igraph-[^-]*-msvc/igraph-${VERSION}-msvc/g" \
	-i msvc/igraphtest/igraphtest.vcproj

# Adjust tools/launchpad_nightly.recipe
sed -e "s/deb-version [^~]*/deb-version ${VERSION}/" \
	-e "s|lp:igraph/[^-]*-main|lp:igraph/${VERSION}-main|" \
	-i tools/launchpad_nightly.recipe

# Adjust debian/changelog
DATE="`date -R`"
cat >debian/changelog.new <<EOF
igraph (${VERSION}) unstable; urgency=low

  * Bumped version number to ${VERSION}.

 -- Tamas Nepusz <ntamas@gmail.com>  ${DATE}

EOF
cat debian/changelog >>debian/changelog.new
mv debian/changelog.new debian/changelog

# Done.
echo "Successfully bumped version number to ${VERSION}."
