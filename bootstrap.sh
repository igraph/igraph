#! /bin/sh

cd "`dirname $0`"

## Find out our version number, need git for this
printf "Finding out version number/string... "
tools/getversion.sh > IGRAPH_VERSION
cat IGRAPH_VERSION

for i in glibtoolize libtoolize; do
  LIBTOOLIZE=`which $i` && break
done
if [ -z "$LIBTOOLIZE" ]; then
  echo libtoolize or glibtoolize not found or not in the path!
  exit 1
fi

mkdir -p m4

set -x

# autoreconf calls all its friends in the order it needs.
LIBTOOLIZE="$LIBTOOLIZE" autoreconf -vfi

# Try to patch ltmain.sh to allow -fsanitize=* linker flags to be passed
# through to the linker. Don't do anything if it fails; maybe libtool has
# been upgraded already.
patch -N -p0 -r- <tools/ltmain.patch >/dev/null || true
