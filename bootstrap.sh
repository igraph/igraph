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

# Order of commands in the next few lines are taken from here:
# https://stackoverflow.com/a/11279735/156771

$LIBTOOLIZE --force --copy

aclocal -I m4 --install
autoheader
autoconf

automake --foreign --add-missing --force-missing --copy

# Try to patch ltmain.sh to allow -fsanitize=* linker flags to be passed
# through to the linker. Don't do anything if it fails; maybe libtool has
# been upgraded already.
patch -N -p0 -r- <tools/ltmain.patch >/dev/null || true
