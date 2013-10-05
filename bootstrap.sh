#! /bin/sh

## Find out our version number, need git for this
printf "Finding out version number/string... "
git describe HEAD --tags | rev | sed 's/g-/./' | sed 's/-/+/' | rev > VERSION
cat VERSION

for i in glibtoolize libtoolize; do
  LIBTOOLIZE=`which $i` && break
done
if [ -z "$LIBTOOLIZE" ]; then
  echo libtoolize or glibtoolize not found or not in the path!
  exit 1
fi

set -x
aclocal
$LIBTOOLIZE --force --copy
autoheader
automake --add-missing --copy
autoconf
