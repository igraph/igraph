#! /bin/sh
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
cd interfaces/R && autoconf 
cd -
cd interfaces/java && autoconf
cd -
