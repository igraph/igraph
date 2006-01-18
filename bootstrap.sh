#! /bin/sh

set -x
aclocal
libtoolize --force --copy
autoheader
automake --add-missing --copy
autoconf
cd interfaces/R/igraph && autoconf 
cd -
