#! /bin/sh

set -x
aclocal
libtoolize --force --copy
gtkdocize || exit 0
autoheader
automake --add-missing --copy
autoconf
