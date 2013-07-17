AC_INIT(igraph, @VERSION@, csardi.gabor@gmail.com)
AC_CONFIG_SRCDIR(src/rinterface.c)
AC_CONFIG_HEADERS(src/config.h)

: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
   echo "could not determine R_HOME"
   exit 1
fi
CC=`"${R_HOME}/bin/R" CMD config CC`
CXX=`"${R_HOME}/bin/R" CMD config CXX`
FC=`"${R_HOME}/bin/R" CMD config FC`
CFLAGS=`"${R_HOME}/bin/R" CMD config CFLAGS`
CXXFLAGS=`"${R_HOME}/bin/R" CMD config CXXFLAGS`
CPPFLAGS=`"${R_HOME}/bin/R" CMD config CPPFLAGS`
FCFLAGS=`"${R_HOME}/bin/R" CMD config FCFLAGS`
FLIBS=`"${R_HOME}/bin/R" CMD config FLIBS`

AC_LANG(C)
AC_PROG_CC

# Fortran compiler, we need to check if it is the GNU compiler
AC_PROG_FC
if test "x$ac_cv_fc_compiler_gnu" == xyes; then
  AC_DEFINE([HAVE_GFORTRAN], [1], [Define to 1 if using the GNU fortran compiler])
fi


# Tricky check for C++ compiler, because Autoconf has a weird bug:
# http://lists.gnu.org/archive/html/autoconf/2006-03/msg00067.html
AC_PROG_CXX
AC_LANG_PUSH([C++])
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
#include <iostream>
const char hw[] = "Hello, World\n";]],
                [[std::cout << hw;]])],
        [AC_PROG_CXXCPP
        cxx_error=no],
        [AC_MSG_ERROR([no C++ compiler found or it cannot create executables])])
AC_LANG_POP([C++])

LIBS_SAVE=$LIBS
LIBS="$LIBS -lm"
AC_CHECK_FUNCS([rintf finite expm1 rint log2 logbl snprintf log1p round fmin stpcpy])
AC_CHECK_DECL([stpcpy],
	[AC_DEFINE([HAVE_STPCPY_SIGNATURE], [1], [Define to 1 if the stpcpy function has a signature])])
LIBS=$LIBS_SAVE

AC_CHECK_HEADER([sys/times.h], 
      [AC_DEFINE([HAVE_TIMES_H], [1], [Define to 1 if you have the sys/times.h header])])

graphml_support=yes
AC_ARG_ENABLE(graphml,
              AC_HELP_STRING([--disable-graphml], [Disable support for GraphML format]),
              [graphml_support=$enableval], [graphml_support=yes])

HAVE_LIBXML=0
if test $graphml_support = yes; then
  AC_PATH_PROG([XML2CONFIG], [xml2-config], [none])
  if test "$XML2CONFIG" = "none"; then
    graphml_support=no
  else
    XML2_LIBS=`$XML2CONFIG --libs`
    XML2_CFLAGS=`$XML2CONFIG --cflags`
    AC_CHECK_LIB([xml2], [xmlSAXUserParseFile], [
      OLDCFLAGS=${CFLAGS}
      OLDCPPFLAGS=${CPPFLAGS}
      CFLAGS=${XML2_CFLAGS}
      CPPFLAGS=${XML2_CFLAGS}
      AC_CHECK_HEADER([libxml/parser.h], [
        HAVE_LIBXML=1
	AC_DEFINE([HAVE_LIBXML], [1], [Define to 1 if you have the libxml2 libraries installed])
        CFLAGS="${OLDCFLAGS} ${XML2_CFLAGS}"
	CPPFLAGS="${OLDCFLAGS} ${XML2_CFLAGS}"
	AC_SUBST(XML2_LIBS)
	AC_SUBST(XML2_CFLAGS)
      ], [
        graphml_support=no
        CFLAGS=${OLDCFLAGS}
	CPPFLAGS=${OLDCPPFLAGS}
      ])
    ], [
      graphml_support=no
    ])
  fi
fi
AC_SUBST(HAVE_LIBXML)

AC_LANG_PUSH([C++])
HAVE_GMP=0
GMP_LIBS=""
gmp_support=no
AC_ARG_ENABLE(gmp, AC_HELP_STRING([--disable-gmp], [Compile without the GMP library]))
if test "x$enable_gmp" != "xno"; then
  AC_CHECK_LIB([gmp], [__gmpz_add], [
    AC_CHECK_HEADER([gmp.h], [
      HAVE_GMP=1
      AC_DEFINE([HAVE_GMP], [1], [Define to 1 if you have the GMP library])
      gmp_support=yes
      GMP_LIBS="-lgmp"      
    ])
  ])
fi
AC_SUBST(HAVE_GMP)
AC_SUBST(GMP_LIBS)
AC_LANG_POP([C++])

glpk_support=yes
AC_DEFINE([HAVE_GLPK], [1], [Define to 1 if you have the GLPK library])
HAVE_GLPK=1
GLPK_LIBS=""
AC_SUBST(HAVE_GLPK)
AC_SUBST(GLPK_LIBS)

AC_DEFINE(IGRAPH_THREAD_LOCAL, [], [We don't care about thread-local storage in R])

AC_CONFIG_FILES([src/Makevars])

AC_OUTPUT
