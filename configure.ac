AC_INIT(igraph, esyscmd([tr -d '\n' < IGRAPH_VERSION]), igraph@igraph.org)
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_SRCDIR(src/games.c)
AM_INIT_AUTOMAKE([foreign subdir-objects])
AC_CONFIG_HEADERS([config.h])

m4_include(tools/autoconf/ax_tls.m4)
m4_include(tools/autoconf/as-version.m4)

AS_VERSION

# Define list of additional libraries that have to be linked to igraph when
# another app tries to link to the static library of igraph. This is substituted
# into igraph.pc later on.
PKGCONFIG_LIBS_PRIVATE="-lxml2 -lz -lm"
AC_SUBST(PKGCONFIG_LIBS_PRIVATE)

# Test suite
AC_CONFIG_TESTDIR(tests)
AC_CONFIG_FILES([tests/Makefile tests/atlocal])

# Don't allow AC_PROG_CC to set a default CFLAGS or CXXFLAGS
: ${CFLAGS=""}
: ${CXXFLAGS=""}

AC_LANG(C)
AC_PROG_CC

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

AM_PROG_LEX
AC_PROG_YACC

AC_CHECK_HEADER([sys/times.h],
      [AC_DEFINE([HAVE_TIMES_H], [1], [Define to 1 if you have the sys/times.h header])],
      [CPPFLAGS="$CPPFLAGS -DMSDOS"],
      )

AC_LIBTOOL_WIN32_DLL
AC_LIBTOOL_DLOPEN
AC_PROG_LIBTOOL
AM_MISSING_PROG([AUTOM4TE], [autom4te])

AC_HEADER_STDC
AC_CHECK_HEADERS([stdarg.h stdlib.h string.h time.h unistd.h stdint.h sys/int_types.h])
LIBS_SAVE=$LIBS
LIBS="$LIBS -lm"
AC_CHECK_FUNCS([expm1 rint rintf finite log2 snprintf log1p round fabsl fmin strcasecmp isfinite isnan strdup _strdup ftruncate stpcpy])
AC_CHECK_DECL([stpcpy],
	[AC_DEFINE([HAVE_STPCPY_SIGNATURE], [1], [Define to 1 if the stpcpy function has a signature])])
LIBS=$LIBS_SAVE

AC_DEFUN([IGRAPH_WARNING],
[AC_MSG_CHECKING(whether compiler accepts $1)
AC_SUBST(WARNING_CFLAGS)
ac_save_CFLAGS="$CFLAGS"
CFLAGS="$CFLAGS $1"
AC_TRY_COMPILE(,
[int x;],
WARNING_CFLAGS="$WARNING_CFLAGS $1"
AC_MSG_RESULT(yes),
AC_MSG_RESULT(no))
CFLAGS="$ac_save_CFLAGS"])

AC_DEFUN([IGRAPH_CC_SWITCH],
[AC_MSG_CHECKING(whether compiler supports $1)
ac_save_CFLAGS="$CFLAGS"
CFLAGS="$CFLAGS $1"
AC_TRY_COMPILE(,
[int x;],
AC_MSG_RESULT(yes)
$2,
AC_MSG_RESULT(no)
$3)
CFLAGS="$ac_save_CFLAGS"])

## Solaris cc does not support -ffloat-store.
FLOATSTORE=
IGRAPH_CC_SWITCH([-ffloat-store], [FLOATSTORE="-ffloat-store"], [])
IGRAPH_CC_SWITCH([-fstore], [FLOATSTORE="$FLOATSTORE -fstore"], [])
AC_SUBST(FLOATSTORE)

AC_ARG_ENABLE(gcc-warnings,
              AC_HELP_STRING([--enable-gcc-warnings],
	                     [turn on lots of GCC warnings (not recommended)]),
[case "${enableval}" in
   yes|no) ;;
   *)      AC_MSG_ERROR([bad value ${enableval} for gcc-warnings option]) ;;
 esac],
              [enableval=no])
if test "${enableval}" = yes; then
  IGRAPH_WARNING(-Werror)
  AC_SUBST([WERROR_CFLAGS], [$WARNING_CFLAGS])
  WARNING_CFLAGS=
  IGRAPH_WARNING(-Wall)
  IGRAPH_WARNING(-W)
  IGRAPH_WARNING(-Wbad-function-cast)
  IGRAPH_WARNING(-Wcast-align)
  IGRAPH_WARNING(-Wcast-qual)
  IGRAPH_WARNING(-Wformat)
  IGRAPH_WARNING(-Wmissing-declarations)
  IGRAPH_WARNING(-Wmissing-prototypes)
  IGRAPH_WARNING(-Wnested-externs)
  IGRAPH_WARNING(-Wshadow)
  IGRAPH_WARNING(-Wstrict-prototypes)
  IGRAPH_WARNING(-Wwrite-strings)
else
  WARNING_CFLAGS=
  IGRAPH_WARNING(-Wall)
fi

use_gprof=no
AC_ARG_ENABLE(profiling,
              AC_HELP_STRING([--enable-profiling], [Enable gprof profiling]),
              [use_gprof=$enableval], [use_gprof=no])

use_asan=no
AC_ARG_ENABLE(asan,
			  AC_HELP_STRING([--enable-asan], [Enable Clang address sanitizer]),
			  [use_asan=$enableval], [use_asan=no])

debug=no
AC_ARG_ENABLE(debug,
              AC_HELP_STRING([--enable-debug], [Enable debug build]),
	      [debug=$enableval])

graphml_support=yes
AC_ARG_ENABLE(graphml,
              AC_HELP_STRING([--disable-graphml], [Disable support for GraphML format]),
              [graphml_support=$enableval], [graphml_support=yes])

HAVE_LIBXML=0
if test $graphml_support = yes; then
  AC_PATH_PROG([XML2CONFIG], [xml2-config], [none])
  if test "$XML2CONFIG" = "none"; then
	# Hmmm, no xml2-config. Older versions of OS X do not have it while still
	# having libxml2, so let's try an educated guess if
	# /usr/include/libxml2/libxml/parser.h exists.
	AC_CHECK_FILE([/usr/include/libxml2/libxml/parser.h], [
	  XML2_LIBS="-lxml2 -lz -lm"
	  XML2_CFLAGS="-I/usr/include/libxml2"
	], [
	  graphml_support=no
	])
  else
    XML2_LIBS=`$XML2CONFIG --libs`
    XML2_CFLAGS=`$XML2CONFIG --cflags`
	if test -f /usr/bin/sw_vers -a `sw_vers -productVersion | grep -c "^10\.11\."` -gt 0 -a `echo $XML2_LIBS | grep -c '/Developer/'` -gt 0; then
      # We are on OS X 10.11, which has a known bug with xml2-config; see
      # https://github.com/igraph/igraph/issues/973 . We work around it
      # here
      XML2_LIBS=`$XML2CONFIG --exec-prefix=/usr --libs`
	fi
  fi
  AC_CHECK_LIB([xml2], [xmlSAXUserParseFile], [
    ac_save_CFLAGS="$CFLAGS"
    ac_save_CPPFLAGS="$CPPFLAGS"
    CFLAGS=${XML2_CFLAGS}
    CPPFLAGS=${XML2_CFLAGS}
    AC_CHECK_HEADER([libxml/parser.h], [
      HAVE_LIBXML=1
      AC_DEFINE([HAVE_LIBXML], [1], [Define to 1 if you have the libxml2 libraries installed])
	  CFLAGS="$ac_save_CFLAGS ${XML2_CFLAGS}"
      CPPFLAGS="$ac_save_CPPFLAGS"
      AC_SUBST(XML2_LIBS)
      AC_SUBST(XML2_CFLAGS)
    ], [
      graphml_support=no
      CFLAGS="$ac_save_CFLAGS"
	  CPPFLAGS="$ac_save_CPPFLAGS"
    ])
  ], [
    graphml_support=no
  ])
fi

AC_LANG_PUSH([C++])
gmp_support=no
AC_ARG_ENABLE(gmp, AC_HELP_STRING([--disable-gmp], [Compile without the GMP library]))
if test "x$enable_gmp" != "xno"; then
  AC_CHECK_LIB([gmp], [__gmpz_add], [
    AC_CHECK_HEADER([gmp.h], [
      AC_DEFINE([HAVE_GMP], [1], [Define to 1 if you have the GMP library])
      gmp_support=yes
      LDFLAGS="${LDFLAGS} -lgmp"
      PKGCONFIG_LIBS_PRIVATE="${PKGCONFIG_LIBS_PRIVATE} -lgmp"
    ])
  ])
fi
AC_LANG_POP([C++])

tls_support=no
HAVE_TLS=0
THREAD_LOCAL=
AC_ARG_ENABLE(tls, AC_HELP_STRING([--enable-tls], [Compile with thread-local storage]))
if test "x$enable_tls" = "xyes"; then
  keywords="__thread __declspec(thread)"
  for kw in $keywords ; do
      AC_TRY_COMPILE([int $kw test;], [], ac_cv_tls=$kw)
      AC_TRY_COMPILE([int $kw test;], [], ac_cv_tls=$kw ; break ;)
  done
  AX_TLS([
    AC_DEFINE([HAVE_TLS], [1], [Define to 1 if you want to use thread-local storage for global igraph structures])
    tls_support=yes
    HAVE_TLS=1
    THREAD_LOCAL=$ac_cv_tls
  ], [])
fi
AC_SUBST(HAVE_TLS)
AC_DEFINE_UNQUOTED([IGRAPH_THREAD_LOCAL], $THREAD_LOCAL,
          [Keyword for thread local storage, or empty if not available])
AC_DEFINE_UNQUOTED([IGRAPH_F77_SAVE], [static IGRAPH_THREAD_LOCAL],
          [Keyword for thread local storage, or just static if not available])

AC_ARG_WITH([external-f2c], [AS_HELP_STRING([--with-external-f2c],
		                      [Use external F2C library [default=no]])],
            [internal_f2c=no],
	    [internal_f2c=yes])
AC_ARG_WITH([external-blas], [AS_HELP_STRING([--with-external-blas],
		                      [Use external BLAS library [default=no]])],
            [internal_blas=no],
	    [internal_blas=yes])
AC_ARG_WITH([external-lapack], [AS_HELP_STRING([--with-external-lapack],
		                      [Use external LAPACK library [default=no]])],
            [internal_lapack=no],
	    [internal_lapack=yes])
AC_ARG_WITH([external-arpack], [AS_HELP_STRING([--with-external-arpack],
		                      [Use external ARPACK library [default=no]])],
            [internal_arpack=no],
	    [internal_arpack=yes])

AC_ARG_WITH([external-glpk], [AS_HELP_STRING([--with-external-glpk],
				      [Use external GLPK library [default=no]])],
            [internal_glpk=no],
            [internal_glpk=yes])

needs_f2c="no"
if test "$internal_blas" = "yes" -o "$internal_lapack" = "yes" -o "$internal_arpack" = "yes"; then
  needs_f2c="yes"
fi

if test "$needs_f2c" = "yes"; then
  if test "$internal_f2c" = "no"; then
    AC_CHECK_LIB([f2c], [f77_alloc_], [],
      AC_CHECK_LIB([f2c], [f77_alloc], [],
        AC_CHECK_LIB([f2c], [F77_ALLOC_], [],
           AC_CHECK_LIB([f2c], [F77_ALLOC], [],
	      [AC_MSG_RESULT(not found, trying to use -lf2c anyway.)]))))
    LDFLAGS="${LDFLAGS} -lf2c"
  else
    AC_DEFINE([INTERNAL_F2C], [1], [Define to 1 if you use the internal F2C library])
  fi
else
  internal_f2c=no
fi

if test "$internal_blas" = "no"; then
  AC_CHECK_LIB([blas], [daxpy_], [],
     AC_CHECK_LIB([blas], [daxpy], [],
        AC_CHECK_LIB([blas], [DAXPY_], [],
           AC_CHECK_LIB([blas], [DAXPY], [],
	      [AC_MSG_RESULT(not found, trying to use -lblas anyway.)]))))
  LDFLAGS="${LDFLAGS} -lblas"
  PKGCONFIG_LIBS_PRIVATE="${PKGCONFIG_LIBS_PRIVATE} -lblas"
else
  AC_DEFINE([INTERNAL_BLAS], [1], [Define to 1 if you use the internal BLAS library])
fi

if test "$internal_lapack" = "no"; then
  AC_CHECK_LIB([lapack], [dlarnv_], [],
     AC_CHECK_LIB([lapack], [dlarnv], [],
        AC_CHECK_LIB([lapack], [DLARNV_], [],
           AC_CHECK_LIB([lapack], [DLARNV], [],
	      [AC_MSG_RESULT(not found, trying to use -llapack anyway.)]))))
  LDFLAGS="${LDFLAGS} -llapack"
  PKGCONFIG_LIBS_PRIVATE="${PKGCONFIG_LIBS_PRIVATE} -llapack"
else
  AC_DEFINE([INTERNAL_LAPACK], [1], [Define to 1 if you use the internal LAPACK library])
fi

if test "$internal_arpack" = "no"; then
  if test "$tls_support" = "yes"; then
    AC_MSG_ERROR([Thread-local storage only supported with internal ARPACK library])
  fi
  AC_CHECK_LIB([arpack], [dsaupd_], [],
     AC_CHECK_LIB([arpack], [dsaupd], [],
        AC_CHECK_LIB([arpack], [DSAUPD_], [],
           AC_CHECK_LIB([arpack], [DSAUPD], [],
	      [AC_MSG_RESULT(not found, trying to use -larpack anyway.)]))))
  LDFLAGS="${LDFLAGS} -larpack"
  PKGCONFIG_LIBS_PRIVATE="${PKGCONFIG_LIBS_PRIVATE} -larpack"
else
  AC_DEFINE([INTERNAL_ARPACK], [1], [Define to 1 if you use the internal ARPACK library])
fi

glpk_support=no
AC_ARG_ENABLE(glpk, AC_HELP_STRING([--disable-glpk], [Compile without the GLPK library]))
if test "x$enable_glpk" != "xno"; then
  if test "$internal_glpk" = "no"; then
    AC_CHECK_LIB([glpk], [glp_read_mps], [
      AC_CHECK_HEADER([glpk.h], [
        AC_EGREP_CPP(yes, [
          #include <glpk.h>
          #if GLP_MAJOR_VERSION > 4 || (GLP_MAJOR_VERSION == 4 && GLP_MINOR_VERSION >= 38)
            yes
          #endif
        ], [
          AC_DEFINE([HAVE_GLPK], [1], [Define to 1 if you have the GLPK library])
          glpk_support=yes
          LDFLAGS="${LDFLAGS} -lglpk"
          PKGCONFIG_LIBS_PRIVATE="${PKGCONFIG_LIBS_PRIVATE} -lglpk"
        ])
      ])
    ])
  else
    AC_DEFINE([HAVE_GLPK], [1], [Define to 1 if you have the GLPK library])
    AC_DEFINE([INTERNAL_GLPK], [1], [Define to 1 if you use the internal GLPK library])
    glpk_support=yes
  fi
else
  internal_glpk=no
fi

AM_CONDITIONAL(INTERNAL_GLPK, test x$internal_glpk = xyes)
AM_CONDITIONAL(INTERNAL_ARPACK, test x$internal_arpack = xyes)
AM_CONDITIONAL(INTERNAL_LAPACK, test x$internal_lapack = xyes)
AM_CONDITIONAL(INTERNAL_BLAS, test x$internal_blas = xyes)
AM_CONDITIONAL(INTERNAL_F2C, test x$internal_f2c = xyes)

if test "$debug" = "yes"; then
  CFLAGS="${CFLAGS} -ggdb -O0"
  CPPFLAGS="${CPPFLAGS} -DRC_DEBUG"
  CXXFLAGS="${CXXFLAGS} -ggdb -O0"
fi

if test "$use_gprof" = "yes"; then
  CFLAGS="${CFLAGS} -pg"
  CXXFLAGS="${CXXFLAGS} -pg"
fi

if test "$use_asan" = "yes"; then
  CFLAGS="${CFLAGS} -g -fsanitize=address -fno-omit-frame-pointer"
  CXXFLAGS="${CXXFLAGS} -g -fsanitize=address -fno-omit-frame-pointer"
fi

if test "$use_asan" != "yes" -a "$use_gprof" != "yes" -a "$debug" != "yes"; then
  CFLAGS="${CFLAGS} -O3"
  CPPFLAGS="${CPPFLAGS} -O3"
  CXXFLAGS="${CXXFLAGS} -O3"
fi

AC_CONFIG_FILES([Makefile src/Makefile igraph.pc igraph_Info.plist doc/Makefile include/igraph_version.h include/igraph_threading.h])
AC_OUTPUT

AC_MSG_RESULT([igraph successfully configured.])
AC_MSG_RESULT([  GraphML format support -- $graphml_support])
AC_MSG_RESULT([  GMP library support    -- $gmp_support])
AC_MSG_RESULT([  GLPK library support   -- $glpk_support])
AC_MSG_RESULT([  Thread-local storage   -- $tls_support])
AC_MSG_RESULT([  Use internal ARPACK    -- $internal_arpack])
AC_MSG_RESULT([  Use internal LAPACK    -- $internal_lapack])
AC_MSG_RESULT([  Use internal BLAS      -- $internal_blas])
if test "$needs_f2c" != "yes"; then
  AC_MSG_RESULT([  Use internal F2C       -- f2c not needed])
else
  AC_MSG_RESULT([  Use internal F2C       -- $internal_f2c])
fi
if test "$glpk_support" != "no"; then
  AC_MSG_RESULT([  Use internal GLPK      -- $internal_glpk])
fi
AC_MSG_RESULT([  Debug build            -- $debug])
AC_MSG_RESULT([  Clang AddressSanitizer -- $use_asan])
AC_MSG_RESULT([  Profiling              -- $use_gprof])

