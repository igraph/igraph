#! /bin/sh
#
# ./getlapack.sh dgeev dsyevr dnaupd dneupd dsaupd dseupd dgemv dgeevx \
#                dgetrf dgetrs dgesv dlapy2 dpotrf dsyrk dtrsv
#

BLAS_VERSION=3.8.0
LAPACK_VERSION=3.5.0

# We can't go any further than LAPACK 3.5.0 because LAPACK 3.6.0 starts using
# recursive functions, which is a Fortran 90 construct and f2c can translate
# Fortran 77 only.

make

origdir=`pwd`
destdir=lapack-new

cd /tmp
rm -rf $destdir
mkdir $destdir

## Download and unpack BLAS

if test ! -f blas.tgz; then
    curl -o blas.tgz http://www.netlib.org/blas/blas-${BLAS_VERSION}.tgz
fi
blasdir=`tar tzf blas.tgz | head -1 | cut -f1 -d"/"`
rm -rf ${blasdir}
tar xzf blas.tgz

## Download, unpack and patch LAPACK

if test ! -f lapack.tgz; then
    curl -o lapack.tgz http://www.netlib.org/lapack/lapack-${LAPACK_VERSION}.tgz
fi
lapackdir=`tar tzf lapack.tgz | head -1 | cut -f1 -d"/"`
rm -rf ${lapackdir}
tar xzf lapack.tgz

cd /tmp/${lapackdir}
patch -p 1 <${origdir}/lapack.patch
cd /tmp

## Download and unpack ARPACK

if test ! -f arpack96.tar.gz; then
    curl -O https://www.caam.rice.edu/software/ARPACK/SRC/arpack96.tar.gz
fi
arpackdir=`tar tzf arpack96.tar.gz | head -1 | cut -f1 -d"/"`
rm -rf ${arpackdir}
tar xzf arpack96.tar.gz

alreadydone=()
lapack=()
arpack=()
blas=()

known() {
    needle=$1
    res=0
    for i in ${alreadydone[@]}; do
	    if [[ $i == ${needle} ]]; then
		return 0
	    fi
    done
    return 1
}

getdeps() {
    name=$1;
    f2c -a ${name}.f >/dev/null 2>/dev/null &&
    gcc -Wno-logical-op-parentheses -Wno-shift-op-parentheses \
		-I/Users/tamas/include \
		-c ${name}.c >/dev/null &&
    nm ${name}.o | grep " U " | awk ' { print $2 }' |
    sed 's/_$//g' | sed 's/^_//g'
}

dofunction() {
    name=$1;

    if known $name; then return 0; fi

    if test -f /tmp/${arpackdir}/SRC/${name}.f; then
	cd /tmp/${arpackdir}/SRC
	arpack[$[${#arpack[@]}+1]]=$name
    elif test -f /tmp/${lapackdir}/SRC/${name}.f; then
	cd /tmp/${lapackdir}/SRC
	lapack[$[${#lapack[@]}+1]]=$name
    elif test -f /tmp/${blasdir}/${name}.f; then
	cd /tmp/${blasdir}
	blas[$[${#blas[@]}+1]]=$name
    elif test -f /tmp/${arpackdir}/UTIL/${name}.f; then
	cd /tmp/${arpackdir}/UTIL
	arpack[$[${#arpack[@]}+1]]=$name
    elif test -f /tmp/${lapackdir}/INSTALL/${name}.f; then
	cd /tmp/${lapackdir}/INSTALL
	lapack[$[${#lapack[@]}+1]]=$name
    elif test -f ${origdir}/extra/${name}.f; then
	cd ${origdir}/extra
	lapack[$[${#lapack[@]}+1]]=$name
    else
	return
    fi

    cp ${name}.f /tmp/${destdir}

    alreadydone[$[${#alreadydone[@]}+1]]=$name

    deps=`getdeps $name`
    for i in $deps; do
	dofunction $i
    done
}

## Collect and copy the needed files

FUNCS="$@"
if [ "x$FUNCS" = x ]; then
	FUNCS="dgeev dsyevr dnaupd dneupd dsaupd dseupd dgemv dgeevx dgetrf dgetrs dgesv dlapy2 dpotrf dsyrk dtrsv"
fi

for i in $FUNCS; do
    dofunction $i
done

## Some more required files

dofunction second
dofunction dvout
dofunction ivout
dofunction dmout
dofunction dlamch
dofunction len_trim

## Polish them

cd /tmp/${destdir}

# debug.h and stat.h contained common data blocks that we want to get rid of
# because it violates encapsulation. Therefore, we replace them with empty
# files, and patch the f2c-translated files later on to initialize the variables
# in these data blocks to zero.
touch debug.h
touch stat.h
trans_dir=${origdir} ${origdir}/CompletePolish *.f

## Remove the .f files.

cd /tmp/${destdir}
rm -f *.f

## Prefix the function calls with 'igraph', this is needed
## if the user wants to link igraph including internal BLAS/LAPACK/ARPACK
## and BLAS/LAPACK/ARPACK for some reason

extrafunctions=(dlamc1 dlamc2 dlamc3 dlamc4 dlamc5)

for name in ${alreadydone[@]} ${extrafunctions[@]}; do
    echo "s/${name}_/igraph${name}_/g"
done > /tmp/lapack-sed.txt

for name in ${alreadydone[@]}; do
    sed -f /tmp/lapack-sed.txt < ${name}.c >/tmp/arpackfun.c
    mv /tmp/arpackfun.c ${name}.c
done

## Update the file that is included into the main Makefile,
## this contains the ARPACK/LAPACK/BLAS source files

blasinc=/tmp/${destdir}/blas.inc
/bin/echo -n "BLAS = " > ${blasinc}
for name in ${blas[@]}; do
    /bin/echo -n "lapack/${name}.c "
done >> ${blasinc}
/bin/echo >> ${blasinc}

lapackinc=/tmp/${destdir}/lapack.inc
/bin/echo -n "LAPACK = " > ${lapackinc}
for name in ${lapack[@]}; do
    /bin/echo -n "lapack/${name}.c "
done | sed 's/lapack\/dlamch\.c//' >> ${lapackinc}
/bin/echo >> ${lapackinc}

arpackinc=/tmp/${destdir}/arpack.inc
/bin/echo -n "ARPACK = " > ${arpackinc}
for name in ${arpack[@]}; do
    /bin/echo -n "lapack/${name}.c "
done >> ${arpackinc}
/bin/echo >> ${arpackinc}

## This is a patch to make BLAS / LAPACK / ARPACK thread-safe

cd /tmp/${destdir}
patch -p2 < ${origdir}/mt.patch

## We are done

echo "Sources are ready, to update your tree please run:

  git rm -rf ${origdir}/../../src/lapack
  mv /tmp/${destdir} ${origdir}/../../src/lapack
  git add ${origdir}/../../src/lapack

"
