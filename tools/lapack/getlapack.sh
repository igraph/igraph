#! /bin/sh
#
# ./getlapack.sh
#

# Troubleshooting note: Remember that this script does not re-download sources if
# it determines that they are already present. When you require a clean slate,
# delete them manually from the temporary directory.

BLAS_VERSION=3.12.0
LAPACK_VERSION=3.5.0
ARPACK_VERSION=3.7.0

# We can't go any further than LAPACK 3.5.0 because LAPACK 3.6.0 starts using
# recursive functions, which is a Fortran 90 construct and f2c can translate
# Fortran 77 only.

# We use ARPACK-NG instead of the original ARPACK. The latest feasible version
# is 3.7.0, as 3.8.0 starts using Fortran 90 features.

make

## Set the location of f2c headers below

f2cinclude=/opt/local/include

tempdir=/tmp
origdir=`pwd`
destdir=lapack-new

cd ${tempdir}
rm -rf $destdir
mkdir $destdir

## Download and unpack BLAS

if test ! -f blas.tgz; then
    echo "Downloading BLAS ${BLAS_VERSION} ..."
    curl -L -o blas.tgz http://www.netlib.org/blas/blas-${BLAS_VERSION}.tgz
fi
# blasdir=`tar tzf blas.tgz | head -1 | cut -f1 -d"/"`
## The following workaround is required for the source archives of some BLAS versions.
blasdir=BLAS-${BLAS_VERSION}
rm -rf ${blasdir}
tar xzf blas.tgz
echo "BLAS directory is ${blasdir}"

## Download, unpack and patch LAPACK

if test ! -f lapack.tgz; then
    echo "Downloading LAPACK ${LAPACK_VERSION} ..."
    curl -L -o lapack.tgz http://www.netlib.org/lapack/lapack-${LAPACK_VERSION}.tgz
fi
lapackdir=`tar tzf lapack.tgz | head -1 | cut -f1 -d"/"`
rm -rf ${lapackdir}
tar xzf lapack.tgz
echo "LAPACK directory is ${lapackdir}"

cd ${tempdir}/${lapackdir}
patch -p 1 <${origdir}/lapack.patch
cd ${tempdir}

## Download and unpack ARPACK

if test ! -f arpack.tar.gz; then
    echo "Downloading ARPACK ${ARPACK_VERSION} ..."
    curl -L -o arpack.tar.gz https://github.com/opencollab/arpack-ng/archive/refs/tags/${ARPACK_VERSION}.tar.gz
fi
arpackdir=`tar tzf arpack.tar.gz | head -1 | cut -f1 -d"/"`
rm -rf ${arpackdir}
tar xzf arpack.tar.gz
echo "ARPACK directory is ${arpackdir}"

## ARPACK-NG renames the original ARPACK function SECOND to ARSCND
## but keeps the file name. This is a timing function with several
## implementations. We use the one with an empty definition, since
## we do not need timing.

mv ${arpackdir}/UTIL/second_NONE.f ${arpackdir}/UTIL/arscnd.f

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
    cc -Wno-logical-op-parentheses -Wno-shift-op-parentheses \
       -I${f2cinclude} \
       -c ${name}.c >/dev/null &&
    nm ${name}.o | grep " U " | awk ' { print $2 }' |
    sed 's/_$//g' | sed 's/^_//g'
}

dofunction() {
    name=$1;

    if known $name; then return 0; fi

    echo "Processing function ${name} ... "

    ## Note: The order in which we look for functions in different packages matters
    ## because some functions are distributed with more than one package.

    if test -f ${tempdir}/${arpackdir}/SRC/${name}.f; then
        echo "Found ${name} in ARPACK."
        cd ${tempdir}/${arpackdir}/SRC
        arpack[$[${#arpack[@]}+1]]=$name
    elif test -f ${tempdir}/${lapackdir}/SRC/${name}.f; then
        echo "Found ${name} in LAPACK."
        cd ${tempdir}/${lapackdir}/SRC
        lapack[$[${#lapack[@]}+1]]=$name
    elif test -f ${tempdir}/${blasdir}/${name}.f; then
        echo "Found ${name} in BLAS."
        cd ${tempdir}/${blasdir}
        blas[$[${#blas[@]}+1]]=$name
    elif test -f ${tempdir}/${arpackdir}/UTIL/${name}.f; then
        echo "Found ${name} in ARPACK utils."
        cd ${tempdir}/${arpackdir}/UTIL
        arpack[$[${#arpack[@]}+1]]=$name
    elif test -f ${tempdir}/${lapackdir}/INSTALL/${name}.f; then
        echo "Found ${name} in LAPACK/INSTALL."
        cd ${tempdir}/${lapackdir}/INSTALL
        lapack[$[${#lapack[@]}+1]]=$name
    elif test -f ${origdir}/extra/${name}.f; then
        echo "Found ${name} in igraph extra."
        cd ${origdir}/extra
        lapack[$[${#lapack[@]}+1]]=$name
    else
        echo "Will skip ${name}."
        return
    fi

    cp ${name}.f ${tempdir}/${destdir}

    alreadydone[$[${#alreadydone[@]}+1]]=$name

    deps=`getdeps $name`
    for i in $deps; do
        dofunction $i
    done

    echo "Done with ${name}."
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

dofunction arscnd
dofunction dvout
dofunction ivout
dofunction dmout
dofunction dlamch
dofunction len_trim

## Polish them

cd ${tempdir}/${destdir}

# debug.h and stat.h contained common data blocks that we want to get rid of
# because it violates encapsulation. Therefore, we replace them with empty
# files, and patch the f2c-translated files later on to initialize the variables
# in these data blocks to zero.
touch debug.h
touch stat.h
trans_dir=${origdir} ${origdir}/CompletePolish *.f

## Remove the .f files.

cd ${tempdir}/${destdir}
rm -f *.f

## Prefix the function calls with 'igraph', this is needed
## if the user wants to link igraph including internal BLAS/LAPACK/ARPACK
## and BLAS/LAPACK/ARPACK for some reason

extrafunctions=(dlamc1 dlamc2 dlamc3 dlamc4 dlamc5 dnrm2)

for name in ${alreadydone[@]} ${extrafunctions[@]}; do
    echo "s/${name}_/igraph${name}_/g"
done > ${tempdir}/lapack-sed.txt

for name in ${alreadydone[@]}; do
    sed -f ${tempdir}/lapack-sed.txt < ${name}.c >${tempdir}/arpackfun.c
    mv ${tempdir}/arpackfun.c ${name}.c
done

## Update the file that is included into the main Makefile,
## this contains the ARPACK/LAPACK/BLAS source files

blasinc=${tempdir}/${destdir}/blas.inc
/bin/echo -n "BLAS = " > ${blasinc}
for name in ${blas[@]}; do
    /bin/echo -n "lapack/${name}.c "
done >> ${blasinc}
/bin/echo >> ${blasinc}

lapackinc=${tempdir}/${destdir}/lapack.inc
/bin/echo -n "LAPACK = " > ${lapackinc}
for name in ${lapack[@]}; do
    /bin/echo -n "lapack/${name}.c "
done | sed 's/lapack\/dlamch\.c//' >> ${lapackinc}
/bin/echo >> ${lapackinc}

arpackinc=${tempdir}/${destdir}/arpack.inc
/bin/echo -n "ARPACK = " > ${arpackinc}
for name in ${arpack[@]}; do
    /bin/echo -n "lapack/${name}.c "
done >> ${arpackinc}
/bin/echo >> ${arpackinc}

## This is a patch to make BLAS / LAPACK / ARPACK thread-safe

cd ${tempdir}/${destdir}
patch -l < ${origdir}/mt.patch

## We are done

echo "Sources are ready, to update your tree please run:

  git rm -rf ${origdir}/../../src/lapack
  mv ${tempdir}/${destdir} ${origdir}/../../src/lapack
  git add ${origdir}/../../src/lapack

"
