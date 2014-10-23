#! /bin/bash

tempdir=$(mktemp -d 2>/dev/null || mktemp -d -t 'mytmpdir')
trap "rm -rf ${tempdir}" EXIT
cp -r igraph/* ${tempdir}/

(
    cd ${tempdir}
    rm -rf src
    Rscript -e 'library(devtools) ; document()'
)

cp ${tempdir}/DESCRIPTION igraph/
cp ${tempdir}/NAMESPACE igraph/
cp ${tempdir}/man/* igraph/man/

rm -rf ${tempdir}
