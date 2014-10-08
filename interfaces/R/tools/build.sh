#! /bin/bash

tempdir=$(mktemp -dt igraph)
trap "rm -rf ${tempdir}" EXIT
cp -r igraph/* ${tempdir}/

(
    cd ${tempdir}
    rm -rf src
    Rscript -e 'library(devtools) ; document()'
)

# cp ${tempdir}/DESCRIPTION igraph/
# cp ${tempdir}/NAMESPACE igraph/
cp ${tempdir}/man/* igraph/man/

rm -rf ${tempdir}

R CMD build igraph
