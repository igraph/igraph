#! /bin/bash

cd igraph/tests
tmpdir=`mktemp -d $TMPDIR/igraph-XXXXXX`

## Copy unknown files to the test directory, just to be sure

otherfiles=`ls | grep -v '\.R$' | grep -v '\.Rout$' | grep -v '\.Rout.save'`
cp $otherfiles $tmpdir

if [ -z "$@" ]; then
    ## If no arguments, then run all tests
    cp * $tmpdir
else
    ## Else copy all test files that match any of the arguments
    for arg in "$@"; do
	cp `ls | grep $arg` $tmpdir
    done
fi

## Run the tests
cd $tmpdir && echo "tools:::.runPackageTestsR()" | R --no-save && echo

## Echo temporary directory name
echo "Tests were run in"
echo $tmpdir
