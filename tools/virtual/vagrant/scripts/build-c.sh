# /bin/sh

## Quit immediately on error
set -e

## Build the C library from a github branch.
## We assume that all the build tools, and the 
## igraph dependencies are already installed.

## If not specified, we build the master branch
branch=${1:-master}

## We freshly clone the repo from github and build igraph from scratch.
builddir=`mktemp -d`
trap "rm -rf $builddir" EXIT
cd $builddir

git clone -b $branch https://github.com/igraph/igraph.git
cd igraph
./bootstrap.sh
./configure
make
cd tests ; make testsuite ; cd ..
if [ "$branch" = "0.5-main" ]; then
## Need the info file for the 0.5 tree
    cd doc ; make info ; cd ..
fi
make dist

## Canonical filename
version=`grep " VERSION " config.h | cut -f3 -d" " | tr -d '"'`
commit=`git rev-parse --short HEAD`
filename=igraph-${version}-${branch}-$commit.tar.gz
mv igraph-$version.tar.gz $filename

## Upload file to igraph.org
eval `ssh-agent -s` 
trap "kill $SSH_AGENT_PID" EXIT
ssh-add
scp -P 2222 ${filename} csardi@igraph.org:www/nightly/files/c/

## Clean up
rm -rf $builddir
kill $SSH_AGENT_PID
