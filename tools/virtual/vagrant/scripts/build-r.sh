#! /bin/sh

## Quit immediately an error
set -e

## Build the R package from a github branch.
## We assume that all the build tools, and the 
## igraph dependencies are already installed, 
## but the R packages we depend on are not.

## If not specified, we build the master branch
branch=${1-master}

## If not specified, we use the system R version
R=${2-R}

## We freshly clone the repo from github and build igraph from scratch.
builddir=`mktemp -d`
trap "rm -rf $builddir" EXIT
cd $builddir

git clone -b $branch https://github.com/igraph/igraph.git
cd igraph
./bootstrap.sh
./configure
./bootstrap.sh
make

version=`grep " VERSION " config.h | cut -f3 -d" " | tr -d '"'`
commit=`git rev-parse --short HEAD`

cd interfaces/R
make

## Canonical filename
if [ $branch = "0.5-main" ]; then
    package=igraph0
else
    package=igraph
fi
filename=${package}_${version}-${branch}-$commit.tar.gz
mv ${package}_${version}.tar.gz $filename

## Upload file to igraph.org
eval `ssh-agent -s` 
trap "kill $SSH_AGENT_PID" EXIT
ssh-add
scp -P 2222 ${filename} csardi@igraph.org:www/nightly/files/r/

## Clean up
rm -rf $builddir
kill $SSH_AGENT_PID

