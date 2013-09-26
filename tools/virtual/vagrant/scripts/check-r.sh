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
make parsersources
cd interfaces/R
make

## A temporary directory for R packages
libdir=`mktemp -d`
trap "rm -rf $libdir" EXIT

## Install dependent packages
${R} -e "
  options(repos=structure(c(CRAN='http://cran.rstudio.com/'))); \
  desc <- read.dcf('igraph/DESCRIPTION');                       \
  depkeys <- c('Depends', 'Imports', 'Suggests', 'LinkingTo');  \
  cn <- intersect(colnames(desc), depkeys);                     \
  pkg <- gsub(' ', '', unlist(strsplit(desc[,cn], ',')));       \
  install.packages(pkg, lib='$libdir', dependencies=NA);        \
"

${R} -e "
  .libPaths('$libdir');                                             \
  source('http://bioconductor.org/biocLite.R');                     \
  biocLite('graph', suppressUpdates=TRUE, suppressAutoUpdate=TRUE); \
"

package=`cat igraph/DESCRIPTION | grep ^Package: | cut -f2 -d" "`
version=`cat igraph/DESCRIPTION  | grep ^Version | cut -f2 -d" "`
commit=`git rev-parse --short HEAD`

## Check R package
R_LIBS=${libdir} ${R} CMD check --as-cran ${package}_${version}.tar.gz || true

## Upload the output
eval `ssh-agent -s` 
trap "kill $SSH_AGENT_PID" EXIT
ssh-add
ssh -p 2222 csardi@igraph.org mkdir -p www/nightly/check/r/${branch}/${commit}
scp -P 2222 ${package}.Rcheck/00check.log ${package}.Rcheck/00install.out \
    csardi@igraph.org:www/nightly/check/r/${branch}/${commit}/

## Clean up
rm -rf $builddir
rm -rf $libdir
kill $SSH_AGENT_PID

