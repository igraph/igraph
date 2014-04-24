#! /bin/sh

## If not specified, we build the master branch
branch=${1-master}

## R version with AddressSanitizer support
R=~vagrant/R/R-3.0.1-asan/bin/R

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
version=`cat igraph/DESCRIPTION | grep ^Version: | cut -f2 -d" "`

echo "
CC = clang -fsanitize=address -fno-omit-frame-pointer
CXX = clang++ -fsanitize=address -fno-omit-frame-pointer
" >${builddir}/Makevars

R_MAKEVARS_USER=${builddir}/Makevars R_LIBS=${libdir} \
    ${R} CMD INSTALL -l ${libdir} ${package}_${version}.tar.gz

## Extract examples and run them
${R} -e "
  library(tools);                                                              \
  rdfiles <- list.files('igraph/man', pattern='.*\\\\.Rd$', full.names=TRUE);  \
  out <- file('igraph-Ex.R', open='w');                                        \
  cat('### Load the package\\n.libPaths(\'${libdir}\');library(graph);library(\'${package}\')\\n\\n', \
      file=out);                                                               \
  sapply(rdfiles, Rd2ex, out=out);                                             \
  close(out)                                                                   \
"
${R} --no-save < igraph-Ex.R

rm -rf $builddir
rm -rf $libdir
