#! /bin/sh -ex

# Exit if error
set -e

echo $RVERSION

## Check R version to build
if [ -z "$RVERSION" ]; then exit 2; fi
version=$RVERSION
majorversion=`echo $version | cut -f1 -d.`
markerfile=~vagrant/R/R-$version/DONE

## Target directory
rdir=~vagrant/R/R-$version
mkdir -p ${rdir}

echo -n "Installing R version $version.... "

## Check if we have anything to do
if [ -e ${markerfile} ]; then echo "Already installed" ; exit 0; fi

## Temporary build directory
tmp=`mktemp -d`
cd ${tmp}

## Download, extract and build
wget http://cran.rstudio.com/src/base/R-${majorversion}/R-${version}.tar.gz
tar xzf R-${version}.tar.gz
cd R-${version}
./configure --prefix=${rdir}
make
make install

## Clean up
cd
rm -rf ${tmp}

## Mark this as done
echo DONE.
touch ${markerfile}
