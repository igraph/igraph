#! /bin/sh

## Quit immediately on error
set -e

## Want to run this as 'vagrant', so rerun if root
if [ "$(id -u)" = "0" ]; then
    sudo -u vagrant bash $0 $@
    exit 0
fi

## Check arguments
if [ ! "x"$# = "x1" ]; then exit 2; fi
version=$@
majorversion=`echo $version | cut -f1 -d.`
markerfile=~vagrant/vagrant/provisions/installR-${version}

echo -n "Installing R version $version.... "

## Check if we have anything to do
if [ -e ${markerfile} ]; then echo "Already installed" ; exit 0; fi

## Target directory
rdir=~vagrant/R/R-$version
mkdir -p ${rdir}

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
