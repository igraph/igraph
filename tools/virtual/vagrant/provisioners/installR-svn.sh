#! /bin/sh

## Quit immediately on error
set -e

## Want to run this as 'vagrant', so rerun if root
if [ "$(id -u)" = "0" ]; then
    sudo -u vagrant bash $0 $@
    exit 0
fi

echo -n "Installing R-devel..."

if [ -e ~vagrant/src/R-devel ]; then echo "already installed" ; exit 0 ; fi

mkdir -p ~vagrant/src/
cd ~vagrant/src/
svn checkout https://svn.r-project.org/R/trunk/ R-devel
cd R-devel
./tools/rsync-recommended
./configure --prefix=$HOME/R/R-devel
make
make install


