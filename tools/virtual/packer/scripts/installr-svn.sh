#! /bin/sh -ex

## Quit immediately on error
set -e

echo -n "Installing R-devel..."

mkdir -p ~vagrant/src/
cd ~vagrant/src/
svn checkout https://svn.r-project.org/R/trunk/ R-devel
cd R-devel
./tools/rsync-recommended
./configure --prefix=$HOME/R/R-devel
make
make install
