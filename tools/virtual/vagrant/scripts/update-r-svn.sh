#! /bin/sh

## Quit immediately an error
set -e

if [ ! -e ~vagrant/src/R-devel ]; then echo "No R SVN folder"; exit 0; fi

cd ~vagrant/src/R-devel
svn update
./tools/rsync-recommended
./configure --prefix=$HOME/R/R-devel
make
make install

