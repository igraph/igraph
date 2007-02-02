#!/bin/sh
#
# Creates documentation for igraph's Python interface using epydoc
#
# Usage: ./mkdoc.sh [--sync]

PWD=`pwd`
if [ ! -f mkdoc.sh ]; then
  cd scripts
  if [ ! -f mkdoc.sh ]; then
    echo "Must be run from the package root or from the scripts/ subdirectory"
    exit 1
  fi
  cd ..
else
  cd ..
fi

echo "Removing existing documentation..."
rm -rf html

echo "Generating HTML documentation..."
epydoc --html -o html --inheritance=included -v --no-sourcecode --name="IGraph library" --url="http://cneurocvs.rmki.kfki.hu/igraph" igraph

echo "Generating PDF documentation..."
epydoc --pdf -o pdf --inheritance=listed -v --no-sourcecode --name="IGraph library" --url="http://cneurocvs.rmki.kfki.hu/igraph" igraph

echo "Moving PDF documentation to the HTML subdirectory..."
mv pdf/api.pdf html/igraph.pdf
rm -rf pdf

if [ "x$1" == "x--sync" ]; then
  echo "Syncing documentation to web"
  rsync --delete -avz html/ root@cneurocvs.rmki.kfki.hu:/var/www/igraph/doc/python/
fi

cd "$PWD"
