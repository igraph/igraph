#!/bin/sh
#
# Creates documentation for igraph's Python interface using epydoc
#
# Usage: ./mkdoc.sh [--sync] [directory]

PWD=`pwd`

SYNC=0
if [ x$1 = x--sync ]; then
  SYNC=1
  shift
fi
if [ x$1 != x ]; then
  cd $1 || exit 1
fi

echo "Removing existing documentation..."
rm -rf html

echo "Generating HTML documentation..."
epydoc --html -o html -v --name="IGraph library" --url="http://cneurocvs.rmki.kfki.hu/igraph" igraph

PDF=0
if `which latex`; then PDF=1; fi

if [ $PDF -eq 1 ]; then
  echo "Generating PDF documentation..."
  epydoc --pdf -o pdf --inheritance=listed -v --name="IGraph library" --url="http://cneurocvs.rmki.kfki.hu/igraph" igraph

  echo "Moving PDF documentation to the HTML subdirectory..."
  mv pdf/api.pdf html/igraph.pdf
  rm -rf pdf
fi

if [ $SYNC -eq 1 ]; then
  echo "Syncing documentation to web"
  rsync --delete -avz html/ root@cneurocvs.rmki.kfki.hu:/var/www/igraph/doc/python/
fi

cd "$PWD"
