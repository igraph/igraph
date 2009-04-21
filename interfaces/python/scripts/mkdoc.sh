#!/bin/sh
#
# Creates documentation for igraph's Python interface using epydoc
#
# Usage: ./mkdoc.sh [--sync] [directory]

PACKAGES="igraph igraph.clustering igraph.colors igraph.configuration igraph.datatypes igraph.drawing igraph.layout igraph.statistics igraph.app igraph.app.shell"

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
epydoc --html -o html -v \
       --name="IGraph library" \
	   --url="http://igraph.sourceforge.net" \
	   --no-private \
	   $PACKAGES

PDF=0
which latex >/dev/null && PDF=1

if [ $PDF -eq 1 ]; then
  echo "Generating PDF documentation..."
  epydoc --pdf -o pdf --inheritance=listed -v --name="IGraph library" --url="http://cneurocvs.rmki.kfki.hu/igraph" $PACKAGES

  echo "Zipping HTML documentation..."
  cd html && zip ../python-igraph-docs.zip -r . && cd ..

  echo "Moving PDF documentation to the HTML subdirectory..."
  mv pdf/api.pdf html/igraph.pdf
  rm -rf pdf
fi

if [ $SYNC -eq 1 ]; then
  echo "Syncing documentation to web"
  rsync --delete -avz html/ ntamas@cneurocvs.rmki.kfki.hu:/var/www/igraph-new/doc/python/
fi

cd "$PWD"
