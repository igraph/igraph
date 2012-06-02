#!/bin/sh
#
# Creates documentation for igraph's Python interface using epydoc
#
# Usage: ./mkdoc.sh [--sync] [directory]

SCRIPTS_FOLDER=`dirname $0`

cd ${SCRIPTS_FOLDER}/..
ROOT_FOLDER=`pwd`
DOC_API_FOLDER=${ROOT_FOLDER}/doc/api

cd ${ROOT_FOLDER}
mkdir -p ${DOC_API_FOLDER}/pdf
mkdir -p ${DOC_API_FOLDER}/html

EPYDOC="${ROOT_FOLDER}/scripts/epydoc-patched"
python -m epydoc.__init__
if [ $? -gt 0 ]; then
  echo "Epydoc not installed, exiting..."
  exit 1
fi

PACKAGES="igraph igraph.statistics igraph.app igraph.app.shell"

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
${EPYDOC} --html -o ${DOC_API_FOLDER}/html -v \
          --name="IGraph library" \
	      --url="http://igraph.sourceforge.net" \
	      --no-private \
	      --exclude=igraph.test \
	      $PACKAGES

PDF=0
which latex >/dev/null && PDF=1

if [ $PDF -eq 1 ]; then
  echo "Generating PDF documentation..."
  ${EPYDOC} --pdf -o ${DOC_API_FOLDER}/pdf --exclude=igraph.test --inheritance=listed -v --name="IGraph library" --url="http://igraph.sourceforge.net" $PACKAGES

fi

if [ $SYNC -eq 1 ]; then
  echo "Syncing documentation to web"
  cp ${DOC_API_FOLDER}/pdf/igraph.pdf ${DOC_API_FOLDER}/html
  rsync --delete -avz ${DOC_API_FOLDER}/html/ ntamas,igraph@web.sourceforge.net:/home/groups/i/ig/igraph/htdocs/doc/python/
  rm ${DOC_API_FOLDER}/html/igraph.pdf
fi

cd "$PWD"
