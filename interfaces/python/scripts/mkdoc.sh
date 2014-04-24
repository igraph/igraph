#!/bin/sh
#
# Creates documentation for igraph's Python interface using epydoc
#
# Usage: ./mkdoc.sh [--sync] [directory]

SCRIPTS_FOLDER=`dirname $0`

cd ${SCRIPTS_FOLDER}/..
ROOT_FOLDER=`pwd`
DOC_API_FOLDER=${ROOT_FOLDER}/doc/api
CONFIG=${ROOT_FOLDER}/scripts/epydoc.cfg

cd ${ROOT_FOLDER}
mkdir -p ${DOC_API_FOLDER}/pdf
mkdir -p ${DOC_API_FOLDER}/html

EPYDOC="${ROOT_FOLDER}/scripts/epydoc-patched"
python -m epydoc.__init__
if [ $? -gt 0 ]; then
  echo "Epydoc not installed, exiting..."
  exit 1
fi

PWD=`pwd`

SYNC=0
if [ x$1 = x--sync ]; then
  SYNC=1
  shift
fi
if [ x$1 != x ]; then
  cd $1 || exit 1
fi

echo "Checking symlinked _igraph.so in ${ROOT_FOLDER}/igraph..."
if [ ! -e ${ROOT_FOLDER}/igraph/_igraph.so -o ! -L ${ROOT_FOLDER}/igraph/_igraph.so ]; then
	rm -f ${ROOT_FOLDER}/igraph/_igraph.so
	cd ${ROOT_FOLDER}/igraph
	ln -s ../build/lib*/igraph/_igraph.so .
	cd ${ROOT_FOLDER}
fi

echo "Removing existing documentation..."
rm -rf html

echo "Generating HTML documentation..."
${EPYDOC} --html -v -o ${DOC_API_FOLDER}/html --config ${CONFIG}

PDF=0
which latex >/dev/null && PDF=1

if [ $PDF -eq 1 ]; then
  echo "Generating PDF documentation..."
${EPYDOC} --pdf -v -o ${DOC_API_FOLDER}/pdf --config ${CONFIG}
fi

if [ $SYNC -eq 1 ]; then
  echo "Syncing documentation to web"
  cp ${DOC_API_FOLDER}/pdf/igraph.pdf ${DOC_API_FOLDER}/html
  rsync --delete -avz ${DOC_API_FOLDER}/html/ csardi@igraph.org:2222:www/doc/python/
  rm ${DOC_API_FOLDER}/html/igraph.pdf
fi

cd "$PWD"
