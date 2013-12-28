
all: nothing

VERSION=$(shell ../tools/getversion.sh)

HTML= index.html r/news.html c/news.html python/news.html \
      r/index.html c/index.html python/index.html         \
      _layouts/default.html _layouts/newstemp.html _layouts/manual.html

CSS= css/affix.css css/manual.css css/other.css fonts/fonts.css

RMAN:= $(wildcard ../interfaces/R/igraph/man/*)

TMP:=$(shell mktemp -d /tmp/.XXXXX)

../doc/jekyll/stamp: ../doc/html/stamp
	cd ../doc && make jekyll

c/doc/stamp: ../doc/jekyll/stamp
	rm -rf doc/c
	mkdir -p doc
	cp -r ../doc/jekyll doc/c

../doc/igraph-docs.pdf: ../doc/igraph-docs.xml
	cd ../doc && make igraph-docs.pdf

c/doc/igraph-docs.pdf: ../doc/igraph-docs.pdf
	mkdir -p doc/c
	cp ../doc/igraph-docs.pdf c/doc/

../doc/igraph.info: ../doc/igraph-docs.xml
	cd ../doc && make igraph.info

c/doc/igraph.info: ../doc/igraph.info
	mkdir -p doc/c
	cp ../doc/igraph.info c/doc/

r/doc/stamp: $(RMAN)
	cd ../interfaces/R && make && \
	R CMD INSTALL --html --no-R --no-configure --no-inst \
	  --no-libs --no-exec --no-test-load -l $(TMP) igraph
	rm -rf doc/r
	mkdir -p doc/r
	../tools/rhtml.sh $(TMP)/igraph/html doc/r
	ln -s 00Index.html r/doc/index.html
	touch r/doc/stamp

r/doc/igraph.pdf: $(RMAN)
	mkdir -p doc/r
	cd ../interfaces/R/ && make
	R CMD Rd2pdf --no-preview --force -o r/doc/igraph.pdf \
	  ../interfaces/R/igraph

../interfaces/python/doc/api/pdf/api.pdf:
	cd ../interfaces/python && python setup.py build \
		--c-core-url=http://igraph.org/nightly/get/c/igraph-$(VERSION).tar.gz
	cd ../interfaces/python && scripts/mkdoc.sh

python/doc/python-igraph.pdf: ../interfaces/python/doc/api/pdf/api.pdf
	mkdir -p doc/python
	cp $< $@

python/doc/stamp: ../interfaces/python/doc/api/html/igraph-module.html
	mkdir -p doc/python
	cp -r ../interfaces/python/doc/api/html/ doc/python
	../tools/pyhtml.sh doc/python
	touch $@

python/doc/tutorial/stamp: ../interfaces/python/doc/source/tutorial.rst
	mkdir -p python/doc/tutorial
	cd ../interfaces/python/doc && sphinx-build source api/tutorial
	cp -r ../interfaces/python/doc/api/tutorial/ python/doc/tutorial/
	touch $@

stamp: $(HTML) $(CSS) c/doc/stamp r/doc/stamp r/doc/igraph.pdf \
               c/doc/igraph.info c/doc/igraph-docs.pdf \
	       python/doc/python-igraph.pdf python/doc/stamp \
	       python/doc/tutorial/stamp
	../tools/getversion.sh > _includes/igraph-version
	../interfaces/R/tools/convertversion.sh > _includes/igraph-rversion
	jekyll build
	touch stamp

deploy: stamp
	rsync -avz --delete _site/ igraph.org:www/

.PHONY: all deploy

