
all: nothing

VERSION=$(shell ../tools/getversion.sh)

HTML= index.html r/news.html c/news.html \
      r/index.html c/index.html \
      _layouts/default.html _layouts/newstemp.html _layouts/manual.html

CSS= css/affix.css css/manual.css css/other.css fonts/fonts.css

RMAN:= $(wildcard ../interfaces/R/igraph/man/*)

TMP:=$(shell mktemp -d /tmp/.XXXXX)

../doc/jekyll/stamp: ../doc/html/stamp
	cd ../doc && make jekyll

../doc/html/stamp: $(wildcard ../src/*.c) ../doc/Makefile
	cd ../doc && make html

../doc/Makefile: ../doc/Makefile.am
	cd .. && ./bootstrap.sh
	cd .. && ./configure

c/doc/stamp: ../doc/jekyll/stamp
	rm -rf c/doc
	mkdir -p c
	cp -r ../doc/jekyll c/doc

../doc/igraph-docs.pdf: ../doc/igraph-docs.xml
	cd ../doc && make igraph-docs.pdf

c/doc/igraph-docs.pdf: ../doc/igraph-docs.pdf c/doc/stamp
	cp ../doc/igraph-docs.pdf c/doc/

../doc/igraph.info: ../doc/igraph-docs.xml c/doc/stamp
	cd ../doc && make igraph.info

c/doc/igraph.info: ../doc/igraph.info c/doc/stamp
	cp ../doc/igraph.info c/doc/

r/doc/stamp: $(RMAN)
	cd ../interfaces/R && make && \
	R CMD INSTALL --html --no-R --no-configure --no-inst \
	  --no-libs --no-exec --no-test-load -l $(TMP) igraph
	rm -rf r/doc
	mkdir -p r/doc
	../tools/rhtml.sh $(TMP)/igraph/html r/doc
	ln -s 00Index.html r/doc/index.html
	touch r/doc/stamp

r/doc/igraph.pdf: $(RMAN)
	mkdir -p r/doc
	cd ../interfaces/R/ && make
	R CMD Rd2pdf --no-preview --force -o r/doc/igraph.pdf \
	  ../interfaces/R/igraph

stamp: $(HTML) $(CSS) c/doc/stamp r/doc/stamp r/doc/igraph.pdf \
               c/doc/igraph.info c/doc/igraph-docs.pdf
	../tools/getversion.sh > _includes/igraph-version
	../interfaces/R/tools/convertversion.sh > _includes/igraph-rversion
	jekyll build
	touch stamp

deploy: stamp
	rsync -avz --delete --exclude python _site/ igraph.org:www2/www-test/

realdeploy: stamp
	rsync -avz --delete --exclude python _site/ igraph.org:www/

.PHONY: all deploy realdeploy

