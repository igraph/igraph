## Process this file with automake to produce Makefile.in

DOCINCLUDES = vector.xml error.xml matrix.xml basicigraph.xml generators.xml \
	      structural.xml iterators.xml attributes.xml layout.xml \
	      foreign.xml nongraph.xml

DOCFIX = fdl.xml gpl.xml igraph-docs.xml installation.xml introduction.xml \
	 tutorial.xml

all: doctemplate

html: html/stamp

info: igraph-docs.info

dvi: igraph-docs.dvi

pdf: igraph-docs.pdf

ps: igraph-docs.ps

doctemplate: $(DOCINCLUDES)

igraph-docs.xml: $(DOCINCLUDES)
	touch igraph-docs.xml

vector.xml: vector.xxml ../src/types.h ../src/vector.c
	./doxrox.py -t $< -e c-docbook.re -o $@ ../src/types.h ../src/vector.c

error.xml: error.xxml ../src/error.h
	./doxrox.py -t $< -e c-docbook.re -o $@ ../src/error.h

matrix.xml: matrix.xxml ../src/types.h ../src/matrix.c
	./doxrox.py -t $< -e c-docbook.re -o $@ ../src/types.h ../src/matrix.c

basicigraph.xml: basicigraph.xxml ../src/type_indexededgelist.c
	./doxrox.py -t $< -e c-docbook.re -o $@ ../src/type_indexededgelist.c

generators.xml: generators.xxml ../src/structure_generators.c ../src/games.c \
		../src/atlas.c
	./doxrox.py -t $< -e c-docbook.re -o $@ \
	../src/structure_generators.c ../src/games.c ../src/atlas.c

structural.xml: structural.xxml ../src/structural_properties.c \
		../src/conversion.c ../src/basic_query.c \
		../src/cocitation.c ../src/components.c
	./doxrox.py -t $< -e c-docbook.re -o $@ \
	../src/structural_properties.c ../src/conversion.c \
	../src/basic_query.c ../src/cocitation.c ../src/components.c

iterators.xml: iterators.xxml ../src/iterators.c
	./doxrox.py -t $< -e c-docbook.re -o $@ ../src/iterators.c

attributes.xml: attributes.xxml ../src/attributes.c
	./doxrox.py -t $< -e c-docbook.re -o $@ ../src/attributes.c

layout.xml: layout.xxml ../src/layout.c
	./doxrox.py -t $< -e c-docbook.re -o $@ ../src/layout.c

foreign.xml: foreign.xxml ../src/foreign.c
	./doxrox.py -t $< -e c-docbook.re -o $@ ../src/foreign.c

nongraph.xml: nongraph.xxml ../src/other.c ../src/random.c
	./doxrox.py -t $< -e c-docbook.re -o $@ ../src/other.c ../src/random.c


html/stamp: igraph-docs.xml $(DOCFIX)
	xmlto -x gtk-doc.xsl -o html xhtml igraph-docs.xml \
	&& touch html/stamp

igraph-docs.info: igraph-docs.xml $(DOCFIX)
	db2x_xsltproc --xinclude -s texi igraph-docs.xml -o igraph-docs.txml \
	&& db2x_texixml igraph-docs.txml \
	&& makeinfo igraph-docs.texi

igraph-docs.dvi: igraph-docs.xml $(DOCFIX)
	xmllint --xinclude -o igraph-docs2.xml igraph-docs.xml \
	&& docbook2dvi igraph-docs2.xml \
	&& mv igraph-docs2.dvi igraph-docs.dvi

igraph-docs.ps: igraph-docs.xml $(DOCFIX)
	xmllint --xinclude -o igraph-docs2.xml igraph-docs.xml \
	&& docbook2ps igraph-docs2.xml \
	&& mv igraph-docs2.ps igraph-docs.ps

igraph-docs.pdf: igraph-docs.xml $(DOCFIX)
	xmllint --xinclude -o igraph-docs2.xml igraph-docs.xml \
	&& docbook2pdf igraph-docs2.xml \
	&& mv igraph-docs2.pdf igraph-docs.pdf

CLEANFILES=$(DOCINCLUDES) html/*.html html/stamp \
	igraph-docs.{dvi,info,pdf,ps,texi,txml} igraph-docs2.xml 
