## Process this file with automake to produce Makefile.in

DOXROX=$(top_srcdir)/doc/doxrox.py
REGEX=$(top_srcdir)/doc/c-docbook.re
SRCDIR=$(top_srcdir)/src

DOCINCLUDES = vector.xml error.xml matrix.xml basicigraph.xml generators.xml \
	      structural.xml iterators.xml attributes.xml layout.xml \
	      foreign.xml nongraph.xml isomorphism.xml motifs.xml \
	      operators.xml flows.xml community.xml cliques.xml

DOCFIX = fdl.xml gpl.xml igraph-docs.xml installation.xml introduction.xml \
	 tutorial.xml licenses.xml

DOCFIX2 = $(patsubst %,$(top_srcdir)/doc/%,$(DOCFIX))

all: doctemplate

html: html/stamp

info: igraph-docs.info

dvi: igraph-docs.dvi

pdf: igraph-docs.pdf

ps: igraph-docs.ps

doctemplate: $(DOCINCLUDES)

igraph-docs.xml: $(DOCINCLUDES)
	touch igraph-docs.xml

vector.xml: vector.xxml $(SRCDIR)/types.h $(SRCDIR)/vector.c $(SRCDIR)/vector_ptr.c
	$(DOXROX) -t $< -e $(REGEX) -o $@ $(SRCDIR)/types.h \
	$(SRCDIR)/vector.c $(SRCDIR)/vector_ptr.c

error.xml: error.xxml $(SRCDIR)/error.h
	$(DOXROX) -t $< -e $(REGEX) -o $@ $(SRCDIR)/error.h

matrix.xml: matrix.xxml $(SRCDIR)/types.h $(SRCDIR)/matrix.c
	$(DOXROX) -t $< -e $(REGEX) -o $@ $(SRCDIR)/types.h $(SRCDIR)/matrix.c

basicigraph.xml: basicigraph.xxml $(SRCDIR)/type_indexededgelist.c
	$(DOXROX) -t $< -e $(REGEX) -o $@ $(SRCDIR)/type_indexededgelist.c

generators.xml: generators.xxml $(SRCDIR)/structure_generators.c $(SRCDIR)/games.c \
		$(SRCDIR)/atlas.c $(SRCDIR)/structural_properties.c
	$(DOXROX) -t $< -e $(REGEX) -o $@ \
	$(SRCDIR)/structure_generators.c $(SRCDIR)/games.c $(SRCDIR)/atlas.c \
	$(SRCDIR)/structural_properties.c

structural.xml: structural.xxml $(SRCDIR)/structural_properties.c \
		$(SRCDIR)/conversion.c $(SRCDIR)/basic_query.c \
		$(SRCDIR)/cocitation.c $(SRCDIR)/components.c \
		$(SRCDIR)/spectral_properties.c $(SRCDIR)/cores.c
	$(DOXROX) -t $< -e $(REGEX) -o $@ \
	$(SRCDIR)/structural_properties.c $(SRCDIR)/conversion.c \
	$(SRCDIR)/basic_query.c $(SRCDIR)/cocitation.c $(SRCDIR)/components.c \
	$(SRCDIR)/spectral_properties.c $(SRCDIR)/cores.c

iterators.xml: iterators.xxml $(SRCDIR)/iterators.c $(SRCDIR)/igraph.h
	$(DOXROX) -t $< -e $(REGEX) -o $@ $(SRCDIR)/iterators.c \
	$(SRCDIR)/igraph.h

attributes.xml: attributes.xxml $(SRCDIR)/attributes.c $(SRCDIR)/attributes.h
	$(DOXROX) -t $< -e $(REGEX) -o $@ $(SRCDIR)/attributes.c \
	$(SRCDIR)/attributes.h

layout.xml: layout.xxml $(SRCDIR)/layout.c
	$(DOXROX) -t $< -e $(REGEX) -o $@ $(SRCDIR)/layout.c

foreign.xml: foreign.xxml $(SRCDIR)/foreign.c $(SRCDIR)/foreign-graphml.c
	$(DOXROX) -t $< -e $(REGEX) -o $@ $(SRCDIR)/foreign.c \
	$(SRCDIR)/foreign-graphml.c

nongraph.xml: nongraph.xxml $(SRCDIR)/other.c $(SRCDIR)/random.c
	$(DOXROX) -t $< -e $(REGEX) -o $@ $(SRCDIR)/other.c $(SRCDIR)/random.c

isomorphism.xml: isomorphism.xxml $(SRCDIR)/topology.c 
	$(DOXROX) -t $< -e $(REGEX) -o $@ $(SRCDIR)/topology.c

motifs.xml: motifs.xxml $(SRCDIR)/motifs.c
	$(DOXROX) -t $< -e $(REGEX) -o $@ $(SRCDIR)/motifs.c

operators.xml: operators.xxml $(SRCDIR)/operators.c
	$(DOXROX) -t $< -e $(REGEX) -o $@ $(SRCDIR)/operators.c

flows.xml: flows.xxml $(SRCDIR)/flow.c
	$(DOXROX) -t $< -e $(REGEX) -o $@ $(SRCDIR)/flow.c

community.xml: community.xxml $(SRCDIR)/community.c $(SRCDIR)/clustertool.cpp $(SRCDIR)/walktrap.cpp $(SRCDIR)/fastcommunity_mh.cc
	$(DOXROX) -t $< -e $(REGEX) -o $@ $(SRCDIR)/community.c $(SRCDIR)/clustertool.cpp $(SRCDIR)/walktrap.cpp $(SRCDIR)/fastcommunity_mh.cc

cliques.xml: cliques.xxml $(SRCDIR)/cliques.c
	$(DOXROX) -t $< -e $(REGEX) -o $@ $(SRCDIR)/cliques.c

html/stamp: igraph-docs.xml $(DOCFIX2) gtk-doc.xsl
	if [ "x$(top_srcdir)" != "x$(top_builddir)" ]; then cp $(DOCFIX2) . ; fi && \
	xmlto -x $(top_srcdir)/doc/gtk-doc.xsl -o html xhtml igraph-docs.xml \
	&& touch html/stamp

igraph-docs.info: igraph-docs.xml $(DOCFIX2) igraph-docs.info.diff
	if [ "x$(top_srcdir)" != "x$(top_builddir)" ]; then cp $(DOCFIX2) . ; fi && \
	db2x_xsltproc --xinclude -s texi igraph-docs.xml -o igraph-docs.txml \
	&& db2x_texixml igraph-docs.txml \
	&& makeinfo --no-split igraph-docs.texi \
	&& patch -p2 <igraph-docs.info.diff

igraph-docs.dvi: igraph-docs.xml $(DOCFIX2)
	if [ "x$(top_srcdir)" != "x$(top_builddir)" ]; then cp $(DOCFIX2) . ; fi && \
	xmllint --xinclude -o igraph-docs2.xml igraph-docs.xml \
	&& docbook2dvi igraph-docs2.xml \
	&& mv igraph-docs2.dvi igraph-docs.dvi

igraph-docs.ps: igraph-docs.xml $(DOCFIX2)
	if [ "x$(top_srcdir)" != "x$(top_builddir)" ]; then cp $(DOCFIX2) . ; fi && \
	xmllint --xinclude -o igraph-docs2.xml igraph-docs.xml \
	&& docbook2ps igraph-docs2.xml \
	&& mv igraph-docs2.ps igraph-docs.ps

igraph-docs.pdf: igraph-docs.xml $(DOCFIX2)
	if [ "x$(top_srcdir)" != "x$(top_builddir)" ]; then cp $(DOCFIX2) . ; fi && \
	xmllint --xinclude -o igraph-docs2.xml igraph-docs.xml \
	&& docbook2pdf igraph-docs2.xml \
	&& mv igraph-docs2.pdf igraph-docs.pdf

CLEANFILES=$(DOCINCLUDES) html/*.html html/stamp \
	igraph-docs.{dvi,info,pdf,ps,texi,txml} igraph-docs2.xml 
