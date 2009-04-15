## Process this file with automake to produce Makefile.in

DOXROX=$(top_srcdir)/doc/doxrox.py
REGEX=$(top_srcdir)/doc/c-docbook.re
SRCDIR=$(top_srcdir)/src

DOCINCLUDES = vector.xml error.xml matrix.xml basicigraph.xml generators.xml \
	      structural.xml iterators.xml attributes.xml layout.xml \
	      foreign.xml nongraph.xml isomorphism.xml motifs.xml \
	      operators.xml flows.xml community.xml cliques.xml \
	      sparsematrix.xml stack.xml dqueue.xml heap.xml strvector.xml \
	      adjlist.xml arpack.xml bipartite.xml visitors.xml

DOCFIX = fdl.xml gpl.xml igraph-docs.xml installation.xml introduction.xml \
	 tutorial.xml licenses.xml pmt.xml

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

vector.xml: vector.xxml $(SRCDIR)/types.h $(SRCDIR)/vector.pmt $(SRCDIR)/vector_ptr.c $(SRCDIR)/vector.h $(SRCDIR)/vector.c
	$(DOXROX) -t $< -e $(REGEX) -o $@ $(SRCDIR)/types.h \
	$(SRCDIR)/vector.pmt $(SRCDIR)/vector_ptr.c $(SRCDIR)/vector.h $(SRCDIR)/vector.c

error.xml: error.xxml $(SRCDIR)/error.h
	$(DOXROX) -t $< -e $(REGEX) -o $@ $(SRCDIR)/error.h

matrix.xml: matrix.xxml $(SRCDIR)/types.h $(SRCDIR)/matrix.pmt
	$(DOXROX) -t $< -e $(REGEX) -o $@ $(SRCDIR)/types.h $(SRCDIR)/matrix.pmt

sparsematrix.xml: sparsematrix.xxml $(SRCDIR)/types.h $(SRCDIR)/spmatrix.c
	$(DOXROX) -t $< -e $(REGEX) -o $@ $(SRCDIR)/types.h $(SRCDIR)/spmatrix.c

basicigraph.xml: basicigraph.xxml $(SRCDIR)/type_indexededgelist.c $(SRCDIR)/structural_properties.c
	$(DOXROX) -t $< -e $(REGEX) -o $@ $(SRCDIR)/type_indexededgelist.c $(SRCDIR)/structural_properties.c

generators.xml: generators.xxml $(SRCDIR)/structure_generators.c $(SRCDIR)/games.c \
		$(SRCDIR)/atlas.c $(SRCDIR)/structural_properties.c $(SRCDIR)/forestfire.c
	$(DOXROX) -t $< -e $(REGEX) -o $@ \
	$(SRCDIR)/structure_generators.c $(SRCDIR)/games.c $(SRCDIR)/atlas.c \
	$(SRCDIR)/structural_properties.c $(SRCDIR)/forestfire.c

structural.xml: structural.xxml $(SRCDIR)/structural_properties.c \
		$(SRCDIR)/conversion.c $(SRCDIR)/basic_query.c \
		$(SRCDIR)/cocitation.c $(SRCDIR)/components.c \
		$(SRCDIR)/spectral_properties.c $(SRCDIR)/cores.c \
		$(SRCDIR)/centrality.c $(SRCDIR)/decomposition.c
	$(DOXROX) -t $< -e $(REGEX) -o $@ \
	$(SRCDIR)/structural_properties.c $(SRCDIR)/conversion.c \
	$(SRCDIR)/basic_query.c $(SRCDIR)/cocitation.c $(SRCDIR)/components.c \
	$(SRCDIR)/spectral_properties.c $(SRCDIR)/cores.c \
	$(SRCDIR)/centrality.c $(SRCDIR)/decomposition.c

iterators.xml: iterators.xxml $(SRCDIR)/iterators.c $(SRCDIR)/igraph.h
	$(DOXROX) -c -t $< -e $(REGEX) -o $@ $(SRCDIR)/iterators.c \
	$(SRCDIR)/igraph.h

attributes.xml: attributes.xxml $(SRCDIR)/attributes.c $(SRCDIR)/attributes.h $(SRCDIR)/cattributes.c
	$(DOXROX) -t $< -e $(REGEX) -o $@ $(SRCDIR)/attributes.c \
	$(SRCDIR)/attributes.h $(SRCDIR)/cattributes.c

layout.xml: layout.xxml $(SRCDIR)/layout.c $(SRCDIR)/igraph.h $(SRCDIR)/drl_layout.cpp $(SRCDIR)/drl_layout_3d.cpp
	$(DOXROX) -t $< -e $(REGEX) -o $@ $(SRCDIR)/layout.c $(SRCDIR)/igraph.h $(SRCDIR)/drl_layout.cpp $(SRCDIR)/drl_layout_3d.cpp

foreign.xml: foreign.xxml $(SRCDIR)/foreign.c $(SRCDIR)/foreign-graphml.c
	$(DOXROX) -t $< -e $(REGEX) -o $@ $(SRCDIR)/foreign.c \
	$(SRCDIR)/foreign-graphml.c

nongraph.xml: nongraph.xxml $(SRCDIR)/other.c $(SRCDIR)/random.c
	$(DOXROX) -t $< -e $(REGEX) -o $@ $(SRCDIR)/other.c $(SRCDIR)/random.c

isomorphism.xml: isomorphism.xxml $(SRCDIR)/topology.c $(SRCDIR)/igraph.h $(SRCDIR)/bliss.cc
	$(DOXROX) -c -t $< -e $(REGEX) -o $@ $(SRCDIR)/topology.c $(SRCDIR)/igraph.h $(SRCDIR)/bliss.cc

motifs.xml: motifs.xxml $(SRCDIR)/motifs.c
	$(DOXROX) -t $< -e $(REGEX) -o $@ $(SRCDIR)/motifs.c

operators.xml: operators.xxml $(SRCDIR)/operators.c
	$(DOXROX) -t $< -e $(REGEX) -o $@ $(SRCDIR)/operators.c

flows.xml: flows.xxml $(SRCDIR)/flow.c
	$(DOXROX) -t $< -e $(REGEX) -o $@ $(SRCDIR)/flow.c

community.xml: community.xxml $(SRCDIR)/community.c $(SRCDIR)/clustertool.cpp $(SRCDIR)/walktrap.cpp $(SRCDIR)/fast_community.c
	$(DOXROX) -t $< -e $(REGEX) -o $@ $(SRCDIR)/community.c $(SRCDIR)/clustertool.cpp $(SRCDIR)/walktrap.cpp $(SRCDIR)/fast_community.c

cliques.xml: cliques.xxml $(SRCDIR)/cliques.c
	$(DOXROX) -t $< -e $(REGEX) -o $@ $(SRCDIR)/cliques.c

stack.xml: stack.xxml $(SRCDIR)/stack.pmt
	$(DOXROX) -t $< -e $(REGEX) -o $@ $(SRCDIR)/stack.pmt

dqueue.xml: dqueue.xxml $(SRCDIR)/dqueue.pmt
	$(DOXROX) -t $< -e $(REGEX) -o $@ $(SRCDIR)/dqueue.pmt

heap.xml: heap.xxml $(SRCDIR)/heap.pmt
	$(DOXROX) -t $< -e $(REGEX) -o $@ $(SRCDIR)/heap.pmt

strvector.xml: strvector.xxml $(SRCDIR)/igraph_strvector.c $(SRCDIR)/types.h
	$(DOXROX) -t $< -e $(REGEX) -o $@ $(SRCDIR)/igraph_strvector.c $(SRCDIR)/types.h

adjlist.xml: adjlist.xxml $(SRCDIR)/adjlist.c $(SRCDIR)/igraph.h
	$(DOXROX) -c -t $< -e $(REGEX) -o $@ $(SRCDIR)/adjlist.c $(SRCDIR)/igraph.h

arpack.xml: arpack.xxml $(SRCDIR)/arpack.h $(SRCDIR)/arpack.c
	$(DOXROX) -t $< -e $(REGEX) -o $@ $(SRCDIR)/arpack.h $(SRCDIR)/arpack.c

bipartite.xml: bipartite.xxml $(SRCDIR)/bipartite.c
	$(DOXROX) -t $< -e $(REGEX) -o $@ $(SRCDIR)/bipartite.c

visitors.xml: visitors.xxml $(SRCDIR)/visitors.c  $(SRCDIR)/igraph.h
	$(DOXROX) -t $< -e $(REGEX) -o $@ $(SRCDIR)/visitors.c  $(SRCDIR)/igraph.h

html/stamp: igraph-docs.xml $(DOCFIX2) gtk-doc.xsl
	if [ "x$(top_srcdir)" != "x$(top_builddir)" ]; then cp $(DOCFIX2) . ; fi && \
	xmlto -x $(top_srcdir)/doc/gtk-doc.xsl -o html xhtml igraph-docs.xml \
	&& touch html/stamp

igraph-docs.info: igraph-docs.xml $(DOCFIX2) igraph-docs.info.diff
	if [ "x$(top_srcdir)" != "x$(top_builddir)" ]; then cp $(DOCFIX2) . ; fi && \
	db2x_xsltproc --xinclude -s texi igraph-docs.xml -o igraph-docs.txml \
	&& db2x_texixml igraph-docs.txml \
	&& mv igraph_reference_manual.texi igraph-docs.texi \
	&& makeinfo --no-split --number-sections igraph-docs.texi \
	&& mv igraph_reference_manual.info igraph-docs.info \
	&& patch <igraph-docs.info.diff

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
