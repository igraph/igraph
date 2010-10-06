## Process this file with automake to produce Makefile.in

DOXROX=$(top_srcdir)/doc/doxrox.py
REGEX=$(top_srcdir)/doc/c-docbook.re
SRCDIR=$(top_srcdir)/src
INCLUDEDIR=$(top_srcdir)/include

DOCINCLUDES = vector.xml error.xml matrix.xml basicigraph.xml generators.xml \
	      structural.xml iterators.xml attributes.xml layout.xml \
	      foreign.xml nongraph.xml isomorphism.xml motifs.xml \
	      operators.xml flows.xml community.xml cliques.xml \
	      sparsematrix.xml stack.xml dqueue.xml heap.xml strvector.xml \
	      adjlist.xml arpack.xml bipartite.xml visitors.xml random.xml \
	      separators.xml memory.xml

DOCFIX = fdl.xml gpl.xml igraph-docs.xml installation.xml introduction.xml \
	 tutorial.xml licenses.xml pmt.xml

DOCFIX2 = $(patsubst %,$(top_srcdir)/doc/%,$(DOCFIX))

all: doctemplate

html: html/stamp

info: igraph.info

dvi: igraph-docs.dvi

pdf: igraph-docs.pdf

ps: igraph-docs.ps

doctemplate: $(DOCINCLUDES)

igraph-docs.xml: $(DOCINCLUDES)
	touch igraph-docs.xml

vector.xml: vector.xxml $(INCLUDEDIR)/igraph_vector_pmt.h $(INCLUDEDIR)/igraph_vector_ptr.h $(SRCDIR)/vector.pmt $(SRCDIR)/vector_ptr.c $(SRCDIR)/vector.c
	$(DOXROX) -t $< -e $(REGEX) -o $@ $(INCLUDEDIR)/igraph_vector_pmt.h \
	$(INCLUDEDIR)/igraph_vector_ptr.h $(SRCDIR)/vector.pmt \
	$(SRCDIR)/vector_ptr.c $(SRCDIR)/vector.c

error.xml: error.xxml $(INCLUDEDIR)/igraph_error.h
	$(DOXROX) -t $< -e $(REGEX) -o $@ $(INCLUDEDIR)/igraph_error.h

matrix.xml: matrix.xxml $(INCLUDEDIR)/igraph_matrix.h $(SRCDIR)/matrix.pmt
	$(DOXROX) -t $< -e $(REGEX) -o $@ $(INCLUDEDIR)/igraph_matrix.h $(SRCDIR)/matrix.pmt

sparsematrix.xml: sparsematrix.xxml $(INCLUDEDIR)/igraph_spmatrix.h $(SRCDIR)/spmatrix.c
	$(DOXROX) -t $< -e $(REGEX) -o $@ $(INCLUDEDIR)/igraph_spmatrix.h $(SRCDIR)/spmatrix.c

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
		$(SRCDIR)/centrality.c $(SRCDIR)/decomposition.c \
		$(SRCDIR)/mixing.c
	$(DOXROX) -t $< -e $(REGEX) -o $@ \
	$(SRCDIR)/structural_properties.c $(SRCDIR)/conversion.c \
	$(SRCDIR)/basic_query.c $(SRCDIR)/cocitation.c $(SRCDIR)/components.c \
	$(SRCDIR)/spectral_properties.c $(SRCDIR)/cores.c \
	$(SRCDIR)/centrality.c $(SRCDIR)/decomposition.c $(SRCDIR)/mixing.c

iterators.xml: iterators.xxml $(SRCDIR)/iterators.c $(INCLUDEDIR)/igraph_iterators.h
	$(DOXROX) -c -t $< -e $(REGEX) -o $@ $(SRCDIR)/iterators.c \
	$(INCLUDEDIR)/igraph_iterators.h

attributes.xml: attributes.xxml $(SRCDIR)/attributes.c $(INCLUDEDIR)/igraph_attributes.h $(SRCDIR)/cattributes.c
	$(DOXROX) -t $< -e $(REGEX) -o $@ $(SRCDIR)/attributes.c \
	$(INCLUDEDIR)/igraph_attributes.h $(SRCDIR)/cattributes.c

layout.xml: layout.xxml $(SRCDIR)/layout.c $(INCLUDEDIR)/igraph_layout.h $(SRCDIR)/drl_layout.cpp $(SRCDIR)/drl_layout_3d.cpp
	$(DOXROX) -t $< -e $(REGEX) -o $@ $(SRCDIR)/layout.c $(INCLUDEDIR)/igraph_layout.h $(SRCDIR)/drl_layout.cpp $(SRCDIR)/drl_layout_3d.cpp

foreign.xml: foreign.xxml $(SRCDIR)/foreign.c $(SRCDIR)/foreign-graphml.c
	$(DOXROX) -t $< -e $(REGEX) -o $@ $(SRCDIR)/foreign.c \
	$(SRCDIR)/foreign-graphml.c

nongraph.xml: nongraph.xxml $(SRCDIR)/other.c $(SRCDIR)/random.c $(SRCDIR)/version.c
	$(DOXROX) -t $< -e $(REGEX) -o $@ $(SRCDIR)/other.c $(SRCDIR)/random.c $(SRCDIR)/version.c

isomorphism.xml: isomorphism.xxml $(SRCDIR)/topology.c $(INCLUDEDIR)/igraph_topology.h $(SRCDIR)/bliss.cc
	$(DOXROX) -c -t $< -e $(REGEX) -o $@ $(SRCDIR)/topology.c $(INCLUDEDIR)/igraph_topology.h $(SRCDIR)/bliss.cc

motifs.xml: motifs.xxml $(SRCDIR)/motifs.c
	$(DOXROX) -t $< -e $(REGEX) -o $@ $(SRCDIR)/motifs.c

operators.xml: operators.xxml $(SRCDIR)/operators.c
	$(DOXROX) -t $< -e $(REGEX) -o $@ $(SRCDIR)/operators.c

flows.xml: flows.xxml $(SRCDIR)/flow.c $(SRCDIR)/st-cuts.c
	$(DOXROX) -t $< -e $(REGEX) -o $@ $(SRCDIR)/flow.c $(SRCDIR)/st-cuts.c

community.xml: community.xxml $(SRCDIR)/community.c $(SRCDIR)/clustertool.cpp $(SRCDIR)/walktrap.cpp $(SRCDIR)/fast_community.c $(SRCDIR)/optimal_modularity.c
	$(DOXROX) -t $< -e $(REGEX) -o $@ $(SRCDIR)/community.c $(SRCDIR)/clustertool.cpp $(SRCDIR)/walktrap.cpp $(SRCDIR)/fast_community.c $(SRCDIR)/optimal_modularity.c

cliques.xml: cliques.xxml $(SRCDIR)/cliques.c
	$(DOXROX) -t $< -e $(REGEX) -o $@ $(SRCDIR)/cliques.c

stack.xml: stack.xxml $(SRCDIR)/stack.pmt
	$(DOXROX) -t $< -e $(REGEX) -o $@ $(SRCDIR)/stack.pmt

dqueue.xml: dqueue.xxml $(SRCDIR)/dqueue.pmt
	$(DOXROX) -t $< -e $(REGEX) -o $@ $(SRCDIR)/dqueue.pmt

heap.xml: heap.xxml $(SRCDIR)/heap.pmt
	$(DOXROX) -t $< -e $(REGEX) -o $@ $(SRCDIR)/heap.pmt

strvector.xml: strvector.xxml $(SRCDIR)/igraph_strvector.c $(INCLUDEDIR)/igraph_strvector.h
	$(DOXROX) -t $< -e $(REGEX) -o $@ $(SRCDIR)/igraph_strvector.c $(INCLUDEDIR)/igraph_strvector.h

adjlist.xml: adjlist.xxml $(SRCDIR)/adjlist.c $(INCLUDEDIR)/igraph_adjlist.h
	$(DOXROX) -c -t $< -e $(REGEX) -o $@ $(SRCDIR)/adjlist.c $(INCLUDEDIR)/igraph_adjlist.h

arpack.xml: arpack.xxml $(INCLUDEDIR)/igraph_arpack.h $(SRCDIR)/arpack.c $(INCLUDEDIR)/igraph_blas.h $(SRCDIR)/blas.c
	$(DOXROX) -t $< -e $(REGEX) -o $@ $(INCLUDEDIR)/igraph_arpack.h $(SRCDIR)/arpack.c $(INCLUDEDIR)/igraph_blas.h $(SRCDIR)/blas.c

bipartite.xml: bipartite.xxml $(SRCDIR)/bipartite.c
	$(DOXROX) -t $< -e $(REGEX) -o $@ $(SRCDIR)/bipartite.c

visitors.xml: visitors.xxml $(SRCDIR)/visitors.c  $(INCLUDEDIR)/igraph_visitor.h
	$(DOXROX) -t $< -e $(REGEX) -o $@ $(SRCDIR)/visitors.c  $(INCLUDEDIR)/igraph_visitor.h

random.xml: random.xxml $(SRCDIR)/random.c
	$(DOXROX) -t $< -e $(REGEX) -o $@ $(SRCDIR)/random.c

separators.xml: separators.xxml $(SRCDIR)/separators.c $(SRCDIR)/st-cuts.c
	$(DOXROX) -t $< -e $(REGEX) -o $@ $(SRCDIR)/separators.c \
		$(SRCDIR)/st-cuts.c

memory.xml: memory.xxml $(SRCDIR)/memory.c
	$(DOXROX) -t $< -e $(REGEX) -o $@ $(SRCDIR)/memory.c

html/stamp: igraph-docs.xml $(DOCFIX2) gtk-doc.xsl
	if [ "x$(top_srcdir)" != "x$(top_builddir)" ]; then cp $(DOCFIX2) . ; fi && \
	xmlto -x $(top_srcdir)/doc/gtk-doc.xsl -o html xhtml igraph-docs.xml \
	&& touch html/stamp

igraph.info: igraph-docs.xml $(DOCFIX2) igraph.info.diff
	if [ "x$(top_srcdir)" != "x$(top_builddir)" ]; then cp $(DOCFIX2) . ; fi && \
	cat tutorial.xml | sed '/<link/{N;s/<link\n*[^>]*>//g}' | \
	  sed 's/<\/link>//g' >tutorial-info.xml \
	&& sed 's/tutorial\.xml/tutorial-info.xml/' igraph-docs.xml > igraph-docs-info.xml \
	&& db2x_xsltproc --xinclude -s texi igraph-docs-info.xml -o igraph.txml \
	&& db2x_texixml igraph.txml \
	&& mv igraph_reference_manual.texi igraph.texi \
	&& makeinfo --no-split --number-sections igraph.texi \
	&& mv igraph_reference_manual.info igraph.info \
	&& patch <$(top_srcdir)/doc/igraph.info.diff

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
