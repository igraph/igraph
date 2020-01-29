## Process this file with automake to produce Makefile.in

DOXROX=$(top_srcdir)/doc/doxrox.py
REGEX=$(top_srcdir)/doc/c-docbook.re
SRCDIR=$(top_srcdir)/src
INCLUDEDIR=$(top_srcdir)/include
BUILD_INCLUDEDIR=$(top_builddir)/include
REMOVEEXAMPLES=$(top_srcdir)/tools/removeexamples.py

DOCINCLUDES = vector.xml error.xml matrix.xml basicigraph.xml generators.xml \
	      structural.xml iterators.xml attributes.xml layout.xml \
	      foreign.xml nongraph.xml isomorphism.xml motifs.xml \
	      operators.xml flows.xml community.xml cliques.xml \
	      sparsematrix.xml stack.xml dqueue.xml heap.xml strvector.xml \
	      adjlist.xml arpack.xml bipartite.xml visitors.xml random.xml \
	      separators.xml memory.xml sparsemat.xml hrg.xml \
	      scg.xml spatialgames.xml threading.xml progress.xml status.xml \
	      graphlets.xml embedding.xml coloring.xml

DOCFIX = fdl.xml gpl.xml installation.xml introduction.xml \
	 tutorial.xml licenses.xml pmt.xml

DOCFIX2 = $(patsubst %,$(top_srcdir)/doc/%,$(DOCFIX))

all: doctemplate

html: html/stamp

jekyll: jekyll/stamp

info: igraph.info

dvi: igraph-docs.dvi

pdf: igraph-docs.pdf

ps: igraph-docs.ps

doctemplate: $(DOCINCLUDES)

tags: $(DOCINCLUDES) $(DOCFIX)
	cat $(DOCINCLUDES) $(DOCFIX)| grep 'id="[^-"]*">' | sed 's/.*id="\([^"]*\)">.*/\1@@/' | tr '@' '\t' | sort > tags

EXAMPLES = \
    $(patsubst $(top_srcdir)/%.c,$(top_builddir)/%.c.xml,$(wildcard $(top_srcdir)/examples/simple/*.c)) \
    $(patsubst $(top_srcdir)/%.c,$(top_builddir)/%.c.xml,$(wildcard $(top_srcdir)/examples/tests/*.c))

$(top_builddir)/examples/simple/%.c.xml: $(top_srcdir)/examples/simple/%.c tags
	if [ ! -d "$(top_builddir)/examples/simple" -a "x$(top_srcdir)" != "x$(top_builddir)" ]; then mkdir -p $(top_builddir)/examples/simple; fi
	source-highlight --src-lang c --out-format docbook --input $< --output $@ --gen-references=inline --ctags="" --outlang-def=$(top_srcdir)/doc/docbook.outlang

$(top_builddir)/examples/tests/%.c.xml: $(top_srcdir)/examples/tests/%.c tags
	if [ ! -d "$(top_builddir)/examples/tests" -a "x$(top_srcdir)" != "x$(top_builddir)" ]; then mkdir -p $(top_builddir)/examples/tests; fi
	source-highlight --src-lang c --out-format docbook --input $< --output $@ --gen-references=inline --ctags="" --outlang-def=$(top_srcdir)/doc/docbook.outlang

igraph-docs.xml: $(DOCINCLUDES) $(EXAMPLES)
	if [ "x$(top_srcdir)" != "x$(top_builddir)" ]; then cp $(top_srcdir)/doc/igraph-docs.xml . ; fi
	touch igraph-docs.xml

igraph-docs0.xml: igraph-docs.xml
	sed "s/@IGRAPH_VERSION@/`cat $(top_srcdir)/IGRAPH_VERSION`/g" $(top_srcdir)/doc/igraph-docs.xml > igraph-docs0.xml

vector.xml: vector.xxml $(INCLUDEDIR)/igraph_vector_pmt.h $(SRCDIR)/vector.pmt $(SRCDIR)/vector_ptr.c $(SRCDIR)/vector.c $(INCLUDEDIR)/igraph_vector_ptr.h
	$(DOXROX) -t $< -e $(REGEX) -o $@ $(INCLUDEDIR)/igraph_vector_pmt.h \
	$(SRCDIR)/vector.pmt $(SRCDIR)/vector_ptr.c $(SRCDIR)/vector.c \
	$(INCLUDEDIR)/igraph_vector_ptr.h

error.xml: error.xxml $(INCLUDEDIR)/igraph_error.h $(SRCDIR)/igraph_error.c
	$(DOXROX) -t $< -e $(REGEX) -o $@ $(INCLUDEDIR)/igraph_error.h $(SRCDIR)/igraph_error.c

matrix.xml: matrix.xxml $(INCLUDEDIR)/igraph_matrix.h $(SRCDIR)/matrix.pmt
	$(DOXROX) -t $< -e $(REGEX) -o $@ $(INCLUDEDIR)/igraph_matrix.h $(SRCDIR)/matrix.pmt

sparsematrix.xml: sparsematrix.xxml $(INCLUDEDIR)/igraph_spmatrix.h $(SRCDIR)/spmatrix.c
	$(DOXROX) -t $< -e $(REGEX) -o $@ $(INCLUDEDIR)/igraph_spmatrix.h $(SRCDIR)/spmatrix.c

sparsemat.xml: sparsemat.xxml $(SRCDIR)/sparsemat.c
	$(DOXROX) -t $< -e $(REGEX) -o $@ $(SRCDIR)/sparsemat.c

hrg.xml: hrg.xxml $(SRCDIR)/igraph_hrg.cc $(INCLUDEDIR)/igraph_hrg.h
	$(DOXROX) -t $< -e $(REGEX) -o $@ $(SRCDIR)/igraph_hrg.cc $(INCLUDEDIR)/igraph_hrg.h

scg.xml: scg.xxml $(SRCDIR)/scg.c
	$(DOXROX) -t $< -e $(REGEX) -o $@ $(SRCDIR)/scg.c

embedding.xml: embedding.xxml $(SRCDIR)/embedding.c
	$(DOXROX) -t $< -e $(REGEX) -o $@ $(SRCDIR)/embedding.c

basicigraph.xml: basicigraph.xxml $(SRCDIR)/type_indexededgelist.c $(SRCDIR)/structural_properties.c
	$(DOXROX) -t $< -e $(REGEX) -o $@ $(SRCDIR)/type_indexededgelist.c $(SRCDIR)/structural_properties.c

generators.xml: generators.xxml \
		$(SRCDIR)/atlas.c \
		$(SRCDIR)/forestfire.c \
		$(SRCDIR)/games.c \
		$(SRCDIR)/structure_generators.c \
		$(SRCDIR)/structural_properties.c \
		$(SRCDIR)/sbm.c \
		$(SRCDIR)/dotproduct.c \
		$(SRCDIR)/degree_sequence.cpp
	$(DOXROX) -t $< -e $(REGEX) -o $@ \
	$(SRCDIR)/atlas.c \
	$(SRCDIR)/forestfire.c \
	$(SRCDIR)/games.c \
	$(SRCDIR)/structure_generators.c \
	$(SRCDIR)/structural_properties.c \
	$(SRCDIR)/sbm.c \
	$(SRCDIR)/dotproduct.c \
	$(SRCDIR)/degree_sequence.cpp

structural.xml: structural.xxml $(SRCDIR)/structural_properties.c \
		$(SRCDIR)/spanning_trees.c \
		$(SRCDIR)/conversion.c $(SRCDIR)/basic_query.c \
		$(SRCDIR)/cocitation.c $(SRCDIR)/components.c \
		$(SRCDIR)/spectral_properties.c $(SRCDIR)/cores.c \
		$(SRCDIR)/centrality.c $(SRCDIR)/decomposition.c \
		$(SRCDIR)/mixing.c $(INCLUDEDIR)/igraph_arpack.h \
		$(SRCDIR)/distances.c $(SRCDIR)/feedback_arc_set.c \
		$(SRCDIR)/matching.c $(SRCDIR)/triangles.c \
		$(SRCDIR)/paths.c $(INCLUDEDIR)/igraph_centrality.h
	$(DOXROX) -t $< -e $(REGEX) -o $@ \
	$(SRCDIR)/structural_properties.c $(SRCDIR)/spanning_trees.c \
	$(SRCDIR)/conversion.c $(SRCDIR)/basic_query.c $(SRCDIR)/cocitation.c \
	$(SRCDIR)/components.c $(SRCDIR)/spectral_properties.c $(SRCDIR)/cores.c \
	$(SRCDIR)/centrality.c $(SRCDIR)/decomposition.c $(SRCDIR)/mixing.c \
	$(INCLUDEDIR)/igraph_arpack.h  $(SRCDIR)/distances.c \
	$(SRCDIR)/feedback_arc_set.c $(SRCDIR)/matching.c $(SRCDIR)/triangles.c \
	$(SRCDIR)/paths.c $(INCLUDEDIR)/igraph_centrality.h \
	$(SRCDIR)/scan.c


iterators.xml: iterators.xxml $(SRCDIR)/iterators.c $(INCLUDEDIR)/igraph_iterators.h
	$(DOXROX) -c -t $< -e $(REGEX) -o $@ $(SRCDIR)/iterators.c \
	$(INCLUDEDIR)/igraph_iterators.h

attributes.xml: attributes.xxml $(SRCDIR)/attributes.c $(INCLUDEDIR)/igraph_attributes.h $(SRCDIR)/cattributes.c
	$(DOXROX) -t $< -e $(REGEX) -o $@ $(SRCDIR)/attributes.c \
	$(INCLUDEDIR)/igraph_attributes.h $(SRCDIR)/cattributes.c

layout.xml: layout.xxml $(SRCDIR)/layout.c $(INCLUDEDIR)/igraph_layout.h $(SRCDIR)/drl_layout.cpp $(SRCDIR)/drl_layout_3d.cpp $(SRCDIR)/sugiyama.c $(SRCDIR)/layout_fr.c $(SRCDIR)/layout_kk.c $(SRCDIR)/layout_gem.c $(SRCDIR)/layout_dh.c
	$(DOXROX) -t $< -e $(REGEX) -o $@ $(SRCDIR)/layout.c $(INCLUDEDIR)/igraph_layout.h $(SRCDIR)/drl_layout.cpp $(SRCDIR)/drl_layout_3d.cpp $(SRCDIR)/sugiyama.c $(SRCDIR)/layout_fr.c $(SRCDIR)/layout_kk.c $(SRCDIR)/layout_gem.c $(SRCDIR)/layout_dh.c

foreign.xml: foreign.xxml $(SRCDIR)/foreign.c $(SRCDIR)/foreign-graphml.c
	$(DOXROX) -t $< -e $(REGEX) -o $@ $(SRCDIR)/foreign.c \
	$(SRCDIR)/foreign-graphml.c

nongraph.xml: nongraph.xxml $(SRCDIR)/other.c $(SRCDIR)/random.c $(SRCDIR)/version.c $(INCLUDEDIR)/igraph_nongraph.h
	$(DOXROX) -t $< -e $(REGEX) -o $@ $(INCLUDEDIR)/igraph_nongraph.h $(SRCDIR)/other.c $(SRCDIR)/random.c $(SRCDIR)/version.c $(SRCDIR)/dotproduct.c

isomorphism.xml: isomorphism.xxml $(SRCDIR)/topology.c $(INCLUDEDIR)/igraph_topology.h $(SRCDIR)/bliss.cc $(SRCDIR)/lad.c
	$(DOXROX) -c -t $< -e $(REGEX) -o $@ $(SRCDIR)/topology.c $(INCLUDEDIR)/igraph_topology.h $(SRCDIR)/bliss.cc $(SRCDIR)/lad.c

coloring.xml: coloring.xxml $(SRCDIR)/coloring.c $(INCLUDEDIR)/igraph_coloring.h
	$(DOXROX) -t $< -e $(REGEX) -o $@ $(SRCDIR)/coloring.c $(INCLUDEDIR)/igraph_coloring.h

motifs.xml: motifs.xxml $(INCLUDEDIR)/igraph_motifs.h $(SRCDIR)/motifs.c \
		$(SRCDIR)/triangles.c
	$(DOXROX) -t $< -e $(REGEX) -o $@ $(INCLUDEDIR)/igraph_motifs.h \
		$(SRCDIR)/motifs.c $(SRCDIR)/triangles.c

operators.xml: operators.xxml $(SRCDIR)/operators.c
	$(DOXROX) -t $< -e $(REGEX) -o $@ $(SRCDIR)/operators.c

flows.xml: flows.xxml $(SRCDIR)/flow.c $(SRCDIR)/st-cuts.c $(SRCDIR)/cohesive_blocks.c $(INCLUDEDIR)/igraph_flow.h
	$(DOXROX) -t $< -e $(REGEX) -o $@ $(SRCDIR)/flow.c $(SRCDIR)/st-cuts.c $(SRCDIR)/cohesive_blocks.c $(INCLUDEDIR)/igraph_flow.h

community.xml: community.xxml $(SRCDIR)/community.c $(SRCDIR)/clustertool.cpp $(SRCDIR)/walktrap.cpp $(SRCDIR)/community_leiden.c $(SRCDIR)/fast_community.c $(SRCDIR)/optimal_modularity.c $(INCLUDEDIR)/igraph_community.h $(SRCDIR)/infomap.cc
	$(DOXROX) -t $< -e $(REGEX) -o $@ $(SRCDIR)/community.c $(SRCDIR)/clustertool.cpp $(SRCDIR)/walktrap.cpp $(SRCDIR)/community_leiden.c $(SRCDIR)/fast_community.c $(SRCDIR)/optimal_modularity.c $(INCLUDEDIR)/igraph_community.h $(SRCDIR)/infomap.cc

cliques.xml: cliques.xxml $(SRCDIR)/cliques.c $(SRCDIR)/maximal_cliques.c
	$(DOXROX) -t $< -e $(REGEX) -o $@ $(SRCDIR)/cliques.c $(SRCDIR)/maximal_cliques.c $(INCLUDEDIR)/igraph_cliques.h

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

arpack.xml: arpack.xxml $(INCLUDEDIR)/igraph_arpack.h $(SRCDIR)/arpack.c $(INCLUDEDIR)/igraph_blas.h $(SRCDIR)/blas.c $(INCLUDEDIR)/igraph_lapack.h $(SRCDIR)/lapack.c
	$(DOXROX) -t $< -e $(REGEX) -o $@ $(INCLUDEDIR)/igraph_arpack.h $(SRCDIR)/arpack.c $(INCLUDEDIR)/igraph_blas.h $(SRCDIR)/blas.c $(INCLUDEDIR)/igraph_lapack.h $(SRCDIR)/lapack.c

bipartite.xml: bipartite.xxml $(SRCDIR)/bipartite.c
	$(DOXROX) -t $< -e $(REGEX) -o $@ $(SRCDIR)/bipartite.c

visitors.xml: visitors.xxml $(SRCDIR)/visitors.c  $(SRCDIR)/random_walk.c $(INCLUDEDIR)/igraph_visitor.h
	$(DOXROX) -t $< -e $(REGEX) -o $@ $(SRCDIR)/visitors.c  $(SRCDIR)/random_walk.c $(INCLUDEDIR)/igraph_visitor.h

random.xml: random.xxml $(SRCDIR)/random.c
	$(DOXROX) -t $< -e $(REGEX) -o $@ $(SRCDIR)/random.c

separators.xml: separators.xxml $(SRCDIR)/separators.c $(SRCDIR)/st-cuts.c
	$(DOXROX) -t $< -e $(REGEX) -o $@ $(SRCDIR)/separators.c \
		$(SRCDIR)/st-cuts.c

memory.xml: memory.xxml $(SRCDIR)/memory.c
	$(DOXROX) -t $< -e $(REGEX) -o $@ $(SRCDIR)/memory.c

spatialgames.xml: spatialgames.xxml \
		   $(SRCDIR)/microscopic_update.c $(SRCDIR)/sir.c $(INCLUDEDIR)/igraph_epidemics.h
	$(DOXROX) -t $< -e $(REGEX) -o $@ \
	$(SRCDIR)/microscopic_update.c $(SRCDIR)/sir.c $(INCLUDEDIR)/igraph_epidemics.h

threading.xml: threading.xxml $(BUILD_INCLUDEDIR)/igraph_threading.h
	$(DOXROX) -t $< -e $(REGEX) -o $@ $(BUILD_INCLUDEDIR)/igraph_threading.h

progress.xml: progress.xxml $(INCLUDEDIR)/igraph_progress.h $(SRCDIR)/progress.c
	$(DOXROX) -t $< -e $(REGEX) -o $@ $(INCLUDEDIR)/igraph_progress.h $(SRCDIR)/progress.c

status.xml: status.xxml $(INCLUDEDIR)/igraph_statusbar.h $(SRCDIR)/statusbar.c
	$(DOXROX) -t $< -e $(REGEX) -o $@ $(INCLUDEDIR)/igraph_statusbar.h $(SRCDIR)/statusbar.c

graphlets.xml: graphlets.xxml $(SRCDIR)/glet.c
	$(DOXROX) -t $< -e $(REGEX) -o $@ $(SRCDIR)/glet.c

html/stamp: igraph-docs0.xml $(DOCFIX2) gtk-doc.xsl
	if [ "x$(top_srcdir)" != "x$(top_builddir)" ]; then cp $(DOCFIX2) . ; fi && \
	xmlto -x $(top_srcdir)/doc/gtk-doc.xsl -o html xhtml igraph-docs0.xml \
	&& touch html/stamp

jekyll/stamp: html/stamp
	rm -rf jekyll && mkdir jekyll && cp html/* jekyll/ && rm jekyll/stamp
	cd html && for i in *.html; do                        \
	  cat $$i | ../../tools/extract_body.sh               \
		  | ../../tools/protect_braces.sh             \
		  | ../../tools/jekyll_header.sh              \
                  > ../jekyll/$$i; done
	touch jekyll/stamp

igraph.info0.diff: igraph.info.diff
	sed "s/@IGRAPH_VERSION@/`cat $(top_srcdir)/IGRAPH_VERSION`/g" $(top_srcdir)/doc/igraph.info.diff > igraph.info0.diff

igraph.info: igraph-docs0.xml $(DOCFIX2) igraph.info0.diff
	if [ "x$(top_srcdir)" != "x$(top_builddir)" ]; then cp $(DOCFIX2) . ; fi && \
	cat tutorial.xml | sed '/<link/{N;s/<link\n*[^>]*>//g;}' | \
	  sed 's/<\/link>//g' >tutorial-info.xml \
	&& sed 's/tutorial\.xml/tutorial-info.xml/' igraph-docs0.xml > igraph-docs-info.xml \
	&& xmllint --xinclude -o igraph-docs-info2.xml igraph-docs-info.xml \
	&& $(REMOVEEXAMPLES) igraph-docs-info2.xml igraph-docs-info3.xml \
	&& sed '/<book/q' igraph-docs-info2.xml | grep -v '<book' >igraph-docs-info4.xml\
	&& cat igraph-docs-info3.xml >> igraph-docs-info4.xml \
	&& db2x_xsltproc --xinclude -s texi igraph-docs-info4.xml -o igraph.txml \
	&& db2x_texixml --encoding=utf8 igraph.txml \
	&& mv igraph_reference_manual.texi igraph.texi \
	&& makeinfo --no-split --number-sections igraph.texi \
	&& mv igraph_reference_manual.info igraph.info \
	&& patch < igraph.info0.diff

igraph-docs.dvi: igraph-docs0.xml $(DOCFIX2)
	if [ "x$(top_srcdir)" != "x$(top_builddir)" ]; then cp $(DOCFIX2) . ; fi && \
	xmllint --xinclude -o igraph-docs2.xml igraph-docs0.xml \
	&& $(REMOVEEXAMPLES) igraph-docs2.xml igraph-docs3.xml \
	&& sed '/<book/q' igraph-docs2.xml | grep -v '<book' >igraph-docs4.xml\
	&& cat igraph-docs3.xml >> igraph-docs4.xml \
	&& env < /dev/null hash_extra=650000 docbook2dvi igraph-docs4.xml \
	&& mv igraph-docs4.dvi igraph-docs.dvi

igraph-docs.ps: igraph-docs0.xml $(DOCFIX2)
	if [ "x$(top_srcdir)" != "x$(top_builddir)" ]; then cp $(DOCFIX2) . ; fi && \
	xmllint --xinclude -o igraph-docs2.xml igraph-docs0.xml \
	&& $(REMOVEEXAMPLES) igraph-docs2.xml igraph-docs3.xml \
	&& sed '/<book/q' igraph-docs2.xml | grep -v '<book' >igraph-docs4.xml\
	&& cat igraph-docs3.xml >> igraph-docs4.xml \
	&& env < /dev/null hash_extra=65000 docbook2s igraph-docs4.xml \
	&& mv igraph-docs4.ps igraph-docs.ps

igraph-docs.pdf: igraph-docs0.xml $(DOCFIX2)
	if [ "x$(top_srcdir)" != "x$(top_builddir)" ]; then cp $(DOCFIX2) . ; fi && \
	xmllint --xinclude -o igraph-docs2.xml igraph-docs0.xml \
	&& $(REMOVEEXAMPLES) igraph-docs2.xml igraph-docs3.xml \
	&& sed '/<book/q' igraph-docs2.xml | grep -v '<book' >igraph-docs4.xml\
	&& cat igraph-docs3.xml >> igraph-docs4.xml \
	&& xsltproc http://docbook.sourceforge.net/release/xsl/current/fo/docbook.xsl \
	   igraph-docs4.xml > igraph-docs.fo \
	&& fop -fo igraph-docs.fo -pdf igraph-docs.pdf

CLEANFILES=$(DOCINCLUDES) html/*.html html/stamp \
	igraph-docs.{dvi,info,pdf,ps,texi,txml} igraph-docs2.xml igraph-docs0.xml
