## Process this file with automake to produce Makefile.in

DOCFIX = datamodel.xml igraph.xml installation.xml introduction.xml \
	 short.xml structures.xml ../licenses.xml creating.xml graph.xml \
	 random.xml bibliography.xml ../bibdatabase.xml foreign.xml \
	 operators.xml motifs.xml pathlength.xml centrality.xml citation.xml \
	 components.xml attributes.xml iterators.xml plotting.xml

DOCFIX2 = $(patsubst %,$(top_srcdir)/doc/book/%,$(DOCFIX))

all: doctemplate

html: html/stamp

info: igraph.info

dvi: igraph.dvi

pdf: igraph.pdf

ps: igraph.ps

doctemplate: $(DOCINCLUDES)

igraph.xml: $(DOCINCLUDES)
	touch igraph.xml

html/stamp: igraph.xml $(DOCFIX2) gtk-doc.xsl
	if [ "x$(top_srcdir)" != "x$(top_builddir)" ]; then cp $(DOCFIX2) . ; fi && \
	xmlto -x $(top_srcdir)/doc/book/gtk-doc.xsl -o html xhtml igraph.xml \
	&& touch html/stamp

igraph.info: igraph.xml $(DOCFIX2)
	if [ "x$(top_srcdir)" != "x$(top_builddir)" ]; then cp $(DOCFIX2) . ; fi && \
	db2x_xsltproc --xinclude -s texi igraph.xml -o igraph.txml \
	&& db2x_texixml igraph.txml \
	&& makeinfo --no-split igraph.texi

igraph.dvi: igraph.xml $(DOCFIX2)
	if [ "x$(top_srcdir)" != "x$(top_builddir)" ]; then cp $(DOCFIX2) . ; fi && \
	xmllint --xinclude -o igraph2.xml igraph.xml \
	&& docbook2dvi igraph2.xml \
	&& mv igraph2.dvi igraph.dvi

igraph.ps: igraph.xml $(DOCFIX2)
	if [ "x$(top_srcdir)" != "x$(top_builddir)" ]; then cp $(DOCFIX2) . ; fi && \
	xmllint --xinclude -o igraph2.xml igraph.xml \
	&& docbook2ps igraph2.xml \
	&& mv igraph2.ps igraph.ps

igraph.pdf: igraph.xml $(DOCFIX2)
	if [ "x$(top_srcdir)" != "x$(top_builddir)" ]; then cp $(DOCFIX2) . ; fi && \
	xmllint --xinclude -o igraph2.xml igraph.xml \
	&& docbook2pdf igraph2.xml \
	&& mv igraph2.pdf igraph.pdf

CLEANFILES=$(DOCINCLUDES) html/*.html html/stamp \
	igraph.{dvi,info,pdf,ps,texi,txml} igraph2.xml 
