
SOURCE=csardi

all: csardi-final.pdf

ALLSVG = arch.pdf orgnet.pdf er.pdf ba.pdf attr.pdf \
	 plot.pdf indd.pdf outdd.pdf schema.pdf igraph1.pdf \
	 igraph2.pdf igraph3.pdf igraph11.pdf igraph0.pdf \
	 mixed.pdf hyper.pdf bipartite.pdf transitivity.pdf \
	 g1.pdf g2.pdf g3.pdf g4.pdf g5.pdf g6.pdf g7.pdf \
	 g8.pdf

OTHER = igraph.pdf

ALLPNG = jgplot.png rglplot.png source_c.png source_cpp.png \
         Rlogo-2.png Libxml2-Logo-180x168.png gmplogo2.png \
	 hist3d_2lights.jpg logo125.png SQLite.png hash_small.png \
	 matrix.jpg tkplot.png

csardi.pdf: csardi.tex $(ALLSVG) $(ALLPNG) $(OTHER)
	pdflatex csardi.tex

csardi-final.pdf: csardi.pdf
	ppower4 csardi.pdf csardi-final.pdf

igraph.pdf : igraph.svg
	inkscape $< -A $@

%.pdf : %.svg.gz
	gzip -dc $< > $(*F).svg
	inkscape $(*F).svg -A $@
	rm $(*F).svg

clean : 
	rm -rf $(ALLSVG) igraph.pdf csardi-final.pdf csardi.aux \
	csardi.log csardi.pdf

.PHONY: clean
