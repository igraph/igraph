SOURCE = csardi2

all: $(SOURCE).pdf 

FIGTEX = ~/bin/figtex

FIG = $(wildcard *.fig)
EPS = $(patsubst %.fig, %.eps, $(FIG)) 
EPST = $(patsubst %.fig, %.eps_t, $(FIG))
ALLEPS = $(EPS) 

$(SOURCE).dvi: $(SOURCE).tex $(ALLEPS) $(EPST) net.bib
	latex $(SOURCE) && bibtex $(SOURCE) &&  \
	latex $(SOURCE) && latex $(SOURCE) 

$(SOURCE).pdf: $(SOURCE).dvi
	dvips -o $(SOURCE).ps -t letter $(SOURCE).dvi && \
	ps2pdf $(SOURCE).ps

# Figures

%.eps : %.fig
	$(FIGTEX) $<

clean: 
	rm -f $(SOURCE){.aux,.bbl,.blg,.dvi,.log,.out,.pdf,.ps,.toc}\
	$(EPS) $(EPST)

.PHONY: clean
