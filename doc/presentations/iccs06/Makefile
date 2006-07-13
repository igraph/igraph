SOURCE = iccs06

all: $(SOURCE).pdf

FIGTEX = ~/bin/figtex

FIG = $(wildcard *.fig images/*.fig)
EPS = $(patsubst %.fig, %.eps, $(FIG)) 
EPST = $(patsubst %.fig, %.eps_t, $(FIG))
ALLEPS = $(EPS) 
INCLUDES = $(wildcard *.R)

$(SOURCE).pdf: $(SOURCE).tex $(ALLEPS) $(EPST) $(INCLUDES)
	latex $(SOURCE) && latex $(SOURCE) && \
	dvips -o $(SOURCE).ps -t landscape $(SOURCE).dvi && \
	gs -q -dSAFER -dNOPAUSE -dBATCH -sDEVICE=pdfwrite \
      	-sOutputFile=$(SOURCE).pdf -dCompatibilityLevel=1.3 \
        -dAutoRotatePages=/None \
	-dPDFSETTINGS=/screen \
	$(SOURCE).ps

# Figures

%.eps : %.fig
	$(FIGTEX) $<

clean: 
	rm -f $(SOURCE){.aux,.bbl,.blg,.dvi,.log,.out,.pdf,.ps,.toc}\
	$(EPS) $(EPST)

.PHONY: clean
