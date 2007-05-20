
all: ps 

ps: netsci07.ps

pdf: netsci07.pdf 

netsci07.ps: netsci07.dvi
	dvips -t a0 netsci07.dvi -o netsci07.ps

netsci07.pdf: netsci07.ps
	ps2pdf netsci07.ps

IMAGES=karate3d.eps arch.eps arch.eps_t RMKI-BME.eps CCSS.eps python3d.eps \
	sna_screenshot.eps tkplot.eps diameter.eps degreedist.eps python3d-2.eps

karate3d.eps: karate3d.png
	convert karate3d.png karate3d.eps

arch.eps arch.eps_t: arch.svg.gz
	gzip -dc arch.svg.gz >arch.svg && \
	svgtex arch.svg

RMKI-BME.eps: RMKI-BME.svg.gz
	gzip -dc RMKI-BME.svg.gz >RMKI-BME.svg && \
	inkscape -E RMKI-BME.eps -B RMKI-BME.svg

CCSS.eps: CCSS.svg.gz
	gzip -dc CCSS.svg.gz >CCSS.svg && \
	inkscape -E CCSS.eps -B CCSS.svg

python3d.eps: python3d.png
	convert python3d.png python3d.eps

python3d-2.eps: python3d-2.png
	convert python3d-2.png python3d-2.eps

sna_screenshot.eps: sna_screenshot.jpg
	convert sna_screenshot.jpg sna_screenshot.eps

tkplot.eps: tkplot.png
	convert tkplot.png tkplot.eps

diameter.eps: diameter.png
	convert diameter.png diameter.eps

degreedist.eps: degreedist.png
	convert degreedist.png degreedist.eps

netsci07.dvi: netsci07.tex $(IMAGES)
	latex netsci07.tex

clean:
	rm -f karate3d.eps screenshots.eps arch.eps arch.eps_t \
	RMKI-BME.eps CCSS.eps netsci07.pdf netsci07.ps \
	python3d-2.eps sna_screenshot.eps \
	netsci07.aux netsci07.dvi netsci07.log \
	screenshots.svg arch.svg RMKI-BME.svg CCSS.svg


.PHONY: all ps pdf clean
