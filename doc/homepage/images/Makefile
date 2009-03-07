
all: browsesource-real.png reportbug-real.png searcharchive-real.png \
	subscribe-real.png listbugs-real.png requestfeature-real.png \
	searchdocs-real.png

%.png2: %.svg
	inkscape -e $@ -w 160 $<

%-2.png3: %.png2
	convert $< -negate $@

%-real.png: %.png2 %-2.png3
	montage -geometry +0+0 $^ $@
