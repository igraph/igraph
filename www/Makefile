
all: nothing

HTML= index.html getstarted.html news.html \
      _layouts/default.html _layouts/newstemp.html _layouts/manual.html

CSS= css/affix.css css/manual.css css/other.css fonts/fonts.css

../doc/jekyll/stamp: ../doc/html/stamp
	cd ../doc && make jekyll

doc/c/stamp: ../doc/jekyll/stamp
	rm -rf doc/c
	mkdir -p doc
	cp -r ../doc/jekyll doc/c

stamp: $(HTML) $(CSS) doc/c/stamp
	jekyll build
	touch stamp

deploy: stamp
	rsync -avz --delete _site/ igraph.org:www/

.PHONY: all deploy

