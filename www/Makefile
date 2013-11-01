
all: gen

gen: 
	jekyll build

deploy: gen
	rsync -avz --delete _site/ igraph.org:www/

.PHONY: all gen deploy

