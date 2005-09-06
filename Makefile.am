SUBDIRS=src 

EXTRA_DIST=include/*

doc: 
	doxygen Doxyfile

CLEANFILES=doc/*

## to make sure make doc will always call doxygen rule
.PHONY: doc
