## Process this file with automake to produce Makefile.in

DOCINCLUDES = vector.xml error.xml matrix.xml basicigraph.xml generators.xml \
	      structural.xml iterators.xml attributes.xml layout.xml \
	      foreign.xml nongraph.xml

all: doctemplate docfix

doctemplate: $(DOCINCLUDES)

docfix: fdl.xml gpl.xml igraph-docs.xml installation.xml introduction.xml \
	tutorial.xml

vector.xml: vector.xxml ../src/types.h ../src/vector.c
	./doxrox.py -t $< -e c-docbook.re -o $@ ../src/types.h ../src/vector.c

error.xml: error.xxml ../src/error.h
	./doxrox.py -t $< -e c-docbook.re -o $@ ../src/error.h

matrix.xml: matrix.xxml ../src/types.h ../src/matrix.c
	./doxrox.py -t $< -e c-docbook.re -o $@ ../src/types.h ../src/matrix.c

basicigraph.xml: basicigraph.xxml ../src/type_indexededgelist.c
	./doxrox.py -t $< -e c-docbook.re -o $@ ../src/type_indexededgelist.c

generators.xml: generators.xxml ../src/structure_generators.c ../src/games.c
	./doxrox.py -t $< -e c-docbook.re -o $@ \
	../src/structure_generators.c ../src/games.c

structural.xml: structural.xxml ../src/structural_properties.c \
		../src/conversion.c ../src/basic_query.c \
		../src/cocitation.c ../src/components.c
	./doxrox.py -t $< -e c-docbook.re -o $@ \
	../src/structural_properties.c ../src/conversion.c \
	../src/basic_query.c ../src/cocitation.c ../src/components.c

iterators.xml: iterators.xxml ../src/iterators.c
	./doxrox.py -t $< -e c-docbook.re -o $@ ../src/iterators.c

attributes.xml: attributes.xxml ../src/attributes.c
	./doxrox.py -t $< -e c-docbook.re -o $@ ../src/attributes.c

layout.xml: layout.xxml ../src/layout.c
	./doxrox.py -t $< -e c-docbook.re -o $@ ../src/layout.c

foreign.xml: foreign.xxml ../src/foreign.c
	./doxrox.py -t $< -e c-docbook.re -o $@ ../src/foreign.c

nongraph.xml: nongraph.xxml ../src/other.c ../src/random.c
	./doxrox.py -t $< -e c-docbook.re -o $@ ../src/other.c ../src/random.c


CLEANFILES=html/*.html $(DOCINCLUDES)
