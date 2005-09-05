
CC=cc
CFLAGS=-O6
CPPFLAGS=
LDFLAGS=-lm

all: library

include depend

library: libigraph.so.0.0.1

LIBRARY_SOURCES=src/basic_query.c \
		src/games.c \
		src/cocitation.c \
		src/iterators.c \
		src/structural_properties.c \
		src/components.c \
		src/layout.c \
		src/structure_generators.c \
		src/conversion.c \
		src/measure_dynamics.c \
		src/type_indexededgelist.c \
		src/error.c \
		src/other.c \
		src/types.c \
		src/foreign.c \
		src/random.c

LIBRARY_OBJECTS=$(patsubst %.c, %.o, $(LIBRARY_SOURCES))

libigraph.so.0.0.1: $(LIBRARY_OBJECTS)
	$(CC) -shared -Wl,-soname,libigraph.so.0 \
	-o libigraph.so.0.0.1 $(LIBRARY_OBJECTS) $(LDFLAGS)

%.o : %.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -c -fPIC $< -o $@

clean: 
	rm -f libigraph.so.0.0.1 src/*.o src/igraph.so depend

.PHONY: clean

depend: 
	$(CC) -MM $(LIBRARY_SOURCES) \
	$(CPPFLAGS) $(CFLAGS) >depend
