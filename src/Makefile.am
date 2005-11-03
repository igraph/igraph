lib_LTLIBRARIES         = libigraph.la

libigraph_la_SOURCES       = basic_query.c games.c cocitation.c iterators.c \
			     structural_properties.c components.c layout.c \
			     structure_generators.c conversion.c \
			     measure_dynamics.c type_indexededgelist.c \
			     error.c other.c types.c foreign.c random.c \
			     attributes.c \
			     memory.h random.h \
			     foreign-ncol-parser.y foreign-ncol-lexer.l
libigraph_la_LDFLAGS       = -no-undefined
libigraph_la_LIBADD        = -lm

includedir = $(prefix)/include/igraph
include_HEADERS = igraph.h types.h attributes.h
