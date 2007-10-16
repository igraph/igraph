if USE_SHELL
   MAYBE_SHELL_INTERFACE = shell
endif

SUBDIRS = $(MAYBE_SHELL_INTERFACE)

MAINTAINERCLEANFILES = Makefile.in

DISTCLEANFILES = R/igraph/src/*.o R/igraph/src/igraph.so \
                 shell/.deps
