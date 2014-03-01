#------------------------------------------------------------------------------
# AMD Makefile (for GNU Make or original make)
#------------------------------------------------------------------------------

VERSION = 2.3.1

default: all

include ../SuiteSparse_config/SuiteSparse_config.mk

demos: all

# Compile all C code.  Do not compile the FORTRAN versions.
all:
	( cd Lib    ; $(MAKE) )
	( cd Demo   ; $(MAKE) )

# compile just the C-callable libraries (not Demos)
library:
	( cd Lib    ; $(MAKE) )

# compile the FORTRAN libraries and demo programs (not compiled by "make all")
fortran:
	( cd Lib    ; $(MAKE) fortran )
	( cd Demo   ; $(MAKE) fortran )

# compile a FORTRAN demo program that calls the C version of AMD
# (not compiled by "make all")
cross:
	( cd Demo   ; $(MAKE) cross )

# remove object files, but keep the compiled programs and library archives
clean:
	( cd Lib    ; $(MAKE) clean )
	( cd Demo   ; $(MAKE) clean )
	( cd MATLAB ; $(RM) $(CLEAN) )
	( cd Doc    ; $(MAKE) clean )

# clean, and then remove compiled programs and library archives
purge:
	( cd Lib    ; $(MAKE) purge )
	( cd Demo   ; $(MAKE) purge )
	( cd MATLAB ; $(RM) $(CLEAN) ; $(RM) *.mex* )
	( cd Doc    ; $(MAKE) purge )

distclean: purge

# create PDF documents for the original distribution
docs:
	( cd Doc    ; $(MAKE) )

# get ready for distribution
dist: purge
	( cd Demo   ; $(MAKE) dist )
	( cd Doc    ; $(MAKE) )

ccode: library

lib: library

# install AMD
install:
	$(CP) Lib/libamd.a $(INSTALL_LIB)/libamd.$(VERSION).a
	( cd $(INSTALL_LIB) ; ln -sf libamd.$(VERSION).a libamd.a )
	$(CP) Include/amd.h $(INSTALL_INCLUDE)
	chmod 644 $(INSTALL_LIB)/libamd*
	chmod 644 $(INSTALL_INCLUDE)/amd.h

# uninstall AMD
uninstall:
	$(RM) $(INSTALL_LIB)/libamd*.a
	$(RM) $(INSTALL_INCLUDE)/amd.h

