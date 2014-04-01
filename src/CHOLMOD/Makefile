#-------------------------------------------------------------------------------
# CHOLMOD Makefile
#-------------------------------------------------------------------------------

VERSION = 2.1.2

# Note: If you do not have METIS, or do not wish to use it in CHOLMOD, you must
# compile CHOLMOD with the -DNPARTITION flag.
# See ../SuiteSparse_config/SuiteSparse_config.mk .

default: all

include ../SuiteSparse_config/SuiteSparse_config.mk

# Compile the C-callable libraries and the Demo programs.
all:
	( cd Demo ; $(MAKE) )

# Compile the C-callable libraries only.
library:
	( cd Lib ; $(MAKE) )

# Remove all files not in the original distribution
purge:
	( cd Tcov ; $(MAKE) purge )
	( cd Lib ; $(MAKE) purge )
	( cd Valgrind ; $(MAKE) dopurge )
	( cd Demo ; $(MAKE) purge )
	( cd Doc ; $(MAKE) purge )
	( cd MATLAB ; $(RM) $(CLEAN) rename.h *.mex* )

# Remove all files not in the original distribution, except keep the 
# compiled libraries.
clean:
	( cd Tcov ; $(MAKE) clean )
	( cd Lib ; $(MAKE) clean )
	( cd Valgrind ; $(MAKE) clean )
	( cd Demo ; $(MAKE) clean )
	( cd MATLAB ; $(RM) $(CLEAN) )

distclean: purge

ccode: all

# Run the test coverage suite.  Takes about 40 minutes on a 3.2GHz Pentium.
# Requires Linux (gcc, gcov).
cov:
	( cd Tcov ; $(MAKE) )

# Run the test coverage suite using Valgrind.  This takes a *** long *** time.
valgrind:
	( cd Valgrind ; $(MAKE) )

# Compile the C-callable libraries and the Demo programs.
demos:
	( cd Demo ; $(MAKE) )

# create PDF documents for the original distribution
docs:
	( cd Doc    ; $(MAKE) )

# install CHOLMOD
install:
	$(CP) Lib/libcholmod.a $(INSTALL_LIB)/libcholmod.$(VERSION).a
	( cd $(INSTALL_LIB) ; ln -sf libcholmod.$(VERSION).a libcholmod.a )
	$(CP) Include/cholmod*.h $(INSTALL_INCLUDE)
	$(RM) $(INSTALL_INCLUDE)/cholmod_internal.h
	chmod 644 $(INSTALL_LIB)/libcholmod*.a
	chmod 644 $(INSTALL_INCLUDE)/cholmod*.h

# uninstall CHOLMOD
uninstall:
	$(RM) $(INSTALL_LIB)/libcholmod*.a
	$(RM) $(INSTALL_INCLUDE)/cholmod*.h

