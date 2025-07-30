# Simple gcc Makefile for reentrant qhull and rbox (default gcc/g++)
#
#   make help
#   See README.txt and ../../Makefile
#       
# Variables
#   DESTDIR        directory for staged installs (GNU Makefile standards)
#   PREFIX         install directory for 'make install' (default /usr/local)
#   BINDIR         directory where to copy executables
#   DOCDIR         directory where to copy html documentation
#   INCDIR         directory where to copy headers
#   LIBDIR         directory where to copy libraries
#   MANDIR         directory where to copy manual pages
#   PCDIR          directory where to copy pkg-config files
#   PRINTMAN       command for printing manual pages
#   PRINTC         command for printing C files
#   CC             ANSI C or C++ compiler
#   CC_OPTS1       options used to compile .c files
#   CC_OPTS2       options used to link .o files
#   CC_OPTS3       options to build shared libraries
#
#   LIBQHULL_OBJS  .o files for linking
#   LIBQHULL_HDRS  .h files for printing
#   CFILES         .c files for printing
#   DOCFILES       documentation files
#   FILES          miscellaneous files for printing
#   TFILES         .txt versions of html files
#   FILES          all other files
#   LIBQHULL_OBJS  specifies the object files of libqhullstatic_r.a
#
# Results
#   rbox           Generates points sets for qhull, qconvex, etc.
#   qhull          Computes convex hulls and related structures
#   qconvex, qdelaunay, qhalf, qvoronoi
#                  Specializations of qhull for each geometric structure
#   libqhullstatic_r.a Static library for reentrant qhull
#   testqset_r     Standalone test of reentrant qset_r.c with mem_r.c
#   user_eg        An example of using qhull (reentrant)
#   user_eg2       An example of using qhull (reentrant)
#
# Make targets
#   make           Build results using gcc or another compiler
#   make all
#   make clean     Remove object files
#   make cleanall  Remove generated files
#   make doc       Print documentation
#   make help
#   make install   Copy results and documentation to BINDIR, DOCDIR, INCDIR, LIBDIR, MANDIR, PCDOC
#   make uninstall Delete Qhull files from BINDIR, DOCDIR, INCDIR, LIBDIR, MANDIR, PCDOC
#   make new       Rebuild qhull and rbox from source
#   make printall  Print all files
#   make qtest     Quick test of qset, rbox, and qhull
#   make test      Quck test of qhull, qconvex, etc.
#
# $Id: //main/2019/qhull/src/libqhull_r/Makefile#13 $

# Do not replace tabs with spaces.  Needed for build rules
# Unix line endings (\n)

PREFIX ?= /usr/local
BINDIR ?= bin
INCDIR ?= include
LIBDIR ?= lib
DOCDIR ?= share/doc/qhull
MANDIR ?= share/man/man1
PCDIR  ?= $(LIBDIR)/pkgconfig

ABS_BINDIR = $(DESTDIR)$(PREFIX)/$(BINDIR)
ABS_INCDIR = $(DESTDIR)$(PREFIX)/$(INCDIR)
ABS_LIBDIR = $(DESTDIR)$(PREFIX)/$(LIBDIR)
ABS_DOCDIR = $(DESTDIR)$(PREFIX)/$(DOCDIR)
ABS_MANDIR = $(DESTDIR)$(PREFIX)/$(MANDIR)
ABS_PCDIR  = $(DESTDIR)$(PREFIX)/$(PCDIR)

qhull_VERSION=$(shell grep 'set.qhull_VERSION ' ../../CMakeLists.txt | grep -o '[0-9.]\+' || echo 0unknown)

# if you do not have enscript, try a2ps or just use lpr.  The files are text.
PRINTMAN = enscript -2rl
PRINTC = enscript -2r
# PRINTMAN = lpr
# PRINTC = lpr

#for Gnu's gcc compiler, -O3 for optimization, -g for debugging, -pg for profiling
# caller may define CC_WARNINGS
# Qhull uses less memory for 32-bit builds on 64-bit hosts
# Enable 32-bit builds with 'make M32=-m32'
# M32     = -m32
# -fpic is required for linking to shared libraries
# -fpic may be slower for 32-bit builds on 64-bit hosts
# Disable -fpic with 'make FPIC=' 
FPIC      = -fpic
CC        = gcc
CC_OPTS1  = -O3 -ansi -I../../src $(CC_WARNINGS) $(M32) $(FPIC)

# for Sun's cc compiler, -fast or O2 for optimization, -g for debugging, -Xc for ANSI
#CC       = cc
#CC_OPTS1 = -Xc -v -fast -I../../src 

# for Silicon Graphics cc compiler, -O2 for optimization, -g for debugging
#CC       = cc
#CC_OPTS1 = -ansi -O2 -I../../src 

# for Next cc compiler with fat executable
#CC       = cc
#CC_OPTS1 = -ansi -O2 -I../../src -arch m68k -arch i386 -arch hppa

# For loader, ld, 
CC_OPTS2 = $(CC_OPTS1)

# Default targets for make

all: qhull_links qhull_all qtest

help:
	head -n 54 Makefile

clean:
	rm -f *.o 
	# Delete linked files from other directories [qhull_links]
	rm -f qconvex_r.c unix_r.c qdelaun_r.c qhalf_r.c qvoronoi_r.c rbox_r.c
	rm -f user_eg_r.c user_eg2_r.c testqset_r.c
	
cleanall: clean
	rm -f qconvex qdelaunay qhalf qvoronoi qhull *.exe
	rm -f core user_eg_r user_eg2_r testqset_r libqhullstatic_r.a

doc: 
	$(PRINTMAN) $(TXTFILES) $(DOCFILES)

install:
	mkdir -p $(ABS_BINDIR)
	mkdir -p $(ABS_DOCDIR)
	mkdir -p $(ABS_DOCDIR)/src
	mkdir -p $(ABS_INCDIR)/libqhull_r
	mkdir -p $(ABS_LIBDIR)
	mkdir -p $(ABS_MANDIR)
	mkdir -p $(ABS_PCDIR)
	cp -p qconvex qdelaunay qhalf qhull qvoronoi rbox $(ABS_BINDIR)
	cp -p libqhullstatic_r.a $(ABS_LIBDIR)
	(cd ../.. && cp -p README.txt REGISTER.txt Announce.txt COPYING.txt index.htm $(ABS_DOCDIR)/)
	(cd ../.. && cp -pr html $(ABS_DOCDIR)/)
	(cd ../.. && cp -p src/Changes.txt $(ABS_DOCDIR)/src/)
	cp -p ../../html/qhull.man $(ABS_MANDIR)/qhull.1
	cp -p ../../html/rbox.man $(ABS_MANDIR)/rbox.1
	cp *.h $(ABS_INCDIR)/libqhull_r
	sed \
		-e 's#@qhull_VERSION@#$(qhull_VERSION)#' \
		-e 's#@CMAKE_INSTALL_PREFIX@#$(PREFIX)#' \
		-e 's#@LIB_INSTALL_DIR@#$(LIBDIR)#' \
		-e 's#@INCLUDE_INSTALL_DIR@#$(INCDIR)#' \
		-e 's#@LIBRARY_NAME@#qhullstatic_r#' \
		-e 's#@LIBRARY_DESCRIPTION@#Qhull reentrant static library#' \
		../../build/qhull.pc.in > $(ABS_PCDIR)/qhullstatic_r.pc

uninstall:
	-(cd $(ABS_BINDIR) && rm -f qconvex qdelaunay qhalf qhull qvoronoi rbox)
	-(cd $(ABS_BINDIR) && rm -f qconvex.exe qdelaunay.exe qhalf.exe qhull.exe qvoronoi.exe rbox.exe)
	-(cd $(ABS_MANDIR) && rm -f qhull.1 rbox.1)
	-(cd $(ABS_DOCDIR) && rm -f README.txt REGISTER.txt Announce.txt COPYING.txt index.htm src/Changes.txt)
	-(cd $(ABS_DOCDIR) && rm -rf html)
	-(cd $(ABS_LIBDIR) && rm -f libqhullstatic_r.a)
	-(cd $(ABS_INCDIR) && rm -rf libqhull_r)
	-(cd $(ABS_PCDIR) && rm -f qhullstatic_r.pc)
	-rmdir $(ABS_DOCDIR)/src
	-rmdir $(ABS_DOCDIR)

new:	cleanall all

printall: doc printh printc printf

printh:
	$(PRINTC) $(LIBQHULL_HDRS)

printc:
	$(PRINTC) $(CFILES)

# LIBQHULL_OBJS_1 ordered by frequency of execution with small files at end.  Better locality.
# Same definitions as ../../Makefile

LIBQHULLS_OBJS_1= global_r.o stat_r.o geom2_r.o poly2_r.o merge_r.o \
	libqhull_r.o geom_r.o poly_r.o qset_r.o mem_r.o random_r.o 

LIBQHULLS_OBJS_2= $(LIBQHULLS_OBJS_1) usermem_r.o userprintf_r.o io_r.o user_r.o

LIBQHULLS_OBJS= $(LIBQHULLS_OBJS_2)  rboxlib_r.o userprintf_rbox_r.o

LIBQHULL_HDRS= user_r.h libqhull_r.h qhull_ra.h geom_r.h \
	io_r.h mem_r.h merge_r.h poly_r.h random_r.h \
	qset_r.h stat_r.h

# CFILES for 'printc', ordered alphabetically after libqhull_r.c 
CFILES= ../qhull/unix_r.c libqhull_r.c geom_r.c geom2_r.c global_r.c io_r.c \
	mem_r.c merge_r.c poly_r.c poly2_r.c random_r.c rboxlib_r.c \
	qset_r.c stat_r.c user_r.c usermem_r.c userprintf_r.c \
	../qconvex/qconvex_r.c ../qdelaunay/qdelaun_r.c \
	../qhalf/qhalf_r.c ../qvoronoi/qvoronoi_r.c

TXTFILES= ../../Announce.txt ../../REGISTER.txt ../../COPYING.txt ../../README.txt ../Changes.txt
DOCFILES= ../../html/rbox.txt ../../html/qhull.txt

.c.o:
	$(CC) -c $(CC_OPTS1) -o $@ $<

# Work around problems with ../ in Red Hat Linux
qhull_links:
	# On MINSYS, 'ln -s' may create a copy instead of a symbolic link
	[ -f qconvex_r.c ]  || ln -s ../qconvex/qconvex_r.c
	[ -f qdelaun_r.c ]  || ln -s ../qdelaunay/qdelaun_r.c
	[ -f qhalf_r.c ]    || ln -s ../qhalf/qhalf_r.c
	[ -f qvoronoi_r.c ] || ln -s ../qvoronoi/qvoronoi_r.c
	[ -f rbox_r.c ]     || ln -s ../rbox/rbox_r.c
	[ -f testqset_r.c ] || ln -s ../testqset_r/testqset_r.c
	[ -f unix_r.c ]     || ln -s ../qhull/unix_r.c
	[ -f user_eg_r.c ]  || ln -s ../user_eg/user_eg_r.c
	[ -f user_eg2_r.c ] || ln -s ../user_eg2/user_eg2_r.c
	# user_eg3_r.c not compiled.  It requires all of libqhullcpp; use qhull/Makefile instead

# compile qhull without using bin/libqhullstatic_r.a
qhull_all: qconvex_r.o qdelaun_r.o qhalf_r.o qvoronoi_r.o unix_r.o user_eg_r.o user_eg2_r.o rbox_r.o testqset_r.o $(LIBQHULLS_OBJS)
	$(CC) -o qconvex $(CC_OPTS2) $(LIBQHULLS_OBJS_2) qconvex_r.o -lm
	$(CC) -o qdelaunay $(CC_OPTS2) $(LIBQHULLS_OBJS_2) qdelaun_r.o -lm
	$(CC) -o qhalf $(CC_OPTS2) $(LIBQHULLS_OBJS_2) qhalf_r.o -lm
	$(CC) -o qvoronoi $(CC_OPTS2) $(LIBQHULLS_OBJS_2) qvoronoi_r.o -lm
	$(CC) -o qhull $(CC_OPTS2) $(LIBQHULLS_OBJS_2) unix_r.o -lm
	$(CC) -o rbox $(CC_OPTS2) $(LIBQHULLS_OBJS) rbox_r.o -lm
	$(CC) -o user_eg $(CC_OPTS2) $(LIBQHULLS_OBJS_2) user_eg_r.o -lm
	$(CC) -o user_eg2 $(CC_OPTS2) $(LIBQHULLS_OBJS_1) user_eg2_r.o  usermem_r.o userprintf_r.o io_r.o -lm
	$(CC) -o testqset_r $(CC_OPTS2) mem_r.o qset_r.o usermem_r.o testqset_r.o -lm
	-ar -rs libqhullstatic_r.a $(LIBQHULLS_OBJS)
	#libqhullstatic_r.a is not needed for qhull
	#If 'ar -rs' fails try using 'ar -s' with 'ranlib'
	#ranlib libqhullstatic_r.a

qtest:
	@echo ============================================
	@echo == make qtest ==============================
	@echo ============================================
	@echo -n "== "
	@date
	@echo
	@echo Testing qset_r.c and mem_r.c with testqset_r
	./testqset_r 10000
	@echo Run the qhull smoketest
	./rbox D4 | ./qhull
	@echo ============================================
	@echo == To smoketest qhull programs
	@echo '==     make test'
	@echo ============================================
	@echo
	@echo ============================================
	@echo == To install qhull or show help
	@echo '==     make help'
	@echo '==     make install'
	@echo ============================================
	@echo

test: qtest
	@echo ==============================
	@echo ========= qconvex ============
	@echo ==============================
	-./rbox 10 | ./qconvex Tv 
	@echo
	@echo ==============================
	@echo ========= qdelaunay ==========
	@echo ==============================
	-./rbox 10 | ./qdelaunay Tv 
	@echo
	@echo ==============================
	@echo ========= qhalf ==============
	@echo ==============================
	-./rbox 10 | ./qconvex FQ FV n Tv | ./qhalf Tv
	@echo
	@echo ==============================
	@echo ========= qvoronoi ===========
	@echo ==============================
	-./rbox 10 | ./qvoronoi Tv
	@echo
	@echo ==============================
	@echo ========= user_eg ============
	@echo == w/o shared library ========
	@echo ==============================
	-./user_eg
	@echo
	@echo ==============================
	@echo ========= user_eg2 ===========
	@echo ==============================
	-./user_eg2
	@echo

# end of Makefile
