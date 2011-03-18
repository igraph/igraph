# Unix makefile: see README.
# For C++, first "make hadd".
# If your compiler does not recognize ANSI C, add
#	-DKR_headers
# to the CFLAGS = line below.
# On Sun and other BSD systems that do not provide an ANSI sprintf, add
#	-DUSE_STRLEN
# to the CFLAGS = line below.
# On Linux systems, add
#	-DNON_UNIX_STDIO
# to the CFLAGS = line below.  For libf2c.so under Linux, also add
#	-fPIC
# to the CFLAGS = line below.

.SUFFIXES: .c .o
CC = cc
SHELL = /bin/sh
CFLAGS = -O

# compile, then strip unnecessary symbols
.c.o:
	$(CC) -c -DSkip_f2c_Undefs $(CFLAGS) $*.c
	ld -r -x -o $*.xxx $*.o
	mv $*.xxx $*.o
## Under Solaris (and other systems that do not understand ld -x),
## omit -x in the ld line above.
## If your system does not have the ld command, comment out
## or remove both the ld and mv lines above.

MISC =	f77vers.o i77vers.o main.o s_rnge.o abort_.o exit_.o getarg_.o iargc_.o\
	getenv_.o signal_.o s_stop.o s_paus.o system_.o cabs.o ctype.o\
	derf_.o derfc_.o erf_.o erfc_.o sig_die.o uninit.o
POW =	pow_ci.o pow_dd.o pow_di.o pow_hh.o pow_ii.o pow_ri.o pow_zi.o pow_zz.o
CX =	c_abs.o c_cos.o c_div.o c_exp.o c_log.o c_sin.o c_sqrt.o
DCX =	z_abs.o z_cos.o z_div.o z_exp.o z_log.o z_sin.o z_sqrt.o
REAL =	r_abs.o r_acos.o r_asin.o r_atan.o r_atn2.o r_cnjg.o r_cos.o\
	r_cosh.o r_dim.o r_exp.o r_imag.o r_int.o\
	r_lg10.o r_log.o r_mod.o r_nint.o r_sign.o\
	r_sin.o r_sinh.o r_sqrt.o r_tan.o r_tanh.o
DBL =	d_abs.o d_acos.o d_asin.o d_atan.o d_atn2.o\
	d_cnjg.o d_cos.o d_cosh.o d_dim.o d_exp.o\
	d_imag.o d_int.o d_lg10.o d_log.o d_mod.o\
	d_nint.o d_prod.o d_sign.o d_sin.o d_sinh.o\
	d_sqrt.o d_tan.o d_tanh.o
INT =	i_abs.o i_dim.o i_dnnt.o i_indx.o i_len.o i_mod.o i_nint.o i_sign.o\
	lbitbits.o lbitshft.o
HALF =	h_abs.o h_dim.o h_dnnt.o h_indx.o h_len.o h_mod.o h_nint.o h_sign.o
CMP =	l_ge.o l_gt.o l_le.o l_lt.o hl_ge.o hl_gt.o hl_le.o hl_lt.o
EFL =	ef1asc_.o ef1cmc_.o
CHAR =	f77_aloc.o s_cat.o s_cmp.o s_copy.o
I77 =	backspac.o close.o dfe.o dolio.o due.o endfile.o err.o\
	fmt.o fmtlib.o ftell_.o iio.o ilnw.o inquire.o lread.o lwrite.o\
	open.o rdfmt.o rewind.o rsfe.o rsli.o rsne.o sfe.o sue.o\
	typesize.o uio.o util.o wref.o wrtfmt.o wsfe.o wsle.o wsne.o xwsne.o
QINT =	pow_qq.o qbitbits.o qbitshft.o ftell64_.o
TIME =	dtime_.o etime_.o

# If you get an error compiling dtime_.c or etime_.c, try adding
# -DUSE_CLOCK to the CFLAGS assignment above; if that does not work,
# omit $(TIME) from OFILES = assignment below.

# To get signed zeros in write statements on IEEE-arithmetic systems,
# add -DSIGNED_ZEROS to the CFLAGS assignment below and add signbit.o
# to the end of the OFILES = assignment below.

# For INTEGER*8 support (which requires system-dependent adjustments to
# f2c.h), add $(QINT) to the OFILES = assignment below...

OFILES = $(MISC) $(POW) $(CX) $(DCX) $(REAL) $(DBL) $(INT) \
	$(HALF) $(CMP) $(EFL) $(CHAR) $(I77) $(TIME)

all: f2c.h signal1.h sysdep1.h $(OFILES)

libf2c.a: $(OFILES)
	ar r libf2c.a $?
	-ranlib libf2c.a

## Shared-library variant: the following rule works on Linux
## systems.  Details are system-dependent.  Under Linux, -fPIC
## must appear in the CFLAGS assignment when making libf2c.so.
## Under Solaris, use -Kpic in CFLAGS and use "ld -G" instead
## of "$(CC) -shared".
## For MacOSX 10.4 and 10.5 (and perhaps other versions >= 10.3), use
## "MACOSX_DEPLOYMENT_TARGET=10.3 libtool -dynamic -undefined dynamic_lookup -single_module"
## instead of "$(CC) -shared", and when running programs linked against libf2c.so,
## arrange for $DYLD_LIBRARY_PATH to include the directory containing libf2c.so.

libf2c.so: $(OFILES)
	$(CC) -shared -o libf2c.so $(OFILES)

### If your system lacks ranlib, you don't need it; see README.

f77vers.o: f77vers.c
	$(CC) -c f77vers.c

i77vers.o: i77vers.c
	$(CC) -c i77vers.c

# To get an "f2c.h" for use with "f2c -C++", first "make hadd"
hadd: f2c.h0 f2ch.add
	cat f2c.h0 f2ch.add >f2c.h

# For use with "f2c" and "f2c -A":
f2c.h: f2c.h0
	cp f2c.h0 f2c.h

# You may need to adjust signal1.h and sysdep1.h suitably for your system...
signal1.h: signal1.h0
	cp signal1.h0 signal1.h

sysdep1.h: sysdep1.h0
	cp sysdep1.h0 sysdep1.h

# If your system lacks onexit() and you are not using an
# ANSI C compiler, then you should uncomment the following
# two lines (for compiling main.o):
#main.o: main.c
#	$(CC) -c -DNO_ONEXIT -DSkip_f2c_Undefs main.c
# On at least some Sun systems, it is more appropriate to
# uncomment the following two lines:
#main.o: main.c
#	$(CC) -c -Donexit=on_exit -DSkip_f2c_Undefs main.c

install: libf2c.a
	cp libf2c.a $(LIBDIR)
	-ranlib $(LIBDIR)/libf2c.a

clean:
	rm -f *.o arith.h signal1.h sysdep1.h

backspac.o:	fio.h
close.o:	fio.h
dfe.o:		fio.h
dfe.o:		fmt.h
due.o:		fio.h
endfile.o:	fio.h rawio.h
err.o:		fio.h rawio.h
fmt.o:		fio.h
fmt.o:		fmt.h
iio.o:		fio.h
iio.o:		fmt.h
ilnw.o:		fio.h
ilnw.o:		lio.h
inquire.o:	fio.h
lread.o:	fio.h
lread.o:	fmt.h
lread.o:	lio.h
lread.o:	fp.h
lwrite.o:	fio.h
lwrite.o:	fmt.h
lwrite.o:	lio.h
open.o:		fio.h rawio.h
rdfmt.o:	fio.h
rdfmt.o:	fmt.h
rdfmt.o:	fp.h
rewind.o:	fio.h
rsfe.o:		fio.h
rsfe.o:		fmt.h
rsli.o:		fio.h
rsli.o:		lio.h
rsne.o:		fio.h
rsne.o:		lio.h
sfe.o:		fio.h
signbit.o:	arith.h
sue.o:		fio.h
uio.o:		fio.h
uninit.o:	arith.h
util.o:		fio.h
wref.o:		fio.h
wref.o:		fmt.h
wref.o:		fp.h
wrtfmt.o:	fio.h
wrtfmt.o:	fmt.h
wsfe.o:		fio.h
wsfe.o:		fmt.h
wsle.o:		fio.h
wsle.o:		fmt.h
wsle.o:		lio.h
wsne.o:		fio.h
wsne.o:		lio.h
xwsne.o:	fio.h
xwsne.o:	lio.h
xwsne.o:	fmt.h

arith.h: arithchk.c
	$(CC) $(CFLAGS) -DNO_FPINIT arithchk.c -lm ||\
	 $(CC) -DNO_LONG_LONG $(CFLAGS) -DNO_FPINIT arithchk.c -lm
	./a.out >arith.h
	rm -f a.out arithchk.o

check:
	xsum Notice README abort_.c arithchk.c backspac.c c_abs.c c_cos.c \
	c_div.c c_exp.c c_log.c c_sin.c c_sqrt.c cabs.c close.c comptry.bat \
	ctype.c ctype.h \
	d_abs.c d_acos.c d_asin.c d_atan.c d_atn2.c d_cnjg.c d_cos.c d_cosh.c \
	d_dim.c d_exp.c d_imag.c d_int.c d_lg10.c d_log.c d_mod.c \
	d_nint.c d_prod.c d_sign.c d_sin.c d_sinh.c d_sqrt.c d_tan.c \
	d_tanh.c derf_.c derfc_.c dfe.c dolio.c dtime_.c due.c ef1asc_.c \
	ef1cmc_.c endfile.c erf_.c erfc_.c err.c etime_.c exit_.c f2c.h0 \
	f2ch.add f77_aloc.c f77vers.c fio.h fmt.c fmt.h fmtlib.c \
	fp.h ftell_.c ftell64_.c \
	getarg_.c getenv_.c h_abs.c h_dim.c h_dnnt.c h_indx.c h_len.c \
	h_mod.c h_nint.c h_sign.c hl_ge.c hl_gt.c hl_le.c hl_lt.c \
	i77vers.c i_abs.c i_dim.c i_dnnt.c i_indx.c i_len.c i_mod.c \
	i_nint.c i_sign.c iargc_.c iio.c ilnw.c inquire.c l_ge.c l_gt.c \
	l_le.c l_lt.c lbitbits.c lbitshft.c libf2c.lbc libf2c.sy lio.h \
	lread.c lwrite.c main.c makefile.sy makefile.u makefile.vc \
	makefile.wat math.hvc mkfile.plan9 open.c pow_ci.c pow_dd.c \
	pow_di.c pow_hh.c pow_ii.c pow_qq.c pow_ri.c pow_zi.c pow_zz.c \
	qbitbits.c qbitshft.c r_abs.c r_acos.c r_asin.c r_atan.c r_atn2.c \
	r_cnjg.c r_cos.c r_cosh.c r_dim.c r_exp.c r_imag.c r_int.c r_lg10.c \
	r_log.c r_mod.c r_nint.c r_sign.c r_sin.c r_sinh.c r_sqrt.c \
	r_tan.c r_tanh.c rawio.h rdfmt.c rewind.c rsfe.c rsli.c rsne.c \
	s_cat.c s_cmp.c s_copy.c s_paus.c s_rnge.c s_stop.c scomptry.bat sfe.c \
	sig_die.c signal1.h0 signal_.c signbit.c sue.c sysdep1.h0 system_.c \
	typesize.c \
	uio.c uninit.c util.c wref.c wrtfmt.c wsfe.c wsle.c wsne.c xwsne.c \
	z_abs.c z_cos.c z_div.c z_exp.c z_log.c z_sin.c z_sqrt.c >xsum1.out
	cmp xsum0.out xsum1.out && mv xsum1.out xsum.out || diff xsum[01].out
