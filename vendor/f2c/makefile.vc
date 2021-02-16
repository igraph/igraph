# For making f2c.lib (here called vcf2c.lib) with Microsoft Visual C++ .
# Invoke with "nmake -f makefile.vc" .

# To get signed zeros in write statements on IEEE-arithmetic systems,
# add -DSIGNED_ZEROS to the CFLAGS assignment below and add signbit.obj
# to the objects in the "w =" list below.

CC = cl
CFLAGS = -DUSE_CLOCK -DMSDOS -DNO_ONEXIT -Ot1 -DNO_My_ctype -DNO_ISATTY

.c.obj:
	$(CC) -c $(CFLAGS) $*.c

w = \
	abort_.obj \
	backspac.obj \
	c_abs.obj \
	c_cos.obj \
	c_div.obj \
	c_exp.obj \
	c_log.obj \
	c_sin.obj \
	c_sqrt.obj \
	cabs.obj \
	close.obj \
	d_abs.obj \
	d_acos.obj \
	d_asin.obj \
	d_atan.obj \
	d_atn2.obj \
	d_cnjg.obj \
	d_cos.obj \
	d_cosh.obj \
	d_dim.obj \
	d_exp.obj \
	d_imag.obj \
	d_int.obj \
	d_lg10.obj \
	d_log.obj \
	d_mod.obj \
	d_nint.obj \
	d_prod.obj \
	d_sign.obj \
	d_sin.obj \
	d_sinh.obj \
	d_sqrt.obj \
	d_tan.obj \
	d_tanh.obj \
	derf_.obj \
	derfc_.obj \
	dfe.obj \
	dolio.obj \
	dtime_.obj \
	due.obj \
	ef1asc_.obj \
	ef1cmc_.obj \
	endfile.obj \
	erf_.obj \
	erfc_.obj \
	err.obj \
	etime_.obj \
	exit_.obj \
	f77_aloc.obj \
	f77vers.obj \
	fmt.obj \
	fmtlib.obj \
	ftell_.obj \
	getarg_.obj \
	getenv_.obj \
	h_abs.obj \
	h_dim.obj \
	h_dnnt.obj \
	h_indx.obj \
	h_len.obj \
	h_mod.obj \
	h_nint.obj \
	h_sign.obj \
	hl_ge.obj \
	hl_gt.obj \
	hl_le.obj \
	hl_lt.obj \
	i77vers.obj \
	i_abs.obj \
	i_dim.obj \
	i_dnnt.obj \
	i_indx.obj \
	i_len.obj \
	i_mod.obj \
	i_nint.obj \
	i_sign.obj \
	iargc_.obj \
	iio.obj \
	ilnw.obj \
	inquire.obj \
	l_ge.obj \
	l_gt.obj \
	l_le.obj \
	l_lt.obj \
	lbitbits.obj \
	lbitshft.obj \
	lread.obj \
	lwrite.obj \
	main.obj \
	open.obj \
	pow_ci.obj \
	pow_dd.obj \
	pow_di.obj \
	pow_hh.obj \
	pow_ii.obj \
	pow_ri.obj \
	pow_zi.obj \
	pow_zz.obj \
	r_abs.obj \
	r_acos.obj \
	r_asin.obj \
	r_atan.obj \
	r_atn2.obj \
	r_cnjg.obj \
	r_cos.obj \
	r_cosh.obj \
	r_dim.obj \
	r_exp.obj \
	r_imag.obj \
	r_int.obj \
	r_lg10.obj \
	r_log.obj \
	r_mod.obj \
	r_nint.obj \
	r_sign.obj \
	r_sin.obj \
	r_sinh.obj \
	r_sqrt.obj \
	r_tan.obj \
	r_tanh.obj \
	rdfmt.obj \
	rewind.obj \
	rsfe.obj \
	rsli.obj \
	rsne.obj \
	s_cat.obj \
	s_cmp.obj \
	s_copy.obj \
	s_paus.obj \
	s_rnge.obj \
	s_stop.obj \
	sfe.obj \
	sig_die.obj \
	signal_.obj \
	sue.obj \
	system_.obj \
	typesize.obj \
	uio.obj \
	uninit.obj \
	util.obj \
	wref.obj \
	wrtfmt.obj \
	wsfe.obj \
	wsle.obj \
	wsne.obj \
	xwsne.obj \
	z_abs.obj \
	z_cos.obj \
	z_div.obj \
	z_exp.obj \
	z_log.obj \
	z_sin.obj \
	z_sqrt.obj

all: f2c.h math.h signal1.h sysdep1.h vcf2c.lib

f2c.h: f2c.h0
	copy f2c.h0 f2c.h

math.h: math.hvc
	copy math.hvc math.h

signal1.h: signal1.h0
	copy signal1.h0 signal1.h

sysdep1.h: sysdep1.h0
	copy sysdep1.h0 sysdep1.h

vcf2c.lib: $w
	lib -out:vcf2c.lib @libf2c.lbc

open.obj: open.c
	$(CC) -c $(CFLAGS) -DMSDOS open.c

signbit.obj uninit.obj: arith.h

arith.h: arithchk.c
	comptry.bat $(CC) $(CFLAGS) -DNO_FPINIT arithchk.c
	arithchk >arith.h
	del arithchk.exe
	del arithchk.obj
