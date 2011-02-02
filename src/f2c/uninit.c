#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "arith.h"

#define TYSHORT 2
#define TYLONG 3
#define TYREAL 4
#define TYDREAL 5
#define TYCOMPLEX 6
#define TYDCOMPLEX 7
#define TYINT1 11
#define TYQUAD 14
#ifndef Long
#define Long long
#endif

#ifdef __mips
#define RNAN	0xffc00000 /* Quiet NaN */
#define DNAN0	0xfff80000 /* Signalling NaN double Big endian */
#define DNAN1	0
#endif

#ifdef _PA_RISC1_1
#define RNAN	0xffc00000 /* Quiet Nan -- big endian */
#define DNAN0	0xfff80000
#define DNAN1	0
#endif

#ifndef RNAN
#define RNAN	0xff800001
#ifdef IEEE_MC68k /* set on PPC*/
#define DNAN0	0xfff00000 /* Quiet NaN big endian */
#define DNAN1	1
#else
#define DNAN0	1   /* LSB, MSB for little endian machines */
#define DNAN1	0xfff00000
#endif
#endif /*RNAN*/

#ifdef KR_headers
#define Void /*void*/
#define FA7UL (unsigned Long) 0xfa7a7a7aL
#else
#define Void void
#define FA7UL 0xfa7a7a7aUL
#endif

#ifdef __cplusplus
extern "C" {
#endif

static void ieee0(Void);

static unsigned Long rnan = RNAN,
	dnan0 = DNAN0,
	dnan1 = DNAN1;

double _0 = 0.;

void unsupported_error()
{
  fprintf(stderr,"Runtime Error: Your Architecture is not supported by the"
                       " -trapuv option of f2c\n");
  exit(-1);
}



 void
#ifdef KR_headers
_uninit_f2c(x, type, len) void *x; int type; long len;
#else
_uninit_f2c(void *x, int type, long len)
#endif
{
	static int first = 1;

	unsigned Long *lx, *lxe;

	if (first) {
		first = 0;
		ieee0();
		}
	if (len == 1)
	 switch(type) {
	  case TYINT1:
		*(char*)x = 'Z';
		return;
	  case TYSHORT:
		*(short*)x = 0xfa7a;
		break;
	  case TYLONG:
		*(unsigned Long*)x = FA7UL;
		return;
	  case TYQUAD:
	  case TYCOMPLEX:
	  case TYDCOMPLEX:
		break;
	  case TYREAL:
		*(unsigned Long*)x = rnan;
		return;
	  case TYDREAL:
		lx = (unsigned Long*)x;
		lx[0] = dnan0;
		lx[1] = dnan1;
		return;
	  default:
		printf("Surprise type %d in _uninit_f2c\n", type);
	  }
	switch(type) {
	  case TYINT1:
		memset(x, 'Z', len);
		break;
	  case TYSHORT:
		*(short*)x = 0xfa7a;
		break;
	  case TYQUAD:
		len *= 2;
		/* no break */
	  case TYLONG:
		lx = (unsigned Long*)x;
		lxe = lx + len;
		while(lx < lxe)
			*lx++ = FA7UL;
		break;
	  case TYCOMPLEX:
		len *= 2;
		/* no break */
	  case TYREAL:
		lx = (unsigned Long*)x;
		lxe = lx + len;
		while(lx < lxe)
			*lx++ = rnan;
		break;
	  case TYDCOMPLEX:
		len *= 2;
		/* no break */
	  case TYDREAL:
		lx = (unsigned Long*)x;
		for(lxe = lx + 2*len; lx < lxe; lx += 2) {
			lx[0] = dnan0;
			lx[1] = dnan1;
			}
	  }
	}
#ifdef __cplusplus
}
#endif

#ifndef MSpc
#ifdef MSDOS
#define MSpc
#else
#ifdef _WIN32
#define MSpc
#endif
#endif
#endif

#ifdef MSpc
#define IEEE0_done
#include "float.h"
#include "signal.h"

 static void
ieee0(Void)
{
#ifndef __alpha
#ifndef EM_DENORMAL
#define EM_DENORMAL _EM_DENORMAL
#endif
#ifndef EM_UNDERFLOW
#define EM_UNDERFLOW _EM_UNDERFLOW
#endif
#ifndef EM_INEXACT
#define EM_INEXACT _EM_INEXACT
#endif
#ifndef MCW_EM
#define MCW_EM _MCW_EM
#endif
	_control87(EM_DENORMAL | EM_UNDERFLOW | EM_INEXACT, MCW_EM);
#endif
	/* With MS VC++, compiling and linking with -Zi will permit */
	/* clicking to invoke the MS C++ debugger, which will show */
	/* the point of error -- provided SIGFPE is SIG_DFL. */
	signal(SIGFPE, SIG_DFL);
	}
#endif /* MSpc */

/* What follows is for SGI IRIX only */
#if defined(__mips) && defined(__sgi)   /* must link with -lfpe */
#define IEEE0_done
/* code from Eric Grosse */
#include <stdlib.h>
#include <stdio.h>
#include "/usr/include/sigfpe.h"	/* full pathname for lcc -N */
#include "/usr/include/sys/fpu.h"

 static void
#ifdef KR_headers
ieeeuserhand(exception, val) unsigned exception[5]; int val[2];
#else
ieeeuserhand(unsigned exception[5], int val[2])
#endif
{
	fflush(stdout);
	fprintf(stderr,"ieee0() aborting because of ");
	if(exception[0]==_OVERFL) fprintf(stderr,"overflow\n");
	else if(exception[0]==_UNDERFL) fprintf(stderr,"underflow\n");
	else if(exception[0]==_DIVZERO) fprintf(stderr,"divide by 0\n");
	else if(exception[0]==_INVALID) fprintf(stderr,"invalid operation\n");
	else fprintf(stderr,"\tunknown reason\n");
	fflush(stderr);
	abort();
}

 static void
#ifdef KR_headers
ieeeuserhand2(j) unsigned int **j;
#else
ieeeuserhand2(unsigned int **j)
#endif
{
	fprintf(stderr,"ieee0() aborting because of confusion\n");
	abort();
}

 static void
ieee0(Void)
{
	int i;
	for(i=1; i<=4; i++){
		sigfpe_[i].count = 1000;
		sigfpe_[i].trace = 1;
		sigfpe_[i].repls = _USER_DETERMINED;
		}
	sigfpe_[1].repls = _ZERO;	/* underflow */
	handle_sigfpes( _ON,
		_EN_UNDERFL|_EN_OVERFL|_EN_DIVZERO|_EN_INVALID,
		ieeeuserhand,_ABORT_ON_ERROR,ieeeuserhand2);
	}
#endif /* IRIX mips */

/*
 * The following is the preferred method but depends upon a GLIBC extension only
 * to be found in GLIBC 2.2 or later.  It is a GNU extension, not included in the
 * C99 extensions which allow the FP status register to be examined in a platform
 * independent way.  It should be used if at all possible  -- AFRB
 */


#if (defined(__GLIBC__)&& ( __GLIBC__>=2) && (__GLIBC_MINOR__>=2) )
#define _GNU_SOURCE 1
#define IEEE0_done
#include <fenv.h>
 static void
  ieee0(Void)
        
{
    /* Clear all exception flags */
    if (fedisableexcept(FE_ALL_EXCEPT)==-1)
         unsupported_error();
    if (feenableexcept(FE_DIVBYZERO|FE_INVALID|FE_OVERFLOW)==-1)
         unsupported_error();
}

#endif /* Glibc control */

/* Many linux cases will be treated through GLIBC.  Note that modern
 * linux runs on many non-i86 plaforms and as a result the following code
 * must be processor dependent rather than simply OS specific */

#if (defined(__linux__)&&(!defined(IEEE0_done)))
#define IEEE0_done
#include <fpu_control.h>


#ifdef __alpha__
#ifndef USE_setfpucw
#define __setfpucw(x) __fpu_control = (x)
#endif
#endif

/* Not all versions of libc define _FPU_SETCW;
 *  * some only provide the __setfpucw() function.
 *   */
#ifndef _FPU_SETCW
#define _FPU_SETCW(cw) __setfpucw(cw)
#endif

/* The exact set of flags we want to set in the FPU control word
 * depends on the architecture.
 * Note also that whether an exception is enabled or disabled when
 * the _FPU_MASK_nn bit is set is architecture dependent!
 * Enabled-when-set: M68k, ARM, MIPS, PowerPC
 * Disabled-when-set: x86, Alpha
 * The state we are after is:
 * exceptions on division by zero, overflow and invalid operation.
 */


#ifdef __alpha__
#ifndef USE_setfpucw
#define __setfpucw(x) __fpu_control = (x)
#endif
#endif


#ifndef _FPU_SETCW
#undef  Can_use__setfpucw
#define Can_use__setfpucw
#endif

#undef RQD_FPU_MASK
#undef RQD_FPU_CLEAR_MASK

#if (defined(__mc68000__) || defined(__mc68020__) || defined(mc68020) || defined (__mc68k__))
/* Reported 20010705 by Alan Bain <alanb@chiark.greenend.org.uk> */
/* Note that IEEE 754 IOP (illegal operation) */
/* = Signaling NAN (SNAN) + operation error (OPERR). */
#define RQD_FPU_STATE (_FPU_IEEE + _FPU_DOUBLE + _FPU_MASK_OPERR + \
                 _FPU_MASK_DZ + _FPU_MASK_SNAN+_FPU_MASK_OVFL)
#define RQD_FPU_MASK (_FPU_MASK_OPERR+_FPU_MASK_DZ+_FPU_MASK_SNAN+_FPU_MASK_OVFL)

#elif (defined(__powerpc__)||defined(_ARCH_PPC)||defined(_ARCH_PWR)) /* !__mc68k__ */
    /* The following is NOT a mistake -- the author of the fpu_control.h
     * for the PPC has erroneously defined IEEE mode to turn on exceptions
     * other than Inexact! Start from default then and turn on only the ones
     * which we want*/

    /* I have changed _FPU_MASK_UM here to _FPU_MASK_ZM, because that is
     * in line with all the other architectures specified here. -- AFRB
     */
#define RQD_FPU_STATE (_FPU_DEFAULT +_FPU_MASK_OM+_FPU_MASK_IM+_FPU_MASK_ZM)
#define RQD_FPU_MASK (_FPU_MASK_OM+_FPU_MASK_IM+_FPU_MASK_ZM)

#elif (defined(__arm__))
    /* On ARM too, IEEE implies all exceptions enabled.
     * -- Peter Maydell <pmaydell@chiark.greenend.org.uk>
     * Unfortunately some version of ARMlinux don't include any
     * flags in the fpu_control.h file
     */
#define RQD_FPU_STATE (_FPU_DEFAULT +_FPU_MASK_OM+_FPU_MASK_IM+_FPU_MASK_ZM)
#define RQD_FPU_MASK (_FPU_MASK_OM+_FPU_MASK_IM+_FPU_MASK_ZM)

#elif (defined(__mips__))
    /* And same again for MIPS; _FPU_IEEE => exceptions seems a common meme.
     *  * MIPS uses different MASK constant names, no idea why -- PMM
     *   */
#define RQD_FPU_STATE (_FPU_DEFAULT +_FPU_MASK_O+_FPU_MASK_V+_FPU_MASK_Z)
#define RQD_FPU_MASK (_FPU_MASK_O+_FPU_MASK_V+_FPU_MASK_Z)

#elif (defined(__sparc__))
#define RQD_FPU_STATE (_FPU_DEFAULT +_FPU_DOUBLE+_FPU_MASK_OM+_FPU_MASK_IM+_FPU_MASK_ZM)
#define RQD_FPU_MASK (_FPU_MASK_OM+_FPU_MASK_IM+_FPU_MASK_ZM)

#elif (defined(__i386__) || defined(__alpha__))
    /* This case is for Intel, and also Alpha, because the Alpha header 
     * purposely emulates x86 flags and meanings for compatibility with
     * stupid programs.
     * We used to try this case for anything defining _FPU_IEEE, but I think
     * that that's a bad idea because it isn't really likely to work.
     * Instead for unknown architectures we just won't allow -trapuv to work.
     * Trying this case was just getting us 
     *  (a) compile errors on archs which didn't know all these constants
     *  (b) silent wrong behaviour on archs (like SPARC) which do know all
     *      constants but have different semantics for them
     */
#define RQD_FPU_STATE (_FPU_IEEE - _FPU_EXTENDED + _FPU_DOUBLE - _FPU_MASK_IM - _FPU_MASK_ZM - _FPU_MASK_OM)
#define RQD_FPU_CLEAR_MASK (_FPU_MASK_IM + _FPU_MASK_ZM + _FPU_MASK_OM)
#endif

static void ieee0(Void)
{
#ifdef RQD_FPU_STATE
        
#ifndef UNINIT_F2C_PRECISION_53 /* 20051004 */
        __fpu_control = RQD_FPU_STATE;
        _FPU_SETCW(__fpu_control);
#else 
	/* unmask invalid, etc., and keep current rounding precision */
	fpu_control_t cw;
	_FPU_GETCW(cw);
#ifdef RQD_FPU_CLEAR_MASK
	cw &= ~ RQD_FPU_CLEAR_MASK;
#else
        cw |= RQD_FPU_MASK;
#endif
	_FPU_SETCW(cw);
#endif

#else /* !_FPU_IEEE */

	fprintf(stderr, "\n%s\n%s\n%s\n%s\n",
		"WARNING:  _uninit_f2c in libf2c does not know how",
		"to enable trapping on this system, so f2c's -trapuv",
		"option will not detect uninitialized variables unless",
		"you can enable trapping manually.");
	fflush(stderr);

#endif /* _FPU_IEEE */
	}
#endif /* __linux__ */

/* Specific to OSF/1 */
#if (defined(__alpha)&&defined(__osf__))
#ifndef IEEE0_done
#define IEEE0_done
#include <machine/fpu.h>
 static void
ieee0(Void)
{
	ieee_set_fp_control(IEEE_TRAP_ENABLE_INV);
	}
#endif /*IEEE0_done*/
#endif /*__alpha OSF/1*/

#ifdef __hpux
#define IEEE0_done
#define _INCLUDE_HPUX_SOURCE
#include <math.h>

#ifndef FP_X_INV
#include <fenv.h>
#define fpsetmask fesettrapenable
#define FP_X_INV FE_INVALID
#endif

 static void
ieee0(Void)
{
	fpsetmask(FP_X_INV);
	}
#endif /*__hpux*/

#ifdef _AIX
#define IEEE0_done
#include <fptrap.h>

 static void
ieee0(Void)
{
	fp_enable(TRP_INVALID);
	fp_trap(FP_TRAP_SYNC);
	}
#endif /*_AIX*/

#ifdef __sun
#define IEEE0_done
#include <ieeefp.h>

 static void
ieee0(Void)
{
	fpsetmask(FP_X_INV);
	}
#endif /*__sparc*/

#ifndef IEEE0_done
 static void
ieee0(Void) {}
#endif
