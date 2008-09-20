/* -*- mode: C -*-  */
/* 
   IGraph library.
   Copyright (C) 2005  Gabor Csardi <csardi@rmki.kfki.hu>
   MTA RMKI, Konkoly-Thege Miklos st. 29-33, Budapest 1121, Hungary
   
   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
   
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
   
   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 
   02110-1301 USA

*/

#include "random.h"
#include "error.h"
#include "config.h"

#include <math.h>
#include "igraph_math.h"

int igraph_rng_inited = 0;

#ifndef HAVE_EXPM1
#ifndef USING_R			/* R provides a replacement */
/* expm1 replacement */
static double expm1 (double x)
{
    if (fabs(x) < M_LN2)
    {
        /* Compute the taylor series S = x + (1/2!) x^2 + (1/3!) x^3 + ... */

        double i = 1.0;
        double sum = x;
        double term = x / 1.0;

        do
        {
            term *= x / ++i;
            sum += term;
        }
        while (fabs(term) > fabs(sum) * 2.22e-16);
      
        return sum;
    }

    return expl(x) - 1.0L;
}
#endif
#endif

#ifndef HAVE_RINT
#ifndef USING_R			/* R provides a replacement */
/* rint replacement */
static double rint (double x)
{
   return ( (x<0.) ? -floor(-x+.5) : floor(x+.5) );
}
#endif
#endif

#if HAVE_RINTF != 1
static float rintf (float x)
{
   return ( (x<(float)0.) ? -(float)floor(-x+.5) : (float)floor(x+.5) );
}
#endif

/*
 * \ingroup internal
 * 
 * This function appends the rest of the needed random number to the 
 * result vector.
 */

int igraph_random_sample_alga(igraph_vector_t *res, igraph_integer_t l, igraph_integer_t h, 
			      igraph_integer_t length) {
  igraph_real_t N=h-l+1;
  igraph_real_t n=length;
  
  igraph_real_t top=N-n;
  igraph_real_t Nreal=N;
  igraph_real_t S=0;
  igraph_real_t V, quot;
  
  l=l-1;

  while (n>=2) {
    V=RNG_UNIF01();
    S=1;
    quot=top/Nreal;
    while (quot>V) {
      S+=1;
      top=-1.0+top;
      Nreal=-1.0+Nreal;
      quot=(quot*top)/Nreal;
    }
    l+=S;
    igraph_vector_push_back(res, l);	/* allocated */
    Nreal=-1.0+Nreal; n=-1+n;
  }
  
  S=floor(round(Nreal)*RNG_UNIF01());
  l+=S+1;
  igraph_vector_push_back(res, l);	/* allocated */
  
  return 0;
}

/**
 * \ingroup nongraph
 * \function igraph_random_sample
 * \brief Generates an increasing random sequence of integers.
 * 
 * </para><para>
 * This function generates an incresing sequence of random integer
 * numbers from a given interval. The algorithm is taken literally
 * from Jeffrey Scott Vitter: 'An Efficient Algorithm for Sequential
 * Random Sampling', ACM Transactions on Mathematical Software, 13/1,
 * 58--67. This method can be used for generating numbers from a
 * \em very large interval, it is primilarly created for randomly
 * selecting some edges from the sometimes huge set of possible edges
 * in a large graph.
 * \param res Pointer to an initialized vector, this will hold the
 *        result. It will be resized to the proper size.
 * \param l The lower limit of the generation interval (inclusive).
 * \param h The upper limit of the generation interval (inclusive).
 * \param length The number of random integers to generate.
 * \return Error code.
 *
 * Time complexity: according to the referenced paper, the expected
 * running time is O(length).
 */

int igraph_random_sample(igraph_vector_t *res, igraph_integer_t l, igraph_integer_t h, 
			 igraph_integer_t length) {
  igraph_real_t N=h-l+1;
  igraph_real_t n=length;
  int retval;

  igraph_real_t nreal=length;
  igraph_real_t ninv=1.0/nreal;
  igraph_real_t Nreal=N;
  igraph_real_t Vprime;
  igraph_real_t qu1=-n+1+N;
  igraph_real_t qu1real=-nreal+1.0+Nreal;
  igraph_real_t negalphainv=-13;
  igraph_real_t threshold=-negalphainv*n;
  igraph_real_t S;
  
  igraph_vector_clear(res);
  IGRAPH_CHECK(igraph_vector_reserve(res, length));  

  RNG_BEGIN();
  
  Vprime=exp(log(RNG_UNIF01())*ninv);

  while (n>1 && threshold < N) {
    igraph_real_t X, U;
    igraph_real_t limit, t;
    igraph_real_t negSreal, y1, y2, top, bottom;
    igraph_real_t nmin1inv=1.0/(-1.0+nreal);
    while (1) {
      while(1) {
	X=Nreal*(-Vprime+1.0);
	S=floor(X);
	if (S==0) { S=1; }
	if (S <qu1) { break; }
	Vprime = exp(log(RNG_UNIF01())*ninv);
      }
      U=RNG_UNIF01();
      negSreal=-S;
      
      y1=exp(log(U*Nreal/qu1real)*nmin1inv);
      Vprime=y1*(-X/Nreal+1.0)*(qu1real/(negSreal+qu1real));
      if (Vprime <= 1.0) { break; }
      
      y2=1.0;
      top=-1.0+Nreal;
      if (-1+n > S) {
	bottom=-nreal+Nreal; 
	limit=-S+N;
      } else {
	bottom=-1.0+negSreal+Nreal;
	limit=qu1;
      }
      for (t=-1+N; t>=limit; t--) {
	y2=(y2*top)/bottom;
	top=-1.0+top;
	bottom=-1.0+bottom;
      }
      if (Nreal/(-X+Nreal) >= y1*exp(log(y2)*nmin1inv)) {
	Vprime=exp(log(RNG_UNIF01())*nmin1inv);
	break;
      }
      Vprime=exp(log(RNG_UNIF01())*ninv);
    }
        
    l+=S;
    igraph_vector_push_back(res, l);	/* allocated */
    N=-S+(-1+N);   Nreal=negSreal+(-1.0+Nreal);
    n=-1+n;   nreal=-1.0+nreal; ninv=nmin1inv;
    qu1=-S+qu1; qu1real=negSreal+qu1real;
    threshold=threshold+negalphainv;
  }
  
  if (n>1) {
    retval=igraph_random_sample_alga(res, l, h, n);
  } else {
    retval=0;
    S=floor(N*Vprime);
    l+=S;
    igraph_vector_push_back(res, l);	/* allocated */
  }

  RNG_END();
  
  return retval;
}
  
#ifdef USING_R

#else

/*
 *  Mathlib : A C Library of Special Functions
 *  Copyright (C) 1998 Ross Ihaka
 *  Copyright (C) 2000 The R Development Core Team
 *  based on AS 111 (C) 1977 Royal Statistical Society
 *  and   on AS 241 (C) 1988 Royal Statistical Society
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA.
 *
 *  SYNOPSIS
 *
 *	double qnorm5(double p, double mu, double sigma,
 *		      int lower_tail, int log_p)
 *            {qnorm (..) is synonymous and preferred inside R}
 *
 *  DESCRIPTION
 *
 *	Compute the quantile function for the normal distribution.
 *
 *	For small to moderate probabilities, algorithm referenced
 *	below is used to obtain an initial approximation which is
 *	polished with a final Newton step.
 *
 *	For very large arguments, an algorithm of Wichura is used.
 *
 *  REFERENCE
 *
 *	Beasley, J. D. and S. G. Springer (1977).
 *	Algorithm AS 111: The percentage points of the normal distribution,
 *	Applied Statistics, 26, 118-121.
 *
 *      Wichura, M.J. (1988).
 *      Algorithm AS 241: The Percentage Points of the Normal Distribution.
 *      Applied Statistics, 37, 477-484.
 */

/*
 *  Mathlib : A C Library of Special Functions
 *  Copyright (C) 1998-2004  The R Development Core Team
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

/* Private header file for use during compilation of Mathlib */
#ifndef MATHLIB_PRIVATE_H
#define MATHLIB_PRIVATE_H

#ifdef _MSC_VER
#  define ML_POSINF IGRAPH_INFINITY
#  define ML_NEGINF -IGRAPH_INFINITY
#  define ML_NAN    IGRAPH_NAN
#else
#  define ML_POSINF	(1.0 / 0.0)
#  define ML_NEGINF	((-1.0) / 0.0)
#  define ML_NAN		(0.0 / 0.0)
#endif

#define ML_ERROR(x)	/* nothing */
#define ML_UNDERFLOW	(DBL_MIN * DBL_MIN)
#define ML_VALID(x)	(!ISNAN(x))

#define ME_NONE		0
/*	no error */
#define ME_DOMAIN	1
/*	argument out of domain */
#define ME_RANGE	2
/*	value out of range */
#define ME_NOCONV	4
/*	process did not converge */
#define ME_PRECISION	8
/*	does not have "full" precision */
#define ME_UNDERFLOW	16
/*	and underflow occured (important for IEEE)*/

#define ML_ERR_return_NAN { ML_ERROR(ME_DOMAIN); return ML_NAN; }

/* Wilcoxon Rank Sum Distribution */

#define WILCOX_MAX 50

/* Wilcoxon Signed Rank Distribution */

#define SIGNRANK_MAX 50

/* Formerly private part of Mathlib.h */

/* always remap internal functions */
#define bd0       	Rf_bd0
#define chebyshev_eval	Rf_chebyshev_eval
#define chebyshev_init	Rf_chebyshev_init
#define i1mach		Rf_i1mach
#define gammalims	Rf_gammalims
#define lfastchoose	Rf_lfastchoose
#define lgammacor	Rf_lgammacor
#define stirlerr       	Rf_stirlerr

	/* Chebyshev Series */

int	chebyshev_init(double*, int, double);
double	chebyshev_eval(double, const double *, const int);

	/* Gamma and Related Functions */

void	gammalims(double*, double*);
double	lgammacor(double); /* log(gamma) correction */
double  stirlerr(double);  /* Stirling expansion "error" */

double	lfastchoose(double, double);

double  bd0(double, double);

/* Consider adding these two to the API (Rmath.h): */
double	dbinom_raw(double, double, double, double, int);
double	dpois_raw (double, double, int);
double  pnchisq_raw(double, double, double, double, double, int);

int	i1mach(int);

/* From toms708.c */
void bratio(double a, double b, double x, double y, 
	    double *w, double *w1, int *ierr);


#endif /* MATHLIB_PRIVATE_H */


	/* Utilities for `dpq' handling (density/probability/quantile) */

/* give_log in "d";  log_p in "p" & "q" : */
#define give_log log_p
							/* "DEFAULT" */
							/* --------- */
#define R_D__0	(log_p ? ML_NEGINF : 0.)		/* 0 */
#define R_D__1	(log_p ? 0. : 1.)			/* 1 */
#define R_DT_0	(lower_tail ? R_D__0 : R_D__1)		/* 0 */
#define R_DT_1	(lower_tail ? R_D__1 : R_D__0)		/* 1 */

#define R_D_Lval(p)	(lower_tail ? (p) : (1 - (p)))	/*  p  */
#define R_D_Cval(p)	(lower_tail ? (1 - (p)) : (p))	/*  1 - p */

#define R_D_val(x)	(log_p	? log(x) : (x))		/*  x  in pF(x,..) */
#define R_D_qIv(p)	(log_p	? exp(p) : (p))		/*  p  in qF(p,..) */
#define R_D_exp(x)	(log_p	?  (x)	 : exp(x))	/* exp(x) */
#define R_D_log(p)	(log_p	?  (p)	 : log(p))	/* log(p) */
#define R_D_Clog(p)	(log_p	? log1p(-(p)) : (1 - (p)))/* [log](1-p) */

/* log(1-exp(x)):  R_D_LExp(x) == (log1p(- R_D_qIv(x))) but even more stable:*/
#define R_D_LExp(x)     (log_p ? R_Log1_Exp(x) : log1p(-x))

/*till 1.8.x:
 * #define R_DT_val(x)	R_D_val(R_D_Lval(x))
 * #define R_DT_Cval(x)	R_D_val(R_D_Cval(x)) */
#define R_DT_val(x)	(lower_tail ? R_D_val(x)  : R_D_Clog(x))
#define R_DT_Cval(x)	(lower_tail ? R_D_Clog(x) : R_D_val(x))

/*#define R_DT_qIv(p)	R_D_Lval(R_D_qIv(p))		 *  p  in qF ! */
#define R_DT_qIv(p)	(log_p ? (lower_tail ? exp(p) : - expm1(p)) \
			       : R_D_Lval(p))

/*#define R_DT_CIv(p)	R_D_Cval(R_D_qIv(p))		 *  1 - p in qF */
#define R_DT_CIv(p)	(log_p ? (lower_tail ? -expm1(p) : exp(p)) \
			       : R_D_Cval(p))

#define R_DT_exp(x)	R_D_exp(R_D_Lval(x))		/* exp(x) */
#define R_DT_Cexp(x)	R_D_exp(R_D_Cval(x))		/* exp(1 - x) */

#define R_DT_log(p)	(lower_tail? R_D_log(p) : R_D_LExp(p))/* log(p) in qF */
#define R_DT_Clog(p)	(lower_tail? R_D_LExp(p): R_D_log(p))/* log(1-p) in qF*/
#define R_DT_Log(p)	(lower_tail? (p) : R_Log1_Exp(p))
/* ==   R_DT_log when we already "know" log_p == TRUE :*/

#define R_Q_P01_check(p)			\
    if ((log_p	&& p > 0) ||			\
	(!log_p && (p < 0 || p > 1)) )		\
	ML_ERR_return_NAN

/* additions for density functions (C.Loader) */
#define R_D_fexp(f,x)     (give_log ? -0.5*log(f)+(x) : exp(x)/sqrt(f))
#define R_D_forceint(x)   floor((x) + 0.5)
#define R_D_nonint(x) 	  (fabs((x) - floor((x)+0.5)) > 1e-7)
/* [neg]ative or [non int]eger : */
#define R_D_negInonint(x) (x < 0. || R_D_nonint(x))

#define R_D_nonint_check(x) 				\
   if(R_D_nonint(x)) {					\
	MATHLIB_WARNING("non-integer x = %f", x);	\
	return R_D__0;					\
   }

double igraph_qnorm5(double p, double mu, double sigma, int lower_tail, int log_p)
{
    double p_, q, r, val;

#ifdef IEEE_754
    if (ISNAN(p) || ISNAN(mu) || ISNAN(sigma))
	return p + mu + sigma;
#endif
    if (p == R_DT_0)	return ML_NEGINF;
    if (p == R_DT_1)	return ML_POSINF;
    R_Q_P01_check(p);

    if(sigma  < 0)	ML_ERR_return_NAN;
    if(sigma == 0)	return mu;

    p_ = R_DT_qIv(p);/* real lower_tail prob. p */
    q = p_ - 0.5;

/*-- use AS 241 --- */
/* double ppnd16_(double *p, long *ifault)*/
/*      ALGORITHM AS241  APPL. STATIST. (1988) VOL. 37, NO. 3

        Produces the normal deviate Z corresponding to a given lower
        tail area of P; Z is accurate to about 1 part in 10**16.

        (original fortran code used PARAMETER(..) for the coefficients
         and provided hash codes for checking them...)
*/
    if (fabs(q) <= .425) {/* 0.075 <= p <= 0.925 */
        r = .180625 - q * q;
	val =
            q * (((((((r * 2509.0809287301226727 +
                       33430.575583588128105) * r + 67265.770927008700853) * r +
                     45921.953931549871457) * r + 13731.693765509461125) * r +
                   1971.5909503065514427) * r + 133.14166789178437745) * r +
                 3.387132872796366608)
            / (((((((r * 5226.495278852854561 +
                     28729.085735721942674) * r + 39307.89580009271061) * r +
                   21213.794301586595867) * r + 5394.1960214247511077) * r +
                 687.1870074920579083) * r + 42.313330701600911252) * r + 1.);
    }
    else { /* closer than 0.075 from {0,1} boundary */

	/* r = min(p, 1-p) < 0.075 */
	if (q > 0)
	    r = R_DT_CIv(p);/* 1-p */
	else
	    r = p_;/* = R_DT_Iv(p) ^=  p */

	r = sqrt(- ((log_p &&
		     ((lower_tail && q <= 0) || (!lower_tail && q > 0))) ?
		    p : /* else */ log(r)));
        /* r = sqrt(-log(r))  <==>  min(p, 1-p) = exp( - r^2 ) */

        if (r <= 5.) { /* <==> min(p,1-p) >= exp(-25) ~= 1.3888e-11 */
            r += -1.6;
            val = (((((((r * 7.7454501427834140764e-4 +
                       .0227238449892691845833) * r + .24178072517745061177) *
                     r + 1.27045825245236838258) * r +
                    3.64784832476320460504) * r + 5.7694972214606914055) *
                  r + 4.6303378461565452959) * r +
                 1.42343711074968357734)
                / (((((((r *
                         1.05075007164441684324e-9 + 5.475938084995344946e-4) *
                        r + .0151986665636164571966) * r +
                       .14810397642748007459) * r + .68976733498510000455) *
                     r + 1.6763848301838038494) * r +
                    2.05319162663775882187) * r + 1.);
        }
        else { /* very close to  0 or 1 */
            r += -5.;
            val = (((((((r * 2.01033439929228813265e-7 +
                       2.71155556874348757815e-5) * r +
                      .0012426609473880784386) * r + .026532189526576123093) *
                    r + .29656057182850489123) * r +
                   1.7848265399172913358) * r + 5.4637849111641143699) *
                 r + 6.6579046435011037772)
                / (((((((r *
                         2.04426310338993978564e-15 + 1.4215117583164458887e-7)*
                        r + 1.8463183175100546818e-5) * r +
                       7.868691311456132591e-4) * r + .0148753612908506148525)
                     * r + .13692988092273580531) * r +
                    .59983220655588793769) * r + 1.);
        }

	if(q < 0.0)
	    val = -val;
        /* return (q >= 0.)? r : -r ;*/
    }
    return mu + sigma * val;
}

double fsign(double x, double y)
{
#ifdef IEEE_754
    if (ISNAN(x) || ISNAN(y))
        return x + y;
#endif
    return ((y >= 0) ? fabs(x) : -fabs(x));
}

int imax2(int x, int y)
{
    return (x < y) ? y : x;
}

int imin2(int x, int y)
{
    return (x < y) ? x : y;
}

#ifdef HAVE_WORKING_ISFINITE
/* isfinite is defined in <math.h> according to C99 */
# define R_FINITE(x)    isfinite(x)
#elif HAVE_WORKING_FINITE
/* include header needed to define finite() */
#  ifdef HAVE_IEEE754_H
#   include <ieee754.h>         /* newer Linuxen */
#  else
#   ifdef HAVE_IEEEFP_H
#    include <ieeefp.h>         /* others [Solaris], .. */
#   endif
#  endif
# define R_FINITE(x)    finite(x)
#else
# define R_FINITE(x)    R_finite(x)
#endif

int R_finite(double x)
{
#ifdef HAVE_WORKING_ISFINITE
    return isfinite(x);
#elif HAVE_WORKING_FINITE
    return finite(x);
#else
/* neither finite nor isfinite work. Do we really need the AIX exception? */
# ifdef _AIX
#  include <fp.h>
     return FINITE(x);
# elif defined(_MSC_VER)
     return _finite(x);
#else
    return (!isnan(x) & (x != 1/0.0) & (x != -1.0/0.0));
# endif
#endif
}

int R_isnancpp(double x)
{
   return (isnan(x)!=0);
}

#ifdef __cplusplus
  int R_isnancpp(double); /* in arithmetic.c */
#  define ISNAN(x)     R_isnancpp(x)
#else
#  define ISNAN(x)     (isnan(x)!=0)
#endif

double igraph_norm_rand(void) {
  
  double u1;

#define BIG 134217728 /* 2^27 */
  /* unif_rand() alone is not of high enough precision */
  u1 = RNG_UNIF(0,1);
  u1 = (int)(BIG*u1) + RNG_UNIF(0,1);
  return igraph_qnorm5(u1/BIG, 0.0, 1.0, 1, 0);
}

/*
 *  Mathlib : A C Library of Special Functions
 *  Copyright (C) 1998 Ross Ihaka
 *  Copyright (C) 2000-2002 the R Development Core Team
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA.
 *
 *  SYNOPSIS
 *
 *    #include <Rmath.h>
 *    double exp_rand(void);
 *
 *  DESCRIPTION
 *
 *    Random variates from the standard exponential distribution.
 *
 *  REFERENCE
 *
 *    Ahrens, J.H. and Dieter, U. (1972).
 *    Computer methods for sampling from the exponential and
 *    normal distributions.
 *    Comm. ACM, 15, 873-882.
 */

double igraph_exp_rand(void)
{
    /* q[k-1] = sum(log(2)^k / k!)  k=1,..,n, */
    /* The highest n (here 8) is determined by q[n-1] = 1.0 */
    /* within standard precision */
    const double q[] =
    {
	0.6931471805599453,
	0.9333736875190459,
	0.9888777961838675,
	0.9984959252914960,
	0.9998292811061389,
	0.9999833164100727,
	0.9999985691438767,
	0.9999998906925558,
	0.9999999924734159,
	0.9999999995283275,
	0.9999999999728814,
	0.9999999999985598,
	0.9999999999999289,
	0.9999999999999968,
	0.9999999999999999,
	1.0000000000000000
    };
    double a, u, ustar, umin;
    int i;

    a = 0.;
    /* precaution if u = 0 is ever returned */
    u = RNG_UNIF01();
    while(u <= 0.0 || u >= 1.0) u = RNG_UNIF01();
    for (;;) {
	u += u;
	if (u > 1.0)
	    break;
	a += q[0];
    }
    u -= 1.;

    if (u <= q[0])
	return a + u;

    i = 0;
    ustar = RNG_UNIF01();
    umin = ustar;
    do {
	ustar = RNG_UNIF01();
	if (ustar < umin)
	    umin = ustar;
	i++;
    } while (u > q[i]);
    return a + umin * q[0];
}

/*
 *  Mathlib : A C Library of Special Functions
 *  Copyright (C) 1998 Ross Ihaka
 *  Copyright (C) 2000-2001 The R Development Core Team
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA.
 *
 *  SYNOPSIS
 *
 *    #include <Rmath.h>
 *    double rpois(double lambda)
 *
 *  DESCRIPTION
 *
 *    Random variates from the Poisson distribution.
 *
 *  REFERENCE
 *
 *    Ahrens, J.H. and Dieter, U. (1982).
 *    Computer generation of Poisson deviates
 *    from modified normal distributions.
 *    ACM Trans. Math. Software 8, 163-179.
 */

#define a0	-0.5
#define a1	 0.3333333
#define a2	-0.2500068
#define a3	 0.2000118
#define a4	-0.1661269
#define a5	 0.1421878
#define a6	-0.1384794
#define a7	 0.1250060

#define one_7	0.1428571428571428571
#define one_12	0.0833333333333333333
#define one_24	0.0416666666666666667

#define repeat for(;;)

#define FALSE 0
#define TRUE  1
#define M_1_SQRT_2PI    0.398942280401432677939946059934     /* 1/sqrt(2pi) */

double igraph_rpois(double mu)
{
    /* Factorial Table (0:9)! */
    const double fact[10] =
    {
	1., 1., 2., 6., 24., 120., 720., 5040., 40320., 362880.
    };

    /* These are static --- persistent between calls for same mu : */
    static int l, m;

    static double b1, b2, c, c0, c1, c2, c3;
    static double pp[36], p0, p, q, s, d, omega;
    static double big_l;/* integer "w/o overflow" */
    static double muprev = 0., muprev2 = 0.;/*, muold	 = 0.*/

    /* Local Vars  [initialize some for -Wall]: */
    double del, difmuk= 0., E= 0., fk= 0., fx, fy, g, px, py, t, u= 0., v, x;
    double pois = -1.;
    int k, kflag, big_mu, new_big_mu = FALSE;

    if (!R_FINITE(mu))
	ML_ERR_return_NAN;

    if (mu <= 0.)
	return 0.;

    big_mu = mu >= 10.;
    if(big_mu)
	new_big_mu = FALSE;

    if (!(big_mu && mu == muprev)) {/* maybe compute new persistent par.s */

	if (big_mu) {
	    new_big_mu = TRUE;
	    /* Case A. (recalculation of s,d,l	because mu has changed):
	     * The poisson probabilities pk exceed the discrete normal
	     * probabilities fk whenever k >= m(mu).
	     */
	    muprev = mu;
	    s = sqrt(mu);
	    d = 6. * mu * mu;
	    big_l = floor(mu - 1.1484);
	    /* = an upper bound to m(mu) for all mu >= 10.*/
	}
	else { /* Small mu ( < 10) -- not using normal approx. */

	    /* Case B. (start new table and calculate p0 if necessary) */

	    /*muprev = 0.;-* such that next time, mu != muprev ..*/
	    if (mu != muprev) {
		muprev = mu;
		m = imax2(1, (int) mu);
		l = 0; /* pp[] is already ok up to pp[l] */
		q = p0 = p = exp(-mu);
	    }

	    repeat {
		/* Step U. uniform sample for inversion method */
		u = RNG_UNIF01();
		if (u <= p0)
		    return 0.;

		/* Step T. table comparison until the end pp[l] of the
		   pp-table of cumulative poisson probabilities
		   (0.458 > ~= pp[9](= 0.45792971447) for mu=10 ) */
		if (l != 0) {
		    for (k = (u <= 0.458) ? 1 : imin2(l, m);  k <= l; k++)
			if (u <= pp[k])
			    return (double)k;
		    if (l == 35) /* u > pp[35] */
			continue;
		}
		/* Step C. creation of new poisson
		   probabilities p[l..] and their cumulatives q =: pp[k] */
		l++;
		for (k = l; k <= 35; k++) {
		    p *= mu / k;
		    q += p;
		    pp[k] = q;
		    if (u <= q) {
			l = k;
			return (double)k;
		    }
		}
		l = 35;
	    } /* end(repeat) */
	}/* mu < 10 */

    } /* end {initialize persistent vars} */

/* Only if mu >= 10 : ----------------------- */

    /* Step N. normal sample */
    g = mu + s * igraph_norm_rand();/* norm_rand() ~ N(0,1), standard normal */

    if (g >= 0.) {
	pois = floor(g);
	/* Step I. immediate acceptance if pois is large enough */
	if (pois >= big_l)
	    return pois;
	/* Step S. squeeze acceptance */
	fk = pois;
	difmuk = mu - fk;
	u = RNG_UNIF01(); /* ~ U(0,1) - sample */
	if (d * u >= difmuk * difmuk * difmuk)
	    return pois;
    }

    /* Step P. preparations for steps Q and H.
       (recalculations of parameters if necessary) */

    if (new_big_mu || mu != muprev2) {
        /* Careful! muprev2 is not always == muprev
	   because one might have exited in step I or S
	   */
        muprev2 = mu;
	omega = M_1_SQRT_2PI / s;
	/* The quantities b1, b2, c3, c2, c1, c0 are for the Hermite
	 * approximations to the discrete normal probabilities fk. */

	b1 = one_24 / mu;
	b2 = 0.3 * b1 * b1;
	c3 = one_7 * b1 * b2;
	c2 = b2 - 15. * c3;
	c1 = b1 - 6. * b2 + 45. * c3;
	c0 = 1. - b1 + 3. * b2 - 15. * c3;
	c = 0.1069 / mu; /* guarantees majorization by the 'hat'-function. */
    }

    if (g >= 0.) {
	/* 'Subroutine' F is called (kflag=0 for correct return) */
	kflag = 0;
	goto Step_F;
    }


    repeat {
	/* Step E. Exponential Sample */

	E = igraph_exp_rand();	/* ~ Exp(1) (standard exponential) */

	/*  sample t from the laplace 'hat'
	    (if t <= -0.6744 then pk < fk for all mu >= 10.) */
	u = 2 * RNG_UNIF01() - 1.;
	t = 1.8 + fsign(E, u);
	if (t > -0.6744) {
	    pois = floor(mu + s * t);
	    fk = pois;
	    difmuk = mu - fk;

	    /* 'subroutine' F is called (kflag=1 for correct return) */
	    kflag = 1;

	  Step_F: /* 'subroutine' F : calculation of px,py,fx,fy. */

	    if (pois < 10) { /* use factorials from table fact[] */
		px = -mu;
		py = pow(mu, pois) / fact[(int)pois];
	    }
	    else {
		/* Case pois >= 10 uses polynomial approximation
		   a0-a7 for accuracy when advisable */
		del = one_12 / fk;
		del = del * (1. - 4.8 * del * del);
		v = difmuk / fk;
		if (fabs(v) <= 0.25)
		    px = fk * v * v * (((((((a7 * v + a6) * v + a5) * v + a4) *
					  v + a3) * v + a2) * v + a1) * v + a0)
			- del;
		else /* |v| > 1/4 */
		    px = fk * log(1. + v) - difmuk - del;
		py = M_1_SQRT_2PI / sqrt(fk);
	    }
	    x = (0.5 - difmuk) / s;
	    x *= x;/* x^2 */
	    fx = -0.5 * x;
	    fy = omega * (((c3 * x + c2) * x + c1) * x + c0);
	    if (kflag > 0) {
		/* Step H. Hat acceptance (E is repeated on rejection) */
		if (c * fabs(u) <= py * exp(px + E) - fy * exp(fx + E))
		    break;
	    } else
		/* Step Q. Quotient acceptance (rare case) */
		if (fy - u * fy <= py * exp(px - fx))
		    break;
	}/* t > -.67.. */
    }
    return pois;
}

double igraph_rgeom(double p) {
    if (ISNAN(p) || p <= 0 || p > 1) ML_ERR_return_NAN;

    return igraph_rpois(igraph_exp_rand() * ((1 - p) / p));
}

/* This is from nmath/rbinom.c */

#define repeat for(;;)

double igraph_rbinom(double nin, double pp)
{
    /* FIXME: These should become THREAD_specific globals : */

    static double c, fm, npq, p1, p2, p3, p4, qn;
    static double xl, xll, xlr, xm, xr;

    static double psave = -1.0;
    static int nsave = -1;
    static int m;

    double f, f1, f2, u, v, w, w2, x, x1, x2, z, z2;
    double p, q, np, g, r, al, alv, amaxp, ffm, ynorm;
    int i,ix,k, n;

    if (!R_FINITE(nin)) ML_ERR_return_NAN;
    n = floor(nin + 0.5);
    if (n != nin) ML_ERR_return_NAN;

    if (!R_FINITE(pp) ||
	/* n=0, p=0, p=1 are not errors <TSL>*/
	n < 0 || pp < 0. || pp > 1.)	ML_ERR_return_NAN;

    if (n == 0 || pp == 0.) return 0;
    if (pp == 1.) return n;

    p = fmin(pp, 1. - pp);
    q = 1. - p;
    np = n * p;
    r = p / q;
    g = r * (n + 1);

    /* Setup, perform only when parameters change [using static (globals): */

    /* FIXING: Want this thread safe
       -- use as little (thread globals) as possible
    */
    if (pp != psave || n != nsave) {
	psave = pp;
	nsave = n;
	if (np < 30.0) {
	    /* inverse cdf logic for mean less than 30 */
	    qn = pow(q, (double) n);
	    goto L_np_small;
	} else {
	    ffm = np + p;
	    m = ffm;
	    fm = m;
	    npq = np * q;
	    p1 = (int)(2.195 * sqrt(npq) - 4.6 * q) + 0.5;
	    xm = fm + 0.5;
	    xl = xm - p1;
	    xr = xm + p1;
	    c = 0.134 + 20.5 / (15.3 + fm);
	    al = (ffm - xl) / (ffm - xl * p);
	    xll = al * (1.0 + 0.5 * al);
	    al = (xr - ffm) / (xr * q);
	    xlr = al * (1.0 + 0.5 * al);
	    p2 = p1 * (1.0 + c + c);
	    p3 = p2 + c / xll;
	    p4 = p3 + c / xlr;
	}
    } else if (n == nsave) {
	if (np < 30.0)
	    goto L_np_small;
    }

    /*-------------------------- np = n*p >= 30 : ------------------- */
    repeat {
      u = RNG_UNIF01() * p4;
      v = RNG_UNIF01();
      /* triangular region */
      if (u <= p1) {
	  ix = xm - p1 * v + u;
	  goto finis;
      }
      /* parallelogram region */
      if (u <= p2) {
	  x = xl + (u - p1) / c;
	  v = v * c + 1.0 - fabs(xm - x) / p1;
	  if (v > 1.0 || v <= 0.)
	      continue;
	  ix = x;
      } else {
	  if (u > p3) {	/* right tail */
	      ix = xr - log(v) / xlr;
	      if (ix > n)
		  continue;
	      v = v * (u - p3) * xlr;
	  } else {/* left tail */
	      ix = xl + log(v) / xll;
	      if (ix < 0)
		  continue;
	      v = v * (u - p2) * xll;
	  }
      }
      /* determine appropriate way to perform accept/reject test */
      k = abs(ix - m);
      if (k <= 20 || k >= npq / 2 - 1) {
	  /* explicit evaluation */
	  f = 1.0;
	  if (m < ix) {
	      for (i = m + 1; i <= ix; i++)
		  f *= (g / i - r);
	  } else if (m != ix) {
	      for (i = ix + 1; i <= m; i++)
		  f /= (g / i - r);
	  }
	  if (v <= f)
	      goto finis;
      } else {
	  /* squeezing using upper and lower bounds on log(f(x)) */
	  amaxp = (k / npq) * ((k * (k / 3. + 0.625) + 0.1666666666666) / npq + 0.5);
	  ynorm = -k * k / (2.0 * npq);
	  alv = log(v);
	  if (alv < ynorm - amaxp)
	      goto finis;
	  if (alv <= ynorm + amaxp) {
	      /* stirling's formula to machine accuracy */
	      /* for the final acceptance/rejection test */
	      x1 = ix + 1;
	      f1 = fm + 1.0;
	      z = n + 1 - fm;
	      w = n - ix + 1.0;
	      z2 = z * z;
	      x2 = x1 * x1;
	      f2 = f1 * f1;
	      w2 = w * w;
	      if (alv <= xm * log(f1 / x1) + (n - m + 0.5) * log(z / w) + (ix - m) * log(w * p / (x1 * q)) + (13860.0 - (462.0 - (132.0 - (99.0 - 140.0 / f2) / f2) / f2) / f2) / f1 / 166320.0 + (13860.0 - (462.0 - (132.0 - (99.0 - 140.0 / z2) / z2) / z2) / z2) / z / 166320.0 + (13860.0 - (462.0 - (132.0 - (99.0 - 140.0 / x2) / x2) / x2) / x2) / x1 / 166320.0 + (13860.0 - (462.0 - (132.0 - (99.0 - 140.0 / w2) / w2) / w2) / w2) / w / 166320.)
		  goto finis;
	  }
      }
  }

 L_np_small:
    /*---------------------- np = n*p < 30 : ------------------------- */

  repeat {
     ix = 0;
     f = qn;
     u = RNG_UNIF01();
     repeat {
	 if (u < f)
	     goto finis;
	 if (ix > 110)
	     break;
	 u -= f;
	 ix++;
	 f *= (g / ix - r);
     }
  }
 finis:
    if (psave > 0.5)
	 ix = n - ix;
  return (double)ix;
}

#endif

/**********************************************************
 * Testing purposes                                       *
 *********************************************************/

/* int main() { */

/*   int i; */

/*   for (i=0; i<100; i++) { */
/*     printf("%i ", RNG_INTEGER(0,10)); */
/*   } */
/*   printf("\n"); */

/*   for (i=0; i<100; i++) { */
/*     printf("%f ", RNG_UNIF(0,1)); */
/*   } */
/*   printf("\n"); */

/*   for (i=0; i<100; i++) { */
/*     printf("%f ", RNG_NORMAL(0,5)); */
/*   } */
/*   printf("\n"); */

/* } */

  
