/* specfunc/zeta.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2004 Gerard Jungman
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Author:  G. Jungman */

/* This file was taken from the GNU Scientific Library. Some modifications
 * were done in order to make it independent from the rest of GSL
 */

/*
#include <config.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_elementary.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_pow_int.h>
#include <gsl/gsl_sf_zeta.h>

#include "error.h"

#include "chebyshev.h"
#include "cheb_eval.c"
*/

#include <math.h>
#include <stdio.h>
#include "error.h"

/*-*-*-*-*-*-*-*-*-*- From gsl_machine.h -*-*-*-*-*-*-*-*-*-*-*-*-*/

#define GSL_LOG_DBL_MIN   (-7.0839641853226408e+02)
#define GSL_LOG_DBL_MAX    7.0978271289338397e+02
#define GSL_DBL_EPSILON        2.2204460492503131e-16

/*-*-*-*-*-*-*-*-*-* From gsl_sf_result.h *-*-*-*-*-*-*-*-*-*-*-*/

struct gsl_sf_result_struct {
  double val;
  double err;
};
typedef struct gsl_sf_result_struct gsl_sf_result;

/*-*-*-*-*-*-*-*-*-*-*-* Private Section *-*-*-*-*-*-*-*-*-*-*-*/

/* coefficients for Maclaurin summation in hzeta()
 * B_{2j}/(2j)!
 */
static double hzeta_c[15] = {
  1.00000000000000000000000000000,
  0.083333333333333333333333333333,
 -0.00138888888888888888888888888889,
  0.000033068783068783068783068783069,
 -8.2671957671957671957671957672e-07,
  2.0876756987868098979210090321e-08,
 -5.2841901386874931848476822022e-10,
  1.3382536530684678832826980975e-11,
 -3.3896802963225828668301953912e-13,
  8.5860620562778445641359054504e-15,
 -2.1748686985580618730415164239e-16,
  5.5090028283602295152026526089e-18,
 -1.3954464685812523340707686264e-19,
  3.5347070396294674716932299778e-21,
 -8.9535174270375468504026113181e-23
};

/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/

static int gsl_sf_hzeta_e(const double s, const double q, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(s <= 1.0 || q <= 0.0) {
	PLFIT_ERROR("s must be larger than 1.0 and q must be larger than zero", PLFIT_EINVAL);
  }
  else {
    const double max_bits = 54.0;
    const double ln_term0 = -s * log(q);  

    if(ln_term0 < GSL_LOG_DBL_MIN + 1.0) {
	  PLFIT_ERROR("underflow", PLFIT_UNDRFLOW);
    }
    else if(ln_term0 > GSL_LOG_DBL_MAX - 1.0) {
	  PLFIT_ERROR("overflow", PLFIT_OVERFLOW);
    }
    else if((s > max_bits && q < 1.0) || (s > 0.5*max_bits && q < 0.25)) {
      result->val = pow(q, -s);
      result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
      return PLFIT_SUCCESS;
    }
    else if(s > 0.5*max_bits && q < 1.0) {
      const double p1 = pow(q, -s);
      const double p2 = pow(q/(1.0+q), s);
      const double p3 = pow(q/(2.0+q), s);
      result->val = p1 * (1.0 + p2 + p3);
      result->err = GSL_DBL_EPSILON * (0.5*s + 2.0) * fabs(result->val);
      return PLFIT_SUCCESS;
    }
    else {
      /* Euler-Maclaurin summation formula 
       * [Moshier, p. 400, with several typo corrections]
       */
      const int jmax = 12;
      const int kmax = 10;
      int j, k;
      const double pmax  = pow(kmax + q, -s);
      double scp = s;
      double pcp = pmax / (kmax + q);
      double ans = pmax*((kmax+q)/(s-1.0) + 0.5);

      for(k=0; k<kmax; k++) {
        ans += pow(k + q, -s);
      }

      for(j=0; j<=jmax; j++) {
        double delta = hzeta_c[j+1] * scp * pcp;
        ans += delta;
        if(fabs(delta/ans) < 0.5*GSL_DBL_EPSILON) break;
        scp *= (s+2*j+1)*(s+2*j+2);
        pcp /= (kmax + q)*(kmax + q);
      }

      result->val = ans;
      result->err = 2.0 * (jmax + 1.0) * GSL_DBL_EPSILON * fabs(ans);
      return PLFIT_SUCCESS;
    }
  }
}

/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/

double gsl_sf_hzeta(const double s, const double a)
{
  gsl_sf_result result;
  gsl_sf_hzeta_e(s, a, &result);
  return result.val;
}

