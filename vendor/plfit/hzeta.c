/* vim:set ts=4 sw=2 sts=2 et: */

/* This file was imported from a private scientific library
 * based on GSL coined Home Scientific Libray (HSL) by its author
 * Jerome Benoit; this very material is itself inspired from the
 * material written by G. Jungan and distributed by GSL.
 * Ultimately, some modifications were done in order to render the
 * imported material independent from the rest of GSL.
 */

/* `hsl/specfunc/hzeta.c' C source file
// HSL - Home Scientific Library
// Copyright (C) 2017-2018  Jerome Benoit
//
// HSL is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 2
// of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/

/*
// The material in this file is mainly inspired by the material written by
// G. Jungan and distributed under GPLv2 by the GNU Scientific Library (GSL)
// ( https://www.gnu.org/software/gsl/ [specfunc/zeta.c]), itself inspired by
// the material written by Moshier and distributed in the Cephes Mathematical
// Library ( http://www.moshier.net/ [zeta.c]).
//
// More specifically, hsl_sf_hzeta_e is a slightly modifed clone of
// gsl_sf_hzeta_e as found in GSL 2.4; the remaining is `inspired by'.
// [Sooner or later a _Working_Note_ may be deposited at ResearchGate
// ( https://www.researchgate.net/profile/Jerome_Benoit )]
*/

/* Author:  Jerome G. Benoit < jgmbenoit _at_ rezozer _dot_ net > */

#ifdef _MSC_VER
#define _USE_MATH_DEFINES
#endif

#include <math.h>
#include <stdio.h>
#include "hzeta.h"
#include "error.h"
#include "platform.h"

/* imported from gsl_machine.h */

#define GSL_LOG_DBL_MIN (-7.0839641853226408e+02)
#define GSL_LOG_DBL_MAX 7.0978271289338397e+02
#define GSL_DBL_EPSILON 2.2204460492503131e-16

/* imported from gsl_math.h */

#ifndef M_LOG2E
#define M_LOG2E    1.44269504088896340735992468100      /* log_2 (e) */
#endif

/* imported from gsl_sf_result.h */

struct gsl_sf_result_struct {
  double val;
  double err;
};
typedef struct gsl_sf_result_struct gsl_sf_result;

/* imported and adapted from hsl/specfunc/specfunc_def.h */

#define HSL_SF_EVAL_RESULT(FnE) \
	gsl_sf_result result; \
	FnE ; \
	return (result.val);

#define HSL_SF_EVAL_TUPLE_RESULT(FnET) \
	gsl_sf_result result0; \
	gsl_sf_result result1; \
	FnET ; \
	*tuple1=result1.val; \
	*tuple0=result0.val; \
	return (result0.val);

/* */


#define HSL_SF_HZETA_EULERMACLAURIN_SERIES_SHIFT 10
#define HSL_SF_HZETA_EULERMACLAURIN_SERIES_ORDER 32

#define HSL_SF_LNHZETA_EULERMACLAURIN_SERIES_SHIFT_MAX 256

// B_{2j}/(2j)
static
double hsl_sf_hzeta_eulermaclaurin_series_coeffs[HSL_SF_HZETA_EULERMACLAURIN_SERIES_ORDER+1]={
	+1.0,
	+1.0/12.0,
	-1.0/720.0,
	+1.0/30240.0,
	-1.0/1209600.0,
	+1.0/47900160.0,
	-691.0/1307674368000.0,
	+1.0/74724249600.0,
	-3.38968029632258286683019539125e-13,
	+8.58606205627784456413590545043e-15,
	-2.17486869855806187304151642387e-16,
	+5.50900282836022951520265260890e-18,
	-1.39544646858125233407076862641e-19,
	+3.53470703962946747169322997780e-21,
	-8.95351742703754685040261131811e-23,
	+2.26795245233768306031095073887e-24,
	-5.74479066887220244526388198761e-26,
	+1.45517247561486490186626486727e-27,
	-3.68599494066531017818178247991e-29,
	+9.33673425709504467203255515279e-31,
	-2.36502241570062993455963519637e-32,
	+5.99067176248213430465991239682e-34,
	-1.51745488446829026171081313586e-35,
	+3.84375812545418823222944529099e-37,
	-9.73635307264669103526762127925e-39,
	+2.46624704420068095710640028029e-40,
	-6.24707674182074369314875679472e-42,
	+1.58240302446449142975108170683e-43,
	-4.00827368594893596853001219052e-45,
	+1.01530758555695563116307139454e-46,
	-2.57180415824187174992481940976e-48,
	+6.51445603523381493155843485864e-50,
	-1.65013099068965245550609878048e-51
	}; // hsl_sf_hzeta_eulermaclaurin_series_coeffs

// 4\zeta(2j)/(2\pi)^(2j)
static
double hsl_sf_hzeta_eulermaclaurin_series_majorantratios[HSL_SF_HZETA_EULERMACLAURIN_SERIES_ORDER+1]={
	-2.0,
	+1.0/6.0,
	+1.0/360.0,
	+1.0/15120.0,
	+1.0/604800.0,
	+1.0/23950080.0,
	+691.0/653837184000.0,
	+1.0/37362124800.0,
	+3617.0/5335311421440000.0,
	+1.71721241125556891282718109009e-14,
	+4.34973739711612374608303284773e-16,
	+1.10180056567204590304053052178e-17,
	+2.79089293716250466814153725281e-19,
	+7.06941407925893494338645995561e-21,
	+1.79070348540750937008052226362e-22,
	+4.53590490467536612062190147774e-24,
	+1.14895813377444048905277639752e-25,
	+2.91034495122972980373252973454e-27,
	+7.37198988133062035636356495982e-29,
	+1.86734685141900893440651103056e-30,
	+4.73004483140125986911927039274e-32,
	+1.19813435249642686093198247936e-33,
	+3.03490976893658052342162627173e-35,
	+7.68751625090837646445889058198e-37,
	+1.94727061452933820705352425585e-38,
	+4.93249408840136191421280056051e-40,
	+1.24941534836414873862975135893e-41,
	+3.16480604892898285950216341362e-43,
	+8.01654737189787193706002438098e-45,
	+2.03061517111391126232614278906e-46,
	+5.14360831648374349984963881946e-48,
	+1.30289120704676298631168697172e-49,
	+3.30026198137930491101219756091e-51
	}; // hsl_sf_hzeta_eulermaclaurin_series_majorantratios


extern
int hsl_sf_hzeta_e(const double s, const double q, gsl_sf_result * result) {

	/* CHECK_POINTER(result) */

	if ((s <= 1.0) || (q <= 0.0)) {
		PLFIT_ERROR("s must be larger than 1.0 and q must be larger than zero", PLFIT_EINVAL);
		}
	else {
		const double max_bits=54.0; // max_bits=\lceil{s}\rceil with \zeta(s,2)=\zeta(s)-1=GSL_DBL_EPSILON
		const double ln_term0=-s*log(q);
		if (ln_term0 < GSL_LOG_DBL_MIN+1.0) {
			PLFIT_ERROR("underflow", PLFIT_UNDRFLOW);
			}
		else if (GSL_LOG_DBL_MAX-1.0 < ln_term0) {
			PLFIT_ERROR("overflow", PLFIT_OVERFLOW);
			}
#if 1
		else if (((max_bits < s) && (q < 1.0)) || ((0.5*max_bits < s) && (q < 0.25))) {
			result->val=pow(q,-s);
			result->err=2.0*GSL_DBL_EPSILON*fabs(result->val);
			return (PLFIT_SUCCESS);
			}
		else if ((0.5*max_bits < s) && (q < 1.0)) {
			const double a0=pow(q,-s);
			const double p1=pow(q/(1.0+q),s);
			const double p2=pow(q/(2.0+q),s);
			const double ans=a0*(1.0+p1+p2);
			result->val=ans;
			result->err=GSL_DBL_EPSILON*(2.0+0.5*s)*fabs(result->val);
			return (PLFIT_SUCCESS);
			}
#endif
		else { // Euler-Maclaurin summation formula
			const double qshift=HSL_SF_HZETA_EULERMACLAURIN_SERIES_SHIFT+q;
			const double inv_qshift=1.0/qshift;
			const double sqr_inv_qshift=inv_qshift*inv_qshift;
			const double inv_sm1=1.0/(s-1.0);
			const double pmax=pow(qshift,-s);
			double terms[HSL_SF_HZETA_EULERMACLAURIN_SERIES_SHIFT+HSL_SF_HZETA_EULERMACLAURIN_SERIES_ORDER+1]={NAN};
			double delta=NAN;
			double tscp=s;
			double scp=tscp;
			double pcp=pmax*inv_qshift;
			double ratio=scp*pcp;
			size_t n=0;
			size_t j=0;
			double ans=0.0;
			double mjr=NAN;

			for(j=0;j<HSL_SF_HZETA_EULERMACLAURIN_SERIES_SHIFT;++j) ans+=(terms[n++]=pow(j+q,-s));
			ans+=(terms[n++]=0.5*pmax);
			ans+=(terms[n++]=pmax*qshift*inv_sm1);
			for(j=1;j<=HSL_SF_HZETA_EULERMACLAURIN_SERIES_ORDER;++j) {
				delta=hsl_sf_hzeta_eulermaclaurin_series_coeffs[j]*ratio;
				ans+=(terms[n++]=delta);
				scp*=++tscp;
				scp*=++tscp;
				pcp*=sqr_inv_qshift;
				ratio=scp*pcp;
				if ((fabs(delta/ans)) < (0.5*GSL_DBL_EPSILON)) break;
				}

			ans=0.0; while (n) ans+=terms[--n];
			mjr=hsl_sf_hzeta_eulermaclaurin_series_majorantratios[j]*ratio;

			result->val=+ans;
			result->err=2.0*((HSL_SF_HZETA_EULERMACLAURIN_SERIES_SHIFT+1.0)*GSL_DBL_EPSILON*fabs(ans)+mjr);
			return (PLFIT_SUCCESS);
			}
		}

}

extern
double hsl_sf_hzeta(const double s, const double q) {
	HSL_SF_EVAL_RESULT(hsl_sf_hzeta_e(s,q,&result)); }

extern
int hsl_sf_hzeta_deriv_e(const double s, const double q, gsl_sf_result * result) {

	/* CHECK_POINTER(result) */

	if ((s <= 1.0) || (q <= 0.0)) {
		PLFIT_ERROR("s must be larger than 1.0 and q must be larger than zero", PLFIT_EINVAL);
		}
	else {
		const double ln_hz_term0=-s*log(q);
		if (ln_hz_term0 < GSL_LOG_DBL_MIN+1.0) {
			PLFIT_ERROR("underflow", PLFIT_UNDRFLOW);
			}
		else if (GSL_LOG_DBL_MAX-1.0 < ln_hz_term0) {
			PLFIT_ERROR("overflow", PLFIT_OVERFLOW);
			}
		else { // Euler-Maclaurin summation formula
			const double qshift=HSL_SF_HZETA_EULERMACLAURIN_SERIES_SHIFT+q;
			const double inv_qshift=1.0/qshift;
			const double sqr_inv_qshift=inv_qshift*inv_qshift;
			const double inv_sm1=1.0/(s-1.0);
			const double pmax=pow(qshift,-s);
			const double lmax=log(qshift);
			double terms[HSL_SF_HZETA_EULERMACLAURIN_SERIES_SHIFT+HSL_SF_HZETA_EULERMACLAURIN_SERIES_ORDER+1]={NAN};
			double delta=NAN;
			double tscp=s;
			double scp=tscp;
			double pcp=pmax*inv_qshift;
			double lcp=lmax-1.0/s;
			double ratio=scp*pcp*lcp;
			double qs=NAN;
			size_t n=0;
			size_t j=0;
			double ans=0.0;
			double mjr=NAN;

			for(j=0,qs=q;j<HSL_SF_HZETA_EULERMACLAURIN_SERIES_SHIFT;++qs,++j) ans+=(terms[n++]=log(qs)*pow(qs,-s));
			ans+=(terms[n++]=0.5*lmax*pmax);
			ans+=(terms[n++]=pmax*qshift*inv_sm1*(lmax+inv_sm1));
			for(j=1;j<=HSL_SF_HZETA_EULERMACLAURIN_SERIES_ORDER;++j) {
				delta=hsl_sf_hzeta_eulermaclaurin_series_coeffs[j]*ratio;
				ans+=(terms[n++]=delta);
				scp*=++tscp; lcp-=1.0/tscp;
				scp*=++tscp; lcp-=1.0/tscp;
				pcp*=sqr_inv_qshift;
				ratio=scp*pcp*lcp;
				if ((fabs(delta/ans)) < (0.5*GSL_DBL_EPSILON)) break;
				}

			ans=0.0; while (n) ans+=terms[--n];
			mjr=hsl_sf_hzeta_eulermaclaurin_series_majorantratios[j]*ratio;

			result->val=-ans;
			result->err=2.0*((HSL_SF_HZETA_EULERMACLAURIN_SERIES_SHIFT+1.0)*GSL_DBL_EPSILON*fabs(ans)+mjr);
			return (PLFIT_SUCCESS);
			}
		}

}

extern
double hsl_sf_hzeta_deriv(const double s, const double q) {
	HSL_SF_EVAL_RESULT(hsl_sf_hzeta_deriv_e(s,q,&result)); }

extern
int hsl_sf_hzeta_deriv2_e(const double s, const double q, gsl_sf_result * result) {

	/* CHECK_POINTER(result) */

	if ((s <= 1.0) || (q <= 0.0)) {
		PLFIT_ERROR("s must be larger than 1.0 and q must be larger than zero", PLFIT_EINVAL);
		}
	else {
		const double ln_hz_term0=-s*log(q);
		if (ln_hz_term0 < GSL_LOG_DBL_MIN+1.0) {
			PLFIT_ERROR("underflow", PLFIT_UNDRFLOW);
			}
		else if (GSL_LOG_DBL_MAX-1.0 < ln_hz_term0) {
			PLFIT_ERROR("overflow", PLFIT_OVERFLOW);
			}
		else { // Euler-Maclaurin summation formula
			const double qshift=HSL_SF_HZETA_EULERMACLAURIN_SERIES_SHIFT+q;
			const double inv_qshift=1.0/qshift;
			const double sqr_inv_qshift=inv_qshift*inv_qshift;
			const double inv_sm1=1.0/(s-1.0);
			const double pmax=pow(qshift,-s);
			const double lmax=log(qshift);
			const double lmax_p_inv_sm1=lmax+inv_sm1;
			const double sqr_inv_sm1=inv_sm1*inv_sm1;
			const double sqr_lmax=lmax*lmax;
			const double sqr_lmax_p_inv_sm1=lmax_p_inv_sm1*lmax_p_inv_sm1;
			double terms[HSL_SF_HZETA_EULERMACLAURIN_SERIES_SHIFT+HSL_SF_HZETA_EULERMACLAURIN_SERIES_ORDER+1]={NAN};
			double delta=NAN;
			double tscp=s;
			double slcp=NAN;
			double plcp=NAN;
			double scp=tscp;
			double pcp=pmax*inv_qshift;
			double lcp=1.0/s-lmax;
			double sqr_lcp=lmax*(lmax-2.0/s);
			double ratio=scp*pcp*sqr_lcp;
			double qs=NAN;
			double lqs=NAN;
			size_t n=0;
			size_t j=0;
			double ans=0.0;
			double mjr=NAN;

			for(j=0,qs=q;j<HSL_SF_HZETA_EULERMACLAURIN_SERIES_SHIFT;++qs,++j) {
				lqs=log(qs);
				ans+=(terms[n++]=lqs*lqs*pow(qs,-s));
				}
			ans+=(terms[n++]=0.5*sqr_lmax*pmax);
			ans+=(terms[n++]=pmax*qshift*inv_sm1*(sqr_lmax_p_inv_sm1+sqr_inv_sm1));
			for(j=1;j<=HSL_SF_HZETA_EULERMACLAURIN_SERIES_ORDER;++j) {
				delta=hsl_sf_hzeta_eulermaclaurin_series_coeffs[j]*ratio;
				ans+=(terms[n++]=delta);
				scp*=++tscp; slcp=plcp=1.0/tscp;
				scp*=++tscp; slcp+=1.0/tscp; plcp/=tscp;
				pcp*=sqr_inv_qshift;
				sqr_lcp+=2.0*(plcp+slcp*lcp);
				ratio=scp*pcp*sqr_lcp;
				if ((fabs(delta/ans)) < (0.5*GSL_DBL_EPSILON)) break;
				lcp+=slcp;
				}

			ans=0.0; while (n) ans+=terms[--n];
			mjr=hsl_sf_hzeta_eulermaclaurin_series_majorantratios[j]*ratio;

			result->val=+ans;
			result->err=2.0*((HSL_SF_HZETA_EULERMACLAURIN_SERIES_SHIFT+1.0)*GSL_DBL_EPSILON*fabs(ans)+mjr);
			return (PLFIT_SUCCESS);
			}
		}

}

extern
double hsl_sf_hzeta_deriv2(const double s, const double q) {
	HSL_SF_EVAL_RESULT(hsl_sf_hzeta_deriv2_e(s,q,&result)); }

static inline
double hsl_sf_hZeta0_zed(const double s, const double q) {
#if 1
	const long double ld_q=(long double)(q);
	const long double ld_s=(long double)(s);
	const long double ld_log1prq=log1pl(1.0L/ld_q);
	const long double ld_epsilon=expm1l(-ld_s*ld_log1prq);
	const long double ld_z=ld_s+(ld_q+0.5L*ld_s+0.5L)*ld_epsilon;
	const double z=(double)(ld_z);
#else
	double z=s+(q+0.5*s+0.5)*expm1(-s*log1p(1.0/q));
#endif
	return (z); }

// Z_{0}(s,a) = a^s \left(\frac{1}{2}+\frac{a}{s-1}\right)^{-1} \zeta(s,a) - 1
// Z_{0}(s,a) = O\left(\frac{(s-1)s}{6a^{2}}\right)
static
int hsl_sf_hZeta0(const double s, const double q, double * value, double * abserror) {
	const double criterion=ceil(10.0*s-q);
	const size_t shift=(criterion<0.0)?0:
		(criterion<HSL_SF_LNHZETA_EULERMACLAURIN_SERIES_SHIFT_MAX)?(size_t)(llrint(criterion)):
			HSL_SF_LNHZETA_EULERMACLAURIN_SERIES_SHIFT_MAX;
	const double qshift=(double)(shift)+q;
	const double inv_qshift=1.0/qshift;
	const double sqr_inv_qshift=inv_qshift*inv_qshift;
	const double sm1=s-1.0;
	double terms[HSL_SF_LNHZETA_EULERMACLAURIN_SERIES_SHIFT_MAX+HSL_SF_HZETA_EULERMACLAURIN_SERIES_ORDER+1]={NAN};
	double delta=NAN;
	double tscp=s;
	double scp=s*sm1;
	double pcp=inv_qshift/(2.0*qshift+sm1);
	double ratio=NAN;
	size_t n=0;
	size_t j=0;
	double ans=0.0;
	double mjr=NAN;

	if (shift) {
		const double hsm1=0.5*sm1;
		const double inv_q=1.0/q;
		const double qphsm1=q+hsm1;
		const double inv_qphsm1=1.0/qphsm1;
		const double qshiftphsm1=qshift+hsm1;
		double qs=q;
		double a=1.0;
		for(j=0;j<shift;) {
			ans+=(terms[n++]=a*hsl_sf_hZeta0_zed(s,qs++)*inv_qphsm1);
			a=exp(-s*log1p((++j)*inv_q));
			}
		pcp*=a*qshiftphsm1*inv_qphsm1;
		}
	ratio=scp*pcp;
	ans+=terms[n++]=ratio/6.0;
	scp*=++tscp;
	scp*=++tscp;
	pcp*=2.0*sqr_inv_qshift;
	ratio=scp*pcp;
	for(j=2;j<=HSL_SF_HZETA_EULERMACLAURIN_SERIES_ORDER;++j) {
		delta=hsl_sf_hzeta_eulermaclaurin_series_coeffs[j]*ratio;
		ans+=(terms[n++]=delta);
		scp*=++tscp;
		scp*=++tscp;
		pcp*=sqr_inv_qshift;
		ratio=scp*pcp;
		if ((fabs(delta/ans)) < (0.5*GSL_DBL_EPSILON)) break;
		}

	ans=0.0; while (n) ans+=terms[--n];
	mjr=hsl_sf_hzeta_eulermaclaurin_series_majorantratios[j]*ratio;

	*value=ans;
	*abserror=2.0*((shift+1)*GSL_DBL_EPSILON*fabs(ans)+mjr);

	return (PLFIT_SUCCESS); }

static inline
double hsl_sf_hZeta1_zed(const double s, const double q) {
#if 1
	const long double ld_q=(long double)(q);
	const long double ld_s=(long double)(s);
	const long double ld_sm1=ld_s-1.0L;
	const long double ld_logq=logl(ld_q);
	const long double ld_log1prq=log1pl(1.0L/ld_q);
	const long double ld_inv_logq=1.0L/ld_logq;
	const long double ld_logratiom1=ld_log1prq*ld_inv_logq;
	const long double ld_powratiom1=expm1l(-ld_s*ld_log1prq);
	const long double ld_varepsilon=expm1l(-ld_sm1*ld_log1prq);
	const long double ld_epsilon=ld_logratiom1+ld_powratiom1+ld_logratiom1*ld_powratiom1;
	const long double ld_z=ld_s+(ld_q+0.5L*ld_s+0.5L)*ld_epsilon+ld_q/ld_sm1*ld_inv_logq*ld_varepsilon;
	const double z=(double)(ld_z);
#else
	const double sm1=s-1.0;
	const double inv_ln_q=1.0/log(q);
	const double log1prq=log1p(1.0/q);
	const double logratiom1=log1prq*inv_ln_q;
	const double powratiom1=expm1(-s*log1prq);
	const double epsilon=logratiom1+powratiom1+logratiom1*powratiom1;
	const double z=s+(q+0.5*s+0.5)*epsilon+q/sm1*inv_ln_q*expm1(-sm1*log1prq);
#endif
	return (z); }

// Z_{1}(s,a) = -\frac{a^s}{\ln(a)} \left(\frac{1}{2}+\frac{a}{s-1}\,\left[1+\frac{1}{(s-1)\,\ln(a)}\right]\right)^{-1} \zeta^{\prime}(s,a) - 1
// Z_{1}(s,a) = O\left(\frac{(s-1)s}{6a^{2}}\right)
static
int hsl_sf_hZeta1(const double s, const double q, const double ln_q, double * value, double * abserror, double * coeff1) {
	const double criterion=ceil(10.0*s-q);
	const size_t shift=(criterion<0.0)?0:
		(criterion<HSL_SF_LNHZETA_EULERMACLAURIN_SERIES_SHIFT_MAX)?(size_t)(llrint(criterion)):
			HSL_SF_LNHZETA_EULERMACLAURIN_SERIES_SHIFT_MAX;
	const double qshift=(double)(shift)+q;
	const double ln_qshift=log(qshift);
	const double inv_qshift=1.0/qshift;
	const double inv_ln_q=1.0/ln_q;
	const double inv_ln_qshift=1.0/ln_qshift;
	const double sqr_inv_qshift=inv_qshift*inv_qshift;
	const double sm1=s-1.0;
	const double hsm1=0.5*sm1;
	const double q_over_ln_q=q*inv_ln_q;
	const double qshift_over_ln_qshift=qshift*inv_ln_qshift;
	const double qphsm1=q+hsm1;
	const double qshiftphsm1=qshift+hsm1;
	double terms[HSL_SF_LNHZETA_EULERMACLAURIN_SERIES_SHIFT_MAX+HSL_SF_HZETA_EULERMACLAURIN_SERIES_ORDER+1]={NAN};
	double delta=NAN;
	double tscp=s;
	double scp=s*sm1;
	double pcp=inv_qshift*sm1/(qshift_over_ln_qshift+sm1*qshiftphsm1);
	double lcp=1.0-inv_ln_qshift/s;
	double ratio=NAN;
	size_t n=0;
	size_t j=0;
	double ans=0.0;
	double mjr=NAN;

	if (shift) {
		const double inv_q=1.0/q;
		const double inv_sm1=1.0/sm1;
		const double w=1.0+inv_sm1*inv_ln_q;
		const double wshift=1.0+inv_sm1*inv_ln_qshift;
		const double qwphsm1=q*w+hsm1;
		const double inv_qwphsm1=1.0/qwphsm1;
		const double qshiftwshiftphsm1=qshift*wshift+hsm1;
		double qs=q;
		double a=1.0;
		for(j=0;j<shift;) {
			ans+=(terms[n++]=a*hsl_sf_hZeta1_zed(s,qs++)*inv_qwphsm1);
			a=log(qs)*inv_ln_q*exp(-s*log1p((++j)*inv_q));
			}
		pcp*=a*qshiftwshiftphsm1*inv_qwphsm1;
		}
	ratio=scp*pcp*lcp;
	ans+=terms[n++]=ratio/12.0;
	scp*=++tscp; lcp-=inv_ln_qshift/tscp;
	scp*=++tscp; lcp-=inv_ln_qshift/tscp;
	pcp*=sqr_inv_qshift;
	ratio=scp*pcp*lcp;
	for(j=2;j<=HSL_SF_HZETA_EULERMACLAURIN_SERIES_ORDER;++j) {
		delta=hsl_sf_hzeta_eulermaclaurin_series_coeffs[j]*ratio;
		ans+=(terms[n++]=delta);
		scp*=++tscp; lcp-=inv_ln_qshift/tscp;
		scp*=++tscp; lcp-=inv_ln_qshift/tscp;
		pcp*=sqr_inv_qshift;
		ratio=scp*pcp*lcp;
		if ((fabs(delta/ans)) < (0.5*GSL_DBL_EPSILON)) break;
		}

	ans=0.0; while (n) ans+=terms[--n];
	mjr=hsl_sf_hzeta_eulermaclaurin_series_majorantratios[j]*ratio;

	*value=ans;
	*abserror=2.0*((shift+1)*GSL_DBL_EPSILON*fabs(ans)+mjr);

	if (coeff1) *coeff1=1.0+q_over_ln_q/qphsm1/sm1;

	return (PLFIT_SUCCESS); }

extern
int hsl_sf_lnhzeta_deriv_tuple_e(const double s, const double q, gsl_sf_result * result, gsl_sf_result * result_deriv) {

	/* CHECK_POINTER(result) */

	if ((s <= 1.0) || (q <= 0.0)) {
		PLFIT_ERROR("s must be larger than 1.0 and q must be larger than zero", PLFIT_EINVAL);
		}
	else if (q == 1.0) {
		const double inv_sm1=1.0/(s-1.0);
		const double inv_qsm1=4.0*inv_sm1;
		const double hz_coeff0=exp2(s+1.0);
		const double hz_coeff1=1.0+inv_qsm1;
		double hZeta0_value=NAN;
		double hZeta0_abserror=NAN;
		hsl_sf_hZeta0(s,2.0,&hZeta0_value,&hZeta0_abserror);
		hZeta0_value+=1.0;
		if (result) {
			const double ln_hz_coeff=hz_coeff1/hz_coeff0;
			const double ln_hZeta0_value=ln_hz_coeff*hZeta0_value;
			result->val=log1p(ln_hZeta0_value);
			result->err=(2.0*GSL_DBL_EPSILON*ln_hz_coeff+hZeta0_abserror)/(1.0+ln_hZeta0_value);
			}

		if (result_deriv) {
			const double ld_hz_coeff2=1.0+inv_sm1*M_LOG2E;
			const double ld_hz_coeff1=1.0+inv_qsm1*ld_hz_coeff2;
			double hZeta1_value=NAN;
			double hZeta1_abserror=NAN;
			hsl_sf_hZeta1(s,2.0,M_LN2,&hZeta1_value,&hZeta1_abserror,NULL);
			hZeta0_value*=hz_coeff1;
			hZeta0_value+=hz_coeff0;
			hZeta1_value+=1.0;
			hZeta1_value*=-M_LN2*ld_hz_coeff1;
			result_deriv->val=hZeta1_value/hZeta0_value;
			result_deriv->err=2.0*GSL_DBL_EPSILON*fabs(result_deriv->val)+(hZeta0_abserror+hZeta1_abserror);
			}
		}
	else {
		const double ln_q=log(q);
		double hZeta0_value=NAN;
		double hZeta0_abserror=NAN;
		hsl_sf_hZeta0(s,q,&hZeta0_value,&hZeta0_abserror);
		if (result) {
			const double ln_hz_term0=-s*ln_q;
			const double ln_hz_term1=log(0.5+q/(s-1.0));
			result->val=ln_hz_term0+ln_hz_term1+log1p(hZeta0_value);
			result->err=2.0*GSL_DBL_EPSILON*(fabs(ln_hz_term0)+fabs(ln_hz_term1))+hZeta0_abserror/(1.0+hZeta0_value);
			}
		if (result_deriv) {
			double hZeta1_value=NAN;
			double hZeta1_abserror=NAN;
			double ld_hz_coeff1=NAN;
			hsl_sf_hZeta1(s,q,ln_q,&hZeta1_value,&hZeta1_abserror,&ld_hz_coeff1);
			result_deriv->val=-ln_q*ld_hz_coeff1*(1.0+hZeta1_value)/(1.0+hZeta0_value);
			result_deriv->err=2.0*GSL_DBL_EPSILON*fabs(result_deriv->val)+(hZeta0_abserror+hZeta1_abserror);
			}
		}

	return (PLFIT_SUCCESS); }

extern
double hsl_sf_lnhzeta_deriv_tuple(const double s, const double q, double * tuple0, double * tuple1) {
	HSL_SF_EVAL_TUPLE_RESULT(hsl_sf_lnhzeta_deriv_tuple_e(s,q,&result0,&result1)); }

extern
int hsl_sf_lnhzeta_e(const double s, const double q, gsl_sf_result * result) {
	return (hsl_sf_lnhzeta_deriv_tuple_e(s,q,result,NULL)); }

extern
double hsl_sf_lnhzeta(const double s, const double q) {
	HSL_SF_EVAL_RESULT(hsl_sf_lnhzeta_e(s,q,&result)); }

extern
int hsl_sf_lnhzeta_deriv_e(const double s, const double q, gsl_sf_result * result) {
	return (hsl_sf_lnhzeta_deriv_tuple_e(s,q,NULL,result)); }

extern
double hsl_sf_lnhzeta_deriv(const double s, const double q) {
	HSL_SF_EVAL_RESULT(hsl_sf_lnhzeta_deriv_e(s,q,&result)); }

//
// End of file `hsl/specfunc/hzeta.c'.
