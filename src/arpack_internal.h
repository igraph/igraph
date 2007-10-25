/* -*- mode: C -*-  */
/* 
   IGraph library.
   Copyright (C) 2007  Gabor Csardi <csardi@rmki.kfki.hu>
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

#ifndef ARPACK_INTERNAL_H
#define ARPACK_INTERNAL_H

/* Note: only files calling the arpack routines directly need to
   include this header.
*/

#include "config.h"

#ifndef INTERNAL_ARPACK
#define igraphdsaupd_	dsaupd_
#define igraphdseupd_	dseupd_
#define igraphdsaup2_	dsaup2_
#define igraphdstats_	dstats_
#define igraphdsesrt_	dsesrt_
#define igraphdsortr_	dsortr_
#define igraphdgetv0_	dgetv0_
#define igraphdsaitr_	dsaitr_
#define igraphdsapps_	dsapps_
#define igraphdsconv_	dsconv_
#define igraphdseigt_	dseigt_
#define igraphdsgets_	dsgets_
#define igraphdstqrb_	dstqrb_
#define igraphdmout_	dmout_
#define igraphivout_	ivout_
#define igraphsecond_	second_
#define igraphdvout_	dvout_
#endif

#ifndef INTERNAL_LAPACK
#define igraphdlarnv_	dlarnv_
#define igraphdlascl_	dlascl_
#define igraphdlartg_	dlartg_
#define igraphdlaset_	dlaset_
#define igraphdlae2_	dlae2_
#define igraphdlaev2_	dlaev2_
#define igraphdlasr_	dlasr_
#define igraphdlasrt_	dlasrt_
#define igraphdgeqr2_	dgeqr2_
#define igraphdlacpy_	dlacpy_
#define igraphdorm2r_	dorm2r_
#define igraphdsteqr_	dsteqr_
#define igraphdlanst_	dlanst_
#define igraphdlapy2_	dlapy2_
#define igraphdlamch_	dlamch_
#define igraphdlaruv_	dlaruv_
#define igraphdlarfg_	dlarfg_
#define igraphdlarf_	dlarf_
#define igraphdlae2_	dlae2_
#define igraphdlassq_	dlassq_
#define igraphdlamc2_	dlamc2_
#define igraphdlamc1_	dlamc1_
#define igraphdlamc2_	dlamc2_
#define igraphdlamc3_	dlamc3_
#define igraphdlamc4_	dlamc4_
#define igraphdlamc5_	dlamc5_
#endif

#ifndef INTERNAL_BLAS
#define igraphdaxpy_	daxpy_
#define igraphdger_	dger_
#define igraphdcopy_	dcopy_
#define igraphdscal_	dscal_
#define igraphdswap_	dswap_
#define igraphdgemv_	dgemv_
#define igraphddot_	ddot_
#define igraphdnrm2_	dnrm2_
#define igraphlsame_	lsame_
#endif

#if 0				/* internal f2c functions always used */
#define igraphd_sign	d_sign
#define igraphetime_	etime_
#define igraphpow_dd	pow_dd
#define igraphpow_di	pow_di
#define igraphs_cmp	s_cmp
#define igraphs_copy	s_copy
#endif

void igraphdsaupd_(int *ido, const char *bmat, int *n,
		   const char *which, int *nev, igraph_real_t *tol,
		   igraph_real_t *resid, int *ncv, igraph_real_t *v,
		   int *ldv, int *iparam, int *ipntr, 
		   igraph_real_t *workd, igraph_real_t *workl,
		   int *lworkl, int *info);

void igraphdseupd_(int *rvec, const char *howmny, int *select,
		   igraph_real_t *d, igraph_real_t *z, int *ldz,
		   igraph_real_t *sigma, const char *bmat, int *n,
		   const char *which, int *nev, igraph_real_t *tol,
		   igraph_real_t *resid, int *ncv, igraph_real_t *v,
		   int *ldv, int *iparam, int *ipntr, 
		   igraph_real_t *workd, igraph_real_t *workl,
		   int *lworkl, int *info);

#endif
