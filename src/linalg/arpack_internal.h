/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2007-2012  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard street, Cambridge, MA 02139 USA

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

#include "igraph_decls.h"
#include "config.h"

__BEGIN_DECLS

#ifndef INTERNAL_ARPACK
    #define igraphdsaupd_   dsaupd_
    #define igraphdseupd_   dseupd_
    #define igraphdsaup2_   dsaup2_
    #define igraphdstats_   dstats_
    #define igraphdsesrt_   dsesrt_
    #define igraphdsortr_   dsortr_
    #define igraphdsortc_   dsortc_
    #define igraphdgetv0_   dgetv0_
    #define igraphdsaitr_   dsaitr_
    #define igraphdsapps_   dsapps_
    #define igraphdsconv_   dsconv_
    #define igraphdseigt_   dseigt_
    #define igraphdsgets_   dsgets_
    #define igraphdstqrb_   dstqrb_
    #define igraphdmout_    dmout_
    #define igraphivout_    ivout_
    #define igraphsecond_   second_
    #define igraphdvout_    dvout_
    #define igraphdnaitr_   dnaitr_
    #define igraphdnapps_   dnapps_
    #define igraphdnaup2_   dnaup2_
    #define igraphdnaupd_   dnaupd_
    #define igraphdnconv_   dnconv_
    #define igraphdlabad_   dlabad_
    #define igraphdlanhs_   dlanhs_
    #define igraphdsortc_   dsortc_
    #define igraphdneigh_   dneigh_
    #define igraphdngets_   dngets_
    #define igraphdstatn_   dstatn_
    #define igraphdlaqrb_   dlaqrb_

    #define igraphdsaupd_   dsaupd_
    #define igraphdseupd_   dseupd_
    #define igraphdnaupd_   dnaupd_
    #define igraphdneupd_   dneupd_
#endif

#ifndef INTERNAL_LAPACK
    #define igraphdlarnv_   dlarnv_
    #define igraphdlascl_   dlascl_
    #define igraphdlartg_   dlartg_
    #define igraphdlaset_   dlaset_
    #define igraphdlae2_    dlae2_
    #define igraphdlaev2_   dlaev2_
    #define igraphdlasr_    dlasr_
    #define igraphdlasrt_   dlasrt_
    #define igraphdgeqr2_   dgeqr2_
    #define igraphdlacpy_   dlacpy_
    #define igraphdorm2r_   dorm2r_
    #define igraphdsteqr_   dsteqr_
    #define igraphdlanst_   dlanst_
    #define igraphdlapy2_   dlapy2_
    #define igraphdlamch_   dlamch_
    #define igraphdlaruv_   dlaruv_
    #define igraphdlarfg_   dlarfg_
    #define igraphdlarf_    dlarf_
    #define igraphdlassq_   dlassq_
    #define igraphdlamc2_   dlamc2_
    #define igraphdlamc1_   dlamc1_
    #define igraphdlamc2_   dlamc2_
    #define igraphdlamc3_   dlamc3_
    #define igraphdlamc4_   dlamc4_
    #define igraphdlamc5_   dlamc5_
    #define igraphdlabad_   dlabad_
    #define igraphdlanhs_   dlanhs_
    #define igraphdtrevc_   dtrevc_
    #define igraphdlanv2_   dlanv2_
    #define igraphdlaln2_   dlaln2_
    #define igraphdladiv_   dladiv_
    #define igraphdtrsen_   dtrsen_
    #define igraphdlahqr_   dlahqr_
    #define igraphdtrsen_   dtrsen_
    #define igraphdlacon_   dlacon_
    #define igraphdtrsyl_   dtrsyl_
    #define igraphdtrexc_   dtrexc_
    #define igraphdlange_   dlange_
    #define igraphdlaexc_   dlaexc_
    #define igraphdlasy2_   dlasy2_
    #define igraphdlarfx_   dlarfx_
#endif

#if 0               /* internal f2c functions always used */
    #define igraphd_sign    d_sign
    #define igraphetime_    etime_
    #define igraphpow_dd    pow_dd
    #define igraphpow_di    pow_di
    #define igraphs_cmp s_cmp
    #define igraphs_copy    s_copy
    #define igraphd_lg10_   d_lg10_
    #define igraphi_dnnt_   i_dnnt_
#endif

#ifdef HAVE_GFORTRAN

/* GFortran-specific calling conventions, used when compiling the R interface.
 * Derived with "gfortran -fc-prototypes-external", applied on the original
 * Fortran sources of these functions.
 *
 * Caveats:
 *
 * 1) gfortran prints size_t for the "_len" arguments, but in fact they must be
 *    long int
 * 2) gofrtran maps Fortran LOGICAL types to int_least32_t, but in fact they
 *    must be void* (anything else doesn't work, not even _Bool*)
 * */

void igraphdsaupd_(int *ido, char *bmat, int *n,
                  char *which, int *nev, double *tol,
                  double *resid, int *ncv, double *v,
                  int *ldv, int *iparam, int *ipntr,
                  double *workd, double *workl,
                  int *lworkl, int *info,
                  long int bmat_len, long int which_len);

void igraphdseupd_(void *rvec, char *howmny, void *select,
                  double *d, double *z, int *ldz,
                  double *sigma, char *bmat, int *n,
                  char *which, int *nev, double *tol,
                  double *resid, int *ncv, double *v,
                  int *ldv, int *iparam, int *ipntr,
                  double *workd, double *workl,
                  int *lworkl, int *info,
                  long int howmny_len, long int bmat_len, long int which_len);

void igraphdnaupd_(int *ido, char *bmat, int *n,
                  char *which, int *nev, double *tol,
                  double *resid, int *ncv, double *v,
                  int *ldv, int *iparam, int *ipntr,
                  double *workd, double *workl,
                  int *lworkl, int *info,
                  long int bmat_len, long int which_len);

void igraphdneupd_(void *rvec, char *howmny, void *select,
                  double *dr, double *di,
                  double *z, int *ldz,
                  double *sigmar, double *sigmai,
                  double *workev, char *bmat, int *n,
                  char *which, int *nev, double *tol,
                  double *resid, int *ncv, double *v,
                  int *ldv, int *iparam, int *ipntr,
                  double *workd, double *workl,
                  int *lworkl, int *info,
                  long int howmny_len, long int bmat_len, long int which_len);

void igraphdsortr_(char *which, void *apply, int* n, double *x1,
                  double *x2, long int which_len);

void igraphdsortc_(char *which, void *apply, int* n, double *xreal,
                  double *ximag, double *y, long int which_len);

#else

int igraphdsaupd_(int *ido, char *bmat, int *n,
                  char *which, int *nev, double *tol,
                  double *resid, int *ncv, double *v,
                  int *ldv, int *iparam, int *ipntr,
                  double *workd, double *workl,
                  int *lworkl, int *info);

int igraphdseupd_(int *rvec, char *howmny, int *select,
                  double *d, double *z, int *ldz,
                  double *sigma, char *bmat, int *n,
                  char *which, int *nev, double *tol,
                  double *resid, int *ncv, double *v,
                  int *ldv, int *iparam, int *ipntr,
                  double *workd, double *workl,
                  int *lworkl, int *info);

int igraphdnaupd_(int *ido, char *bmat, int *n,
                  char *which, int *nev, double *tol,
                  double *resid, int *ncv, double *v,
                  int *ldv, int *iparam, int *ipntr,
                  double *workd, double *workl,
                  int *lworkl, int *info);

int igraphdneupd_(int *rvec, char *howmny, int *select,
                  double *dr, double *di,
                  double *z, int *ldz,
                  double *sigmar, double *sigmai,
                  double *workev, char *bmat, int *n,
                  char *which, int *nev, double *tol,
                  double *resid, int *ncv, double *v,
                  int *ldv, int *iparam, int *ipntr,
                  double *workd, double *workl,
                  int *lworkl, int *info);

int igraphdsortr_(char *which, int *apply, int* n, double *x1,
                  double *x2);

int igraphdsortc_(char *which, int *apply, int* n, double *xreal,
                  double *ximag, double *y);

#endif

__END_DECLS

#endif  /* ARPACK_INTERNAL_H */
