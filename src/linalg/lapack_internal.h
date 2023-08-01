/* -*- mode: C -*-  */
/* vim:set ts=4 sw=4 sts=4 et: */
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

#ifndef LAPACK_INTERNAL_H
#define LAPACK_INTERNAL_H

/* Note: only files calling the LAPACK routines directly need to
   include this header.
*/

#include "igraph_decls.h"
#include "config.h"

__BEGIN_DECLS

#ifndef INTERNAL_LAPACK
    #define igraphdgeevx_   dgeevx_
    #define igraphdgeev_    dgeev_
    #define igraphdgebak_   dgebak_
    #define igraphxerbla_   xerbla_
    #define igraphdgebal_   dgebal_
    #define igraphdisnan_   disnan_
    #define igraphdlaisnan_ dlaisnan_
    #define igraphdgehrd_   dgehrd_
    #define igraphdgehd2_   dgehd2_
    #define igraphdlarf_    dlarf_
    #define igraphiladlc_   iladlc_
    #define igraphiladlr_   iladlr_
    #define igraphdlarfg_   dlarfg_
    #define igraphdlapy2_   dlapy2_
    #define igraphdlahr2_   dlahr2_
    #define igraphdlacpy_   dlacpy_
    #define igraphdlarfb_   dlarfb_
    #define igraphilaenv_   ilaenv_
    #define igraphieeeck_   ieeeck_
    #define igraphiparmq_   iparmq_
    #define igraphdhseqr_   dhseqr_
    #define igraphdlahqr_   dlahqr_
    #define igraphdlabad_   dlabad_
    #define igraphdlanv2_   dlanv2_
    #define igraphdlaqr0_   dlaqr0_
    #define igraphdlaqr3_   dlaqr3_
    #define igraphdlaqr4_   dlaqr4_
    #define igraphdlaqr2_   dlaqr2_
    #define igraphdlaset_   dlaset_
    #define igraphdormhr_   dormhr_
    #define igraphdormqr_   dormqr_
    #define igraphdlarft_   dlarft_
    #define igraphdorm2r_   dorm2r_
    #define igraphdtrexc_   dtrexc_
    #define igraphdlaexc_   dlaexc_
    #define igraphdlange_   dlange_
    #define igraphdlassq_   dlassq_
    #define igraphdlarfx_   dlarfx_
    #define igraphdlartg_   dlartg_
    #define igraphdlasy2_   dlasy2_
    #define igraphdlaqr5_   dlaqr5_
    #define igraphdlaqr1_   dlaqr1_
    #define igraphdlascl_   dlascl_
    #define igraphdorghr_   dorghr_
    #define igraphdorgqr_   dorgqr_
    #define igraphdorg2r_   dorg2r_
    #define igraphdtrevc_   dtrevc_
    #define igraphdlaln2_   dlaln2_
    #define igraphdladiv_   dladiv_
    #define igraphdsyevr_   dsyevr_
    #define igraphdsyrk_    dsyrk_
    #define igraphdlansy_   dlansy_
    #define igraphdormtr_   dormtr_
    #define igraphdormql_   dormql_
    #define igraphdorm2l_   dorm2l_
    #define igraphdstebz_   dstebz_
    #define igraphdlaebz_   dlaebz_
    #define igraphdstein_   dstein_
    #define igraphdlagtf_   dlagtf_
    #define igraphdlagts_   dlagts_
    #define igraphdlarnv_   dlarnv_
    #define igraphdlaruv_   dlaruv_
    #define igraphdstemr_   dstemr_
    #define igraphdlae2_    dlae2_
    #define igraphdlaev2_   dlaev2_
    #define igraphdlanst_   dlanst_
    #define igraphdlarrc_   dlarrc_
    #define igraphdlarre_   dlarre_
    #define igraphdlarra_   dlarra_
    #define igraphdlarrb_   dlarrb_
    #define igraphdlaneg_   dlaneg_
    #define igraphdlarrd_   dlarrd_
    #define igraphdlarrk_   dlarrk_
    #define igraphdlasq2_   dlasq2_
    #define igraphdlasq3_   dlasq3_
    #define igraphdlasq4_   dlasq4_
    #define igraphdlasq5_   dlasq5_
    #define igraphdlasq6_   dlasq6_
    #define igraphdlasrt_   dlasrt_
    #define igraphdlarrj_   dlarrj_
    #define igraphdlarrr_   dlarrr_
    #define igraphdlarrv_   dlarrv_
    #define igraphdlar1v_   dlar1v_
    #define igraphdlarrf_   dlarrf_
    #define igraphdpotrf_   dpotrf_
    #define igraphdsterf_   dsterf_
    #define igraphdsytrd_   dsytrd_
    #define igraphdlatrd_   dlatrd_
    #define igraphdsytd2_   dsytd2_
    #define igraphdlanhs_   dlanhs_
    #define igraphdgeqr2_   dgeqr2_
    #define igraphdtrsen_   dtrsen_
    #define igraphdlacn2_   dlacn2_
    #define igraphdtrsyl_   dtrsyl_
    #define igraphdlasr_    dlasr_
    #define igraphdsteqr_   dsteqr_
    #define igraphdgesv_    dgesv_
    #define igraphdgetrf_   dgetrf_
    #define igraphdgetf2_   dgetf2_
    #define igraphdlaswp_   dlaswp_
    #define igraphdgetrs_   dgetrs_
    #define igraphlen_trim_ len_trim_
    #define igraph_dlamc1_  dlamc1_
    #define igraph_dlamc2_  dlamc2_
    #define igraph_dlamc3_  dlamc3_
    #define igraph_dlamc4_  dlamc4_
    #define igraph_dlamc5_  dlamc5_
#endif

int igraphdgetrf_(int *m, int *n, double *a, int *lda, int *ipiv,
                  int *info);
int igraphdgetrs_(char *trans, int *n, int *nrhs, double *a,
                  int *lda, int *ipiv, double *b, int *ldb,
                  int *info);
int igraphdgesv_(int *n, int *nrhs, double *a, int *lda,
                 int *ipiv, double *b, int *ldb, int *info);

double igraphdlapy2_(double *x, double *y);

int igraphdsyevr_(char *jobz, char *range, char *uplo, int *n,
                  double *a, int *lda, double *vl,
                  double *vu, int * il, int *iu,
                  double *abstol, int *m, double *w,
                  double *z, int *ldz, int *isuppz,
                  double *work, int *lwork, int *iwork,
                  int *liwork, int *info);

int igraphdgeev_(char *jobvl, char *jobvr, int *n, double *a,
                 int *lda, double *wr, double *wi,
                 double *vl, int *ldvl, double *vr, int *ldvr,
                 double *work, int *lwork, int *info);

int igraphdgeevx_(char *balanc, char *jobvl, char *jobvr, char *sense,
                  int *n, double *a, int *lda, double *wr,
                  double *wi, double *vl, int *ldvl,
                  double *vr, int *ldvr, int *ilo, int *ihi,
                  double *scale, double *abnrm,
                  double *rconde, double *rcondv,
                  double *work, int *lwork, int *iwork, int *info);

int igraphdgehrd_(int *n, int *ilo, int *ihi, double *A, int *lda,
                  double *tau, double *work, int *lwork,
                  int *info);

double igraphddot_(int *n, double *dx, int *incx, double *dy, int *incy);

__END_DECLS

#endif
