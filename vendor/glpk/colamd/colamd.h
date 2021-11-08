/* colamd.h */

/* Written by Andrew Makhorin <mao@gnu.org>. */

#ifndef COLAMD_H
#define COLAMD_H

#include "env.h"

#define COLAMD_DATE "Nov 1, 2007"
#define COLAMD_VERSION_CODE(main, sub) ((main) * 1000 + (sub))
#define COLAMD_MAIN_VERSION 2
#define COLAMD_SUB_VERSION 7
#define COLAMD_SUBSUB_VERSION 1
#define COLAMD_VERSION \
        COLAMD_VERSION_CODE(COLAMD_MAIN_VERSION, COLAMD_SUB_VERSION)

#define COLAMD_KNOBS 20
#define COLAMD_STATS 20
#define COLAMD_DENSE_ROW 0
#define COLAMD_DENSE_COL 1
#define COLAMD_AGGRESSIVE 2
#define COLAMD_DEFRAG_COUNT 2
#define COLAMD_STATUS 3
#define COLAMD_INFO1 4
#define COLAMD_INFO2 5
#define COLAMD_INFO3 6

#define COLAMD_OK                            (0)
#define COLAMD_OK_BUT_JUMBLED                (1)
#define COLAMD_ERROR_A_not_present           (-1)
#define COLAMD_ERROR_p_not_present           (-2)
#define COLAMD_ERROR_nrow_negative           (-3)
#define COLAMD_ERROR_ncol_negative           (-4)
#define COLAMD_ERROR_nnz_negative            (-5)
#define COLAMD_ERROR_p0_nonzero              (-6)
#define COLAMD_ERROR_A_too_small             (-7)
#define COLAMD_ERROR_col_length_negative     (-8)
#define COLAMD_ERROR_row_index_out_of_bounds (-9)
#define COLAMD_ERROR_out_of_memory           (-10)
#define COLAMD_ERROR_internal_error          (-999)

#define colamd_recommended _glp_colamd_recommended
size_t colamd_recommended(int nnz, int n_row, int n_col);

#define colamd_set_defaults _glp_colamd_set_defaults
void colamd_set_defaults(double knobs [COLAMD_KNOBS]);

#define colamd _glp_colamd
int colamd(int n_row, int n_col, int Alen, int A[], int p[],
      double knobs[COLAMD_KNOBS], int stats[COLAMD_STATS]);

#define symamd _glp_symamd
int symamd(int n, int A[], int p[], int perm[],
      double knobs[COLAMD_KNOBS], int stats[COLAMD_STATS],
      void *(*allocate)(size_t, size_t), void(*release)(void *));

#define colamd_report _glp_colamd_report
void colamd_report(int stats[COLAMD_STATS]);

#define symamd_report _glp_symamd_report
void symamd_report(int stats[COLAMD_STATS]);

#define colamd_printf xprintf

#endif

/* eof */
