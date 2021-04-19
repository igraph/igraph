/* amd.h */

/* Written by Andrew Makhorin <mao@gnu.org>. */

#ifndef GLPAMD_H
#define GLPAMD_H

#define AMD_DATE "May 31, 2007"
#define AMD_VERSION_CODE(main, sub) ((main) * 1000 + (sub))
#define AMD_MAIN_VERSION 2
#define AMD_SUB_VERSION 2
#define AMD_SUBSUB_VERSION 0
#define AMD_VERSION AMD_VERSION_CODE(AMD_MAIN_VERSION, AMD_SUB_VERSION)

#define AMD_CONTROL 5
#define AMD_INFO 20

#define AMD_DENSE 0
#define AMD_AGGRESSIVE 1

#define AMD_DEFAULT_DENSE 10.0
#define AMD_DEFAULT_AGGRESSIVE 1

#define AMD_STATUS         0
#define AMD_N              1
#define AMD_NZ             2
#define AMD_SYMMETRY       3
#define AMD_NZDIAG         4
#define AMD_NZ_A_PLUS_AT   5
#define AMD_NDENSE         6
#define AMD_MEMORY         7
#define AMD_NCMPA          8
#define AMD_LNZ            9
#define AMD_NDIV           10
#define AMD_NMULTSUBS_LDL  11
#define AMD_NMULTSUBS_LU   12
#define AMD_DMAX           13

#define AMD_OK             0
#define AMD_OUT_OF_MEMORY  (-1)
#define AMD_INVALID        (-2)
#define AMD_OK_BUT_JUMBLED 1

#define amd_order _glp_amd_order
int amd_order(int n, const int Ap[], const int Ai[], int P[],
      double Control[], double Info[]);

#define amd_2 _glp_amd_2
void amd_2(int n, int Pe[], int Iw[], int Len[], int iwlen, int pfree,
      int Nv[], int Next[], int Last[], int Head[], int Elen[],
      int Degree[], int W[], double Control[], double Info[]);

#define amd_valid _glp_amd_valid
int amd_valid(int n_row, int n_col, const int Ap[], const int Ai[]);

#define amd_defaults _glp_amd_defaults
void amd_defaults(double Control[]);

#define amd_control _glp_amd_control
void amd_control(double Control[]);

#define amd_info _glp_amd_info
void amd_info(double Info[]);

#endif

/* eof */
