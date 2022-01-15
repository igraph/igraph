/* mc21a.c (permutations for zero-free diagonal) */

/***********************************************************************
*  This code is part of GLPK (GNU Linear Programming Kit).
*
*  This code is the result of translation of the Fortran subroutines
*  MC21A and MC21B associated with the following paper:
*
*  I.S.Duff, Algorithm 575: Permutations for zero-free diagonal, ACM
*  Trans. on Math. Softw. 7 (1981), 387-390.
*
*  Use of ACM Algorithms is subject to the ACM Software Copyright and
*  License Agreement. See <http://www.acm.org/publications/policies>.
*
*  The translation was made by Andrew Makhorin <mao@gnu.org>.
*
*  GLPK is free software: you can redistribute it and/or modify it
*  under the terms of the GNU General Public License as published by
*  the Free Software Foundation, either version 3 of the License, or
*  (at your option) any later version.
*
*  GLPK is distributed in the hope that it will be useful, but WITHOUT
*  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
*  or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public
*  License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with GLPK. If not, see <http://www.gnu.org/licenses/>.
***********************************************************************/

#include "mc21a.h"

/***********************************************************************
*  NAME
*
*  mc21a - permutations for zero-free diagonal
*
*  SYNOPSIS
*
*  #include "mc21a.h"
*  int mc21a(int n, const int icn[], const int ip[], const int lenr[],
*     int iperm[], int pr[], int arp[], int cv[], int out[]);
*
*  DESCRIPTION
*
*  Given the pattern of nonzeros of a sparse matrix, the routine mc21a
*  attempts to find a permutation of its rows that makes the matrix have
*  no zeros on its diagonal.
*
*  INPUT PARAMETERS
*
*  n     order of matrix.
*
*  icn   array containing the column indices of the non-zeros. Those
*        belonging to a single row must be contiguous but the ordering
*        of column indices within each row is unimportant and wasted
*        space between rows is permitted.
*
*  ip    ip[i], i = 1,2,...,n, is the position in array icn of the
*        first column index of a non-zero in row i.
*
*  lenr  lenr[i], i = 1,2,...,n, is the number of non-zeros in row i.
*
*  OUTPUT PARAMETER
*
*  iperm contains permutation to make diagonal have the smallest
*        number of zeros on it. Elements (iperm[i], i), i = 1,2,...,n,
*        are non-zero at the end of the algorithm unless the matrix is
*        structurally singular. In this case, (iperm[i], i) will be
*        zero for n - numnz entries.
*
*  WORKING ARRAYS
*
*  pr    working array of length [1+n], where pr[0] is not used.
*        pr[i] is the previous row to i in the depth first search.
*
*  arp   working array of length [1+n], where arp[0] is not used.
*        arp[i] is one less than the number of non-zeros in row i which
*        have not been scanned when looking for a cheap assignment.
*
*  cv    working array of length [1+n], where cv[0] is not used.
*        cv[i] is the most recent row extension at which column i was
*        visited.
*
*  out   working array of length [1+n], where out[0] is not used.
*        out[i] is one less than the number of non-zeros in row i
*        which have not been scanned during one pass through the main
*        loop.
*
*  RETURNS
*
*  The routine mc21a returns numnz, the number of non-zeros on diagonal
*  of permuted matrix. */

int mc21a(int n, const int icn[], const int ip[], const int lenr[],
      int iperm[], int pr[], int arp[], int cv[], int out[])
{     int i, ii, in1, in2, j, j1, jord, k, kk, numnz;
      /* Initialization of arrays. */
      for (i = 1; i <= n; i++)
      {  arp[i] = lenr[i] - 1;
         cv[i] = iperm[i] = 0;
      }
      numnz = 0;
      /* Main loop. */
      /* Each pass round this loop either results in a new assignment
       * or gives a row with no assignment. */
      for (jord = 1; jord <= n; jord++)
      {  j = jord;
         pr[j] = -1;
         for (k = 1; k <= jord; k++)
         {  /* Look for a cheap assignment. */
            in1 = arp[j];
            if (in1 >= 0)
            {  in2 = ip[j] + lenr[j] - 1;
               in1 = in2 - in1;
               for (ii = in1; ii <= in2; ii++)
               {  i = icn[ii];
                  if (iperm[i] == 0) goto L110;
               }
               /* No cheap assignment in row. */
               arp[j] = -1;
            }
            /* Begin looking for assignment chain starting with row j.*/
            out[j] = lenr[j] - 1;
            /* Inner loop. Extends chain by one or backtracks. */
            for (kk = 1; kk <= jord; kk++)
            {  in1 = out[j];
               if (in1 >= 0)
               {  in2 = ip[j] + lenr[j] - 1;
                  in1 = in2 - in1;
                  /* Forward scan. */
                  for (ii = in1; ii <= in2; ii++)
                  {  i = icn[ii];
                     if (cv[i] != jord)
                     {  /* Column i has not yet been accessed during
                         * this pass. */
                        j1 = j;
                        j = iperm[i];
                        cv[i] = jord;
                        pr[j] = j1;
                        out[j1] = in2 - ii - 1;
                        goto L100;
                     }
                  }
               }
               /* Backtracking step. */
               j = pr[j];
               if (j == -1) goto L130;
            }
L100:       ;
         }
L110:    /* New assignment is made. */
         iperm[i] = j;
         arp[j] = in2 - ii - 1;
         numnz++;
         for (k = 1; k <= jord; k++)
         {  j = pr[j];
            if (j == -1) break;
            ii = ip[j] + lenr[j] - out[j] - 2;
            i = icn[ii];
            iperm[i] = j;
         }
L130:    ;
      }
      /* If matrix is structurally singular, we now complete the
       * permutation iperm. */
      if (numnz < n)
      {  for (i = 1; i <= n; i++)
            arp[i] = 0;
         k = 0;
         for (i = 1; i <= n; i++)
         {  if (iperm[i] == 0)
               out[++k] = i;
            else
               arp[iperm[i]] = i;
         }
         k = 0;
         for (i = 1; i <= n; i++)
         {  if (arp[i] == 0)
               iperm[out[++k]] = i;
         }
      }
      return numnz;
}

/**********************************************************************/

#ifdef GLP_TEST
#include "env.h"

int sing;

void ranmat(int m, int n, int icn[], int iptr[], int nnnp1, int *knum,
      int iw[]);

void fa01bs(int max, int *nrand);

int main(void)
{     /* test program for the routine mc21a */
      /* these runs on random matrices cause all possible statements in
       * mc21a to be executed */
      int i, iold, j, j1, j2, jj, knum, l, licn, n, nov4, num, numnz;
      int ip[1+21], icn[1+1000], iperm[1+20], lenr[1+20], iw1[1+80];
      licn = 1000;
      /* run on random matrices of orders 1 through 20 */
      for (n = 1; n <= 20; n++)
      {  nov4 = n / 4;
         if (nov4 < 1) nov4 = 1;
L10:     fa01bs(nov4, &l);
         knum = l * n;
         /* knum is requested number of non-zeros in random matrix */
         if (knum > licn) goto L10;
         /* if sing is false, matrix is guaranteed structurally
          * non-singular */
         sing = ((n / 2) * 2 == n);
         /* call to subroutine to generate random matrix */
         ranmat(n, n, icn, ip, n+1, &knum, iw1);
         /* knum is now actual number of non-zeros in random matrix */
         if (knum > licn) goto L10;
         xprintf("n = %2d; nz = %4d; sing = %d\n", n, knum, sing);
         /* set up array of row lengths */
         for (i = 1; i <= n; i++)
            lenr[i] = ip[i+1] - ip[i];
         /* call to mc21a */
         numnz = mc21a(n, icn, ip, lenr, iperm, &iw1[0], &iw1[n],
            &iw1[n+n], &iw1[n+n+n]);
         /* testing to see if there are numnz non-zeros on the diagonal
          * of the permuted matrix. */
         num = 0;
         for (i = 1; i <= n; i++)
         {  iold = iperm[i];
            j1 = ip[iold];
            j2 = j1 + lenr[iold] - 1;
            if (j2 < j1) continue;
            for (jj = j1; jj <= j2; jj++)
            {  j = icn[jj];
               if (j == i)
               {  num++;
                  break;
               }
            }
         }
         if (num != numnz)
            xprintf("Failure in mc21a, numnz = %d instead of %d\n",
               numnz, num);
      }
      return 0;
}

void ranmat(int m, int n, int icn[], int iptr[], int nnnp1, int *knum,
      int iw[])
{     /* subroutine to generate random matrix */
      int i, ii, inum, j, lrow, matnum;
      inum = (*knum / n) * 2;
      if (inum > n-1) inum = n-1;
      matnum = 1;
      /* each pass through this loop generates a row of the matrix */
      for (j = 1; j <= m; j++)
      {  iptr[j] = matnum;
         if (!(sing || j > n))
            icn[matnum++] = j;
         if (n == 1) continue;
         for (i = 1; i <= n; i++) iw[i] = 0;
         if (!sing) iw[j] = 1;
         fa01bs(inum, &lrow);
         lrow--;
         if (lrow == 0) continue;
         /* lrow off-diagonal non-zeros in row j of the matrix */
         for (ii = 1; ii <= lrow; ii++)
         {  for (;;)
            {  fa01bs(n, &i);
               if (iw[i] != 1) break;
            }
            iw[i] = 1;
            icn[matnum++] = i;
         }
      }
      for (i = m+1; i <= nnnp1; i++)
         iptr[i] = matnum;
      *knum = matnum - 1;
      return;
}

double g = 1431655765.0;

double fa01as(int i)
{     /* random number generator */
      g = fmod(g * 9228907.0, 4294967296.0);
      if (i >= 0)
         return g / 4294967296.0;
      else
         return 2.0 * g / 4294967296.0 - 1.0;
}

void fa01bs(int max, int *nrand)
{     *nrand = (int)(fa01as(1) * (double)max) + 1;
      return;
}
#endif

/* eof */
