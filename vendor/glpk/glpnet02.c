/* glpnet02.c (permutations to block triangular form) */

/***********************************************************************
*  This code is part of GLPK (GNU Linear Programming Kit).
*
*  This code is the result of translation of the Fortran subroutines
*  MC13D and MC13E associated with the following paper:
*
*  I.S.Duff, J.K.Reid, Algorithm 529: Permutations to block triangular
*  form, ACM Trans. on Math. Softw. 4 (1978), 189-192.
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

#include "glpnet.h"

/***********************************************************************
*  NAME
*
*  mc13d - permutations to block triangular form
*
*  SYNOPSIS
*
*  #include "glpnet.h"
*  int mc13d(int n, const int icn[], const int ip[], const int lenr[],
*     int ior[], int ib[], int lowl[], int numb[], int prev[]);
*
*  DESCRIPTION
*
*  Given the column numbers of the nonzeros in each row of the sparse
*  matrix, the routine mc13d finds a symmetric permutation that makes
*  the matrix block lower triangular.
*
*  INPUT PARAMETERS
*
*  n     order of the matrix.
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
*  OUTPUT PARAMETERS
*
*  ior   ior[i], i = 1,2,...,n, gives the position on the original
*        ordering of the row or column which is in position i in the
*        permuted form.
*
*  ib    ib[i], i = 1,2,...,num, is the row number in the permuted
*        matrix of the beginning of block i, 1 <= num <= n.
*
*  WORKING ARRAYS
*
*  arp   working array of length [1+n], where arp[0] is not used.
*        arp[i] is one less than the number of unsearched edges leaving
*        node i. At the end of the algorithm it is set to a permutation
*        which puts the matrix in block lower triangular form.
*
*  ib    working array of length [1+n], where ib[0] is not used.
*        ib[i] is the position in the ordering of the start of the ith
*        block. ib[n+1-i] holds the node number of the ith node on the
*        stack.
*
*  lowl  working array of length [1+n], where lowl[0] is not used.
*        lowl[i] is the smallest stack position of any node to which a
*        path from node i has been found. It is set to n+1 when node i
*        is removed from the stack.
*
*  numb  working array of length [1+n], where numb[0] is not used.
*        numb[i] is the position of node i in the stack if it is on it,
*        is the permuted order of node i for those nodes whose final
*        position has been found and is otherwise zero.
*
*  prev  working array of length [1+n], where prev[0] is not used.
*        prev[i] is the node at the end of the path when node i was
*        placed on the stack.
*
*  RETURNS
*
*  The routine mc13d returns num, the number of blocks found. */

int mc13d(int n, const int icn[], const int ip[], const int lenr[],
      int ior[], int ib[], int lowl[], int numb[], int prev[])
{     int *arp = ior;
      int dummy, i, i1, i2, icnt, ii, isn, ist, ist1, iv, iw, j, lcnt,
         nnm1, num, stp;
      /* icnt is the number of nodes whose positions in final ordering
         have been found. */
      icnt = 0;
      /* num is the number of blocks that have been found. */
      num = 0;
      nnm1 = n + n - 1;
      /* Initialization of arrays. */
      for (j = 1; j <= n; j++)
      {  numb[j] = 0;
         arp[j] = lenr[j] - 1;
      }
      for (isn = 1; isn <= n; isn++)
      {  /* Look for a starting node. */
         if (numb[isn] != 0) continue;
         iv = isn;
         /* ist is the number of nodes on the stack ... it is the stack
            pointer. */
         ist = 1;
         /* Put node iv at beginning of stack. */
         lowl[iv] = numb[iv] = 1;
         ib[n] = iv;
         /* The body of this loop puts a new node on the stack or
            backtracks. */
         for (dummy = 1; dummy <= nnm1; dummy++)
         {  i1 = arp[iv];
            /* Have all edges leaving node iv been searched? */
            if (i1 >= 0)
            {  i2 = ip[iv] + lenr[iv] - 1;
               i1 = i2 - i1;
               /* Look at edges leaving node iv until one enters a new
                  node or all edges are exhausted. */
               for (ii = i1; ii <= i2; ii++)
               {  iw = icn[ii];
                  /* Has node iw been on stack already? */
                  if (numb[iw] == 0) goto L70;
                  /* Update value of lowl[iv] if necessary. */
                  if (lowl[iw] < lowl[iv]) lowl[iv] = lowl[iw];
               }
               /* There are no more edges leaving node iv. */
               arp[iv] = -1;
            }
            /* Is node iv the root of a block? */
            if (lowl[iv] < numb[iv]) goto L60;
            /* Order nodes in a block. */
            num++;
            ist1 = n + 1 - ist;
            lcnt = icnt + 1;
            /* Peel block off the top of the stack starting at the top
               and working down to the root of the block. */
            for (stp = ist1; stp <= n; stp++)
            {  iw = ib[stp];
               lowl[iw] = n + 1;
               numb[iw] = ++icnt;
               if (iw == iv) break;
            }
            ist = n - stp;
            ib[num] = lcnt;
            /* Are there any nodes left on the stack? */
            if (ist != 0) goto L60;
            /* Have all the nodes been ordered? */
            if (icnt < n) break;
            goto L100;
L60:        /* Backtrack to previous node on path. */
            iw = iv;
            iv = prev[iv];
            /* Update value of lowl[iv] if necessary. */
            if (lowl[iw] < lowl[iv]) lowl[iv] = lowl[iw];
            continue;
L70:        /* Put new node on the stack. */
            arp[iv] = i2 - ii - 1;
            prev[iw] = iv;
            iv = iw;
            lowl[iv] = numb[iv] = ++ist;
            ib[n+1-ist] = iv;
         }
      }
L100: /* Put permutation in the required form. */
      for (i = 1; i <= n; i++)
         arp[numb[i]] = i;
      return num;
}

/**********************************************************************/

#if 0
#include "glplib.h"

void test(int n, int ipp);

int main(void)
{     /* test program for routine mc13d */
      test( 1,   0);
      test( 2,   1);
      test( 2,   2);
      test( 3,   3);
      test( 4,   4);
      test( 5,  10);
      test(10,  10);
      test(10,  20);
      test(20,  20);
      test(20,  50);
      test(50,  50);
      test(50, 200);
      return 0;
}

void fa01bs(int max, int *nrand);

void setup(int n, char a[1+50][1+50], int ip[], int icn[], int lenr[]);

void test(int n, int ipp)
{     int ip[1+50], icn[1+1000], ior[1+50], ib[1+51], iw[1+150],
         lenr[1+50];
      char a[1+50][1+50], hold[1+100];
      int i, ii, iblock, ij, index, j, jblock, jj, k9, num;
      xprintf("\n\n\nMatrix is of order %d and has %d off-diagonal non-"
         "zeros\n", n, ipp);
      for (j = 1; j <= n; j++)
      {  for (i = 1; i <= n; i++)
            a[i][j] = 0;
         a[j][j] = 1;
      }
      for (k9 = 1; k9 <= ipp; k9++)
      {  /* these statements should be replaced by calls to your
            favorite random number generator to place two pseudo-random
            numbers between 1 and n in the variables i and j */
         for (;;)
         {  fa01bs(n, &i);
            fa01bs(n, &j);
            if (!a[i][j]) break;
         }
         a[i][j] = 1;
      }
      /* setup converts matrix a[i,j] to required sparsity-oriented
         storage format */
      setup(n, a, ip, icn, lenr);
      num = mc13d(n, icn, ip, lenr, ior, ib, &iw[0], &iw[n], &iw[n+n]);
      /* output reordered matrix with blocking to improve clarity */
      xprintf("\nThe reordered matrix which has %d block%s is of the fo"
         "rm\n", num, num == 1 ? "" : "s");
      ib[num+1] = n + 1;
      index = 100;
      iblock = 1;
      for (i = 1; i <= n; i++)
      {  for (ij = 1; ij <= index; ij++)
            hold[ij] = ' ';
         if (i == ib[iblock])
         {  xprintf("\n");
            iblock++;
         }
         jblock = 1;
         index = 0;
         for (j = 1; j <= n; j++)
         {  if (j == ib[jblock])
            {  hold[++index] = ' ';
               jblock++;
            }
            ii = ior[i];
            jj = ior[j];
            hold[++index] = (char)(a[ii][jj] ? 'X' : '0');
         }
         xprintf("%.*s\n", index, &hold[1]);
      }
      xprintf("\nThe starting point for each block is given by\n");
      for (i = 1; i <= num; i++)
      {  if ((i - 1) % 12 == 0) xprintf("\n");
         xprintf(" %4d", ib[i]);
      }
      xprintf("\n");
      return;
}

void setup(int n, char a[1+50][1+50], int ip[], int icn[], int lenr[])
{     int i, j, ind;
      for (i = 1; i <= n; i++)
         lenr[i] = 0;
      ind = 1;
      for (i = 1; i <= n; i++)
      {  ip[i] = ind;
         for (j = 1; j <= n; j++)
         {  if (a[i][j])
            {  lenr[i]++;
               icn[ind++] = j;
            }
         }
      }
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
