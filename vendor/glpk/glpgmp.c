/* glpgmp.c */

/***********************************************************************
*  This code is part of GLPK (GNU Linear Programming Kit).
*
*  Copyright (C) 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008,
*  2009, 2010 Andrew Makhorin, Department for Applied Informatics,
*  Moscow Aviation Institute, Moscow, Russia. All rights reserved.
*  E-mail: <mao@gnu.org>.
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

#ifdef __clang__
#pragma clang diagnostic ignored "-Wlogical-op-parentheses"
#pragma clang diagnostic ignored "-Wsign-conversion"
#endif

#define _GLPSTD_STDIO
#include "glpdmp.h"
#include "glpgmp.h"
#define xfault xerror

#ifdef HAVE_GMP               /* use GNU MP bignum library */

int gmp_pool_count(void) { return 0; }

void gmp_free_mem(void) { return; }

#else                         /* use GLPK bignum module */

static DMP *gmp_pool = NULL;
static int gmp_size = 0;
static unsigned short *gmp_work = NULL;

void *gmp_get_atom(int size)
{     if (gmp_pool == NULL)
         gmp_pool = dmp_create_pool();
      return dmp_get_atom(gmp_pool, size);
}

void gmp_free_atom(void *ptr, int size)
{     xassert(gmp_pool != NULL);
      dmp_free_atom(gmp_pool, ptr, size);
      return;
}

int gmp_pool_count(void)
{     if (gmp_pool == NULL)
         return 0;
      else
         return dmp_in_use(gmp_pool).lo;
}

unsigned short *gmp_get_work(int size)
{     xassert(size > 0);
      if (gmp_size < size)
      {  if (gmp_size == 0)
         {  xassert(gmp_work == NULL);
            gmp_size = 100;
         }
         else
         {  xassert(gmp_work != NULL);
            xfree(gmp_work);
         }
         while (gmp_size < size) gmp_size += gmp_size;
         gmp_work = xcalloc(gmp_size, sizeof(unsigned short));
      }
      return gmp_work;
}

void gmp_free_mem(void)
{     if (gmp_pool != NULL) dmp_delete_pool(gmp_pool);
      if (gmp_work != NULL) xfree(gmp_work);
      gmp_pool = NULL;
      gmp_size = 0;
      gmp_work = NULL;
      return;
}

/*====================================================================*/

mpz_t _mpz_init(void)
{     /* initialize x, and set its value to 0 */
      mpz_t x;
      x = gmp_get_atom(sizeof(struct mpz));
      x->val = 0;
      x->ptr = NULL;
      return x;
}

void mpz_clear(mpz_t x)
{     /* free the space occupied by x */
      mpz_set_si(x, 0);
      xassert(x->ptr == NULL);
      /* free the number descriptor */
      gmp_free_atom(x, sizeof(struct mpz));
      return;
}

void mpz_set(mpz_t z, mpz_t x)
{     /* set the value of z from x */
      struct mpz_seg *e, *ee, *es;
      if (z != x)
      {  mpz_set_si(z, 0);
         z->val = x->val;
         xassert(z->ptr == NULL);
         for (e = x->ptr, es = NULL; e != NULL; e = e->next)
         {  ee = gmp_get_atom(sizeof(struct mpz_seg));
            memcpy(ee->d, e->d, 12);
            ee->next = NULL;
            if (z->ptr == NULL)
               z->ptr = ee;
            else
               es->next = ee;
            es = ee;
         }
      }
      return;
}

void mpz_set_si(mpz_t x, int val)
{     /* set the value of x to val */
      struct mpz_seg *e;
      /* free existing segments, if any */
      while (x->ptr != NULL)
      {  e = x->ptr;
         x->ptr = e->next;
         gmp_free_atom(e, sizeof(struct mpz_seg));
      }
      /* assign new value */
      if (val == 0x80000000)
      {  /* long format is needed */
         x->val = -1;
         x->ptr = e = gmp_get_atom(sizeof(struct mpz_seg));
         memset(e->d, 0, 12);
         e->d[1] = 0x8000;
         e->next = NULL;
      }
      else
      {  /* short format is enough */
         x->val = val;
      }
      return;
}

double mpz_get_d(mpz_t x)
{     /* convert x to a double, truncating if necessary */
      struct mpz_seg *e;
      int j;
      double val, deg;
      if (x->ptr == NULL)
         val = (double)x->val;
      else
      {  xassert(x->val != 0);
         val = 0.0;
         deg = 1.0;
         for (e = x->ptr; e != NULL; e = e->next)
         {  for (j = 0; j <= 5; j++)
            {  val += deg * (double)((int)e->d[j]);
               deg *= 65536.0;
            }
         }
         if (x->val < 0) val = - val;
      }
      return val;
}

double mpz_get_d_2exp(int *exp, mpz_t x)
{     /* convert x to a double, truncating if necessary (i.e. rounding
         towards zero), and returning the exponent separately;
         the return value is in the range 0.5 <= |d| < 1 and the
         exponent is stored to *exp; d*2^exp is the (truncated) x value;
         if x is zero, the return is 0.0 and 0 is stored to *exp;
         this is similar to the standard C frexp function */
      struct mpz_seg *e;
      int j, n, n1;
      double val;
      if (x->ptr == NULL)
         val = (double)x->val, n = 0;
      else
      {  xassert(x->val != 0);
         val = 0.0, n = 0;
         for (e = x->ptr; e != NULL; e = e->next)
         {  for (j = 0; j <= 5; j++)
            {  val += (double)((int)e->d[j]);
               val /= 65536.0, n += 16;
            }
         }
         if (x->val < 0) val = - val;
      }
      val = frexp(val, &n1);
      *exp = n + n1;
      return val;
}

void mpz_swap(mpz_t x, mpz_t y)
{     /* swap the values x and y efficiently */
      int val;
      void *ptr;
      val = x->val, ptr = x->ptr;
      x->val = y->val, x->ptr = y->ptr;
      y->val = val, y->ptr = ptr;
      return;
}

static void normalize(mpz_t x)
{     /* normalize integer x that includes removing non-significant
         (leading) zeros and converting to short format, if possible */
      struct mpz_seg *es, *e;
      /* if the integer is in short format, it remains unchanged */
      if (x->ptr == NULL)
      {  xassert(x->val != 0x80000000);
         goto done;
      }
      xassert(x->val == +1 || x->val == -1);
      /* find the last (most significant) non-zero segment */
      es = NULL;
      for (e = x->ptr; e != NULL; e = e->next)
      {  if (e->d[0] || e->d[1] || e->d[2] ||
             e->d[3] || e->d[4] || e->d[5]) es = e;
      }
      /* if all segments contain zeros, the integer is zero */
      if (es == NULL)
      {  mpz_set_si(x, 0);
         goto done;
      }
      /* remove non-significant (leading) zero segments */
      while (es->next != NULL)
      {  e = es->next;
         es->next = e->next;
         gmp_free_atom(e, sizeof(struct mpz_seg));
      }
      /* convert the integer to short format, if possible */
      e = x->ptr;
      if (e->next == NULL && e->d[1] <= 0x7FFF &&
         !e->d[2] && !e->d[3] && !e->d[4] && !e->d[5])
      {  int val;
         val = (int)e->d[0] + ((int)e->d[1] << 16);
         if (x->val < 0) val = - val;
         mpz_set_si(x, val);
      }
done: return;
}

void mpz_add(mpz_t z, mpz_t x, mpz_t y)
{     /* set z to x + y */
      static struct mpz_seg zero = { { 0, 0, 0, 0, 0, 0 }, NULL };
      struct mpz_seg dumx, dumy, *ex, *ey, *ez, *es, *ee;
      int k, sx, sy, sz;
      unsigned int t;
      /* if [x] = 0 then [z] = [y] */
      if (x->val == 0)
      {  xassert(x->ptr == NULL);
         mpz_set(z, y);
         goto done;
      }
      /* if [y] = 0 then [z] = [x] */
      if (y->val == 0)
      {  xassert(y->ptr == NULL);
         mpz_set(z, x);
         goto done;
      }
      /* special case when both [x] and [y] are in short format */
      if (x->ptr == NULL && y->ptr == NULL)
      {  int xval = x->val, yval = y->val, zval = x->val + y->val;
         xassert(xval != 0x80000000 && yval != 0x80000000);
         if (!(xval > 0 && yval > 0 && zval <= 0 ||
               xval < 0 && yval < 0 && zval >= 0))
         {  mpz_set_si(z, zval);
            goto done;
         }
      }
      /* convert [x] to long format, if necessary */
      if (x->ptr == NULL)
      {  xassert(x->val != 0x80000000);
         if (x->val >= 0)
         {  sx = +1;
            t = (unsigned int)(+ x->val);
         }
         else
         {  sx = -1;
            t = (unsigned int)(- x->val);
         }
         ex = &dumx;
         ex->d[0] = (unsigned short)t;
         ex->d[1] = (unsigned short)(t >> 16);
         ex->d[2] = ex->d[3] = ex->d[4] = ex->d[5] = 0;
         ex->next = NULL;
      }
      else
      {  sx = x->val;
         xassert(sx == +1 || sx == -1);
         ex = x->ptr;
      }
      /* convert [y] to long format, if necessary */
      if (y->ptr == NULL)
      {  xassert(y->val != 0x80000000);
         if (y->val >= 0)
         {  sy = +1;
            t = (unsigned int)(+ y->val);
         }
         else
         {  sy = -1;
            t = (unsigned int)(- y->val);
         }
         ey = &dumy;
         ey->d[0] = (unsigned short)t;
         ey->d[1] = (unsigned short)(t >> 16);
         ey->d[2] = ey->d[3] = ey->d[4] = ey->d[5] = 0;
         ey->next = NULL;
      }
      else
      {  sy = y->val;
         xassert(sy == +1 || sy == -1);
         ey = y->ptr;
      }
      /* main fragment */
      sz = sx;
      ez = es = NULL;
      if (sx > 0 && sy > 0 || sx < 0 && sy < 0)
      {  /* [x] and [y] have identical signs -- addition */
         t = 0;
         for (; ex || ey; ex = ex->next, ey = ey->next)
         {  if (ex == NULL) ex = &zero;
            if (ey == NULL) ey = &zero;
            ee = gmp_get_atom(sizeof(struct mpz_seg));
            for (k = 0; k <= 5; k++)
            {  t += (unsigned int)ex->d[k];
               t += (unsigned int)ey->d[k];
               ee->d[k] = (unsigned short)t;
               t >>= 16;
            }
            ee->next = NULL;
            if (ez == NULL)
               ez = ee;
            else
               es->next = ee;
            es = ee;
         }
         if (t)
         {  /* overflow -- one extra digit is needed */
            ee = gmp_get_atom(sizeof(struct mpz_seg));
            ee->d[0] = 1;
            ee->d[1] = ee->d[2] = ee->d[3] = ee->d[4] = ee->d[5] = 0;
            ee->next = NULL;
            xassert(es != NULL);
            es->next = ee;
         }
      }
      else
      {  /* [x] and [y] have different signs -- subtraction */
         t = 1;
         for (; ex || ey; ex = ex->next, ey = ey->next)
         {  if (ex == NULL) ex = &zero;
            if (ey == NULL) ey = &zero;
            ee = gmp_get_atom(sizeof(struct mpz_seg));
            for (k = 0; k <= 5; k++)
            {  t += (unsigned int)ex->d[k];
               t += (0xFFFF - (unsigned int)ey->d[k]);
               ee->d[k] = (unsigned short)t;
               t >>= 16;
            }
            ee->next = NULL;
            if (ez == NULL)
               ez = ee;
            else
               es->next = ee;
            es = ee;
         }
         if (!t)
         {  /* |[x]| < |[y]| -- result in complement coding */
            sz = - sz;
            t = 1;
            for (ee = ez; ee != NULL; ee = ee->next)
            for (k = 0; k <= 5; k++)
            {  t += (0xFFFF - (unsigned int)ee->d[k]);
               ee->d[k] = (unsigned short)t;
               t >>= 16;
            }
         }
      }
      /* contruct and normalize result */
      mpz_set_si(z, 0);
      z->val = sz;
      z->ptr = ez;
      normalize(z);
done: return;
}

void mpz_sub(mpz_t z, mpz_t x, mpz_t y)
{     /* set z to x - y */
      if (x == y)
         mpz_set_si(z, 0);
      else
      {  y->val = - y->val;
         mpz_add(z, x, y);
         if (y != z) y->val = - y->val;
      }
      return;
}

void mpz_mul(mpz_t z, mpz_t x, mpz_t y)
{     /* set z to x * y */
      struct mpz_seg dumx, dumy, *ex, *ey, *es, *e;
      int sx, sy, k, nx, ny, n;
      unsigned int t;
      unsigned short *work, *wx, *wy;
      /* if [x] = 0 then [z] = 0 */
      if (x->val == 0)
      {  xassert(x->ptr == NULL);
         mpz_set_si(z, 0);
         goto done;
      }
      /* if [y] = 0 then [z] = 0 */
      if (y->val == 0)
      {  xassert(y->ptr == NULL);
         mpz_set_si(z, 0);
         goto done;
      }
      /* special case when both [x] and [y] are in short format */
      if (x->ptr == NULL && y->ptr == NULL)
      {  int xval = x->val, yval = y->val, sz = +1;
         xassert(xval != 0x80000000 && yval != 0x80000000);
         if (xval < 0) xval = - xval, sz = - sz;
         if (yval < 0) yval = - yval, sz = - sz;
         if (xval <= 0x7FFFFFFF / yval)
         {  mpz_set_si(z, sz * (xval * yval));
            goto done;
         }
      }
      /* convert [x] to long format, if necessary */
      if (x->ptr == NULL)
      {  xassert(x->val != 0x80000000);
         if (x->val >= 0)
         {  sx = +1;
            t = (unsigned int)(+ x->val);
         }
         else
         {  sx = -1;
            t = (unsigned int)(- x->val);
         }
         ex = &dumx;
         ex->d[0] = (unsigned short)t;
         ex->d[1] = (unsigned short)(t >> 16);
         ex->d[2] = ex->d[3] = ex->d[4] = ex->d[5] = 0;
         ex->next = NULL;
      }
      else
      {  sx = x->val;
         xassert(sx == +1 || sx == -1);
         ex = x->ptr;
      }
      /* convert [y] to long format, if necessary */
      if (y->ptr == NULL)
      {  xassert(y->val != 0x80000000);
         if (y->val >= 0)
         {  sy = +1;
            t = (unsigned int)(+ y->val);
         }
         else
         {  sy = -1;
            t = (unsigned int)(- y->val);
         }
         ey = &dumy;
         ey->d[0] = (unsigned short)t;
         ey->d[1] = (unsigned short)(t >> 16);
         ey->d[2] = ey->d[3] = ey->d[4] = ey->d[5] = 0;
         ey->next = NULL;
      }
      else
      {  sy = y->val;
         xassert(sy == +1 || sy == -1);
         ey = y->ptr;
      }
      /* determine the number of digits of [x] */
      nx = n = 0;
      for (e = ex; e != NULL; e = e->next)
      for (k = 0; k <= 5; k++)
      {  n++;
         if (e->d[k]) nx = n;
      }
      xassert(nx > 0);
      /* determine the number of digits of [y] */
      ny = n = 0;
      for (e = ey; e != NULL; e = e->next)
      for (k = 0; k <= 5; k++)
      {  n++;
         if (e->d[k]) ny = n;
      }
      xassert(ny > 0);
      /* we need working array containing at least nx+ny+ny places */
      work = gmp_get_work(nx+ny+ny);
      /* load digits of [x] */
      wx = &work[0];
      for (n = 0; n < nx; n++) wx[ny+n] = 0;
      for (n = 0, e = ex; e != NULL; e = e->next)
         for (k = 0; k <= 5; k++, n++)
            if (e->d[k]) wx[ny+n] = e->d[k];
      /* load digits of [y] */
      wy = &work[nx+ny];
      for (n = 0; n < ny; n++) wy[n] = 0;
      for (n = 0, e = ey; e != NULL; e = e->next)
         for (k = 0; k <= 5; k++, n++)
            if (e->d[k]) wy[n] = e->d[k];
      /* compute [x] * [y] */
      bigmul(nx, ny, wx, wy);
      /* construct and normalize result */
      mpz_set_si(z, 0);
      z->val = sx * sy;
      es = NULL;
      k = 6;
      for (n = 0; n < nx+ny; n++)
      {  if (k > 5)
         {  e = gmp_get_atom(sizeof(struct mpz_seg));
            e->d[0] = e->d[1] = e->d[2] = 0;
            e->d[3] = e->d[4] = e->d[5] = 0;
            e->next = NULL;
            if (z->ptr == NULL)
               z->ptr = e;
            else
               es->next = e;
            es = e;
            k = 0;
         }
         es->d[k++] = wx[n];
      }
      normalize(z);
done: return;
}

void mpz_neg(mpz_t z, mpz_t x)
{     /* set z to 0 - x */
      mpz_set(z, x);
      z->val = - z->val;
      return;
}

void mpz_abs(mpz_t z, mpz_t x)
{     /* set z to the absolute value of x */
      mpz_set(z, x);
      if (z->val < 0) z->val = - z->val;
      return;
}

void mpz_div(mpz_t q, mpz_t r, mpz_t x, mpz_t y)
{     /* divide x by y, forming quotient q and/or remainder r
         if q = NULL then quotient is not stored; if r = NULL then
         remainder is not stored
         the sign of quotient is determined as in algebra while the
         sign of remainder is the same as the sign of dividend:
         +26 : +7 = +3, remainder is +5
         -26 : +7 = -3, remainder is -5
         +26 : -7 = -3, remainder is +5
         -26 : -7 = +3, remainder is -5 */
      struct mpz_seg dumx, dumy, *ex, *ey, *es, *e;
      int sx, sy, k, nx, ny, n;
      unsigned int t;
      unsigned short *work, *wx, *wy;
      /* divide by zero is not allowed */
      if (y->val == 0)
      {  xassert(y->ptr == NULL);
         xfault("mpz_div: divide by zero not allowed\n");
      }
      /* if [x] = 0 then [q] = [r] = 0 */
      if (x->val == 0)
      {  xassert(x->ptr == NULL);
         if (q != NULL) mpz_set_si(q, 0);
         if (r != NULL) mpz_set_si(r, 0);
         goto done;
      }
      /* special case when both [x] and [y] are in short format */
      if (x->ptr == NULL && y->ptr == NULL)
      {  int xval = x->val, yval = y->val;
         xassert(xval != 0x80000000 && yval != 0x80000000);
         if (q != NULL) mpz_set_si(q, xval / yval);
         if (r != NULL) mpz_set_si(r, xval % yval);
         goto done;
      }
      /* convert [x] to long format, if necessary */
      if (x->ptr == NULL)
      {  xassert(x->val != 0x80000000);
         if (x->val >= 0)
         {  sx = +1;
            t = (unsigned int)(+ x->val);
         }
         else
         {  sx = -1;
            t = (unsigned int)(- x->val);
         }
         ex = &dumx;
         ex->d[0] = (unsigned short)t;
         ex->d[1] = (unsigned short)(t >> 16);
         ex->d[2] = ex->d[3] = ex->d[4] = ex->d[5] = 0;
         ex->next = NULL;
      }
      else
      {  sx = x->val;
         xassert(sx == +1 || sx == -1);
         ex = x->ptr;
      }
      /* convert [y] to long format, if necessary */
      if (y->ptr == NULL)
      {  xassert(y->val != 0x80000000);
         if (y->val >= 0)
         {  sy = +1;
            t = (unsigned int)(+ y->val);
         }
         else
         {  sy = -1;
            t = (unsigned int)(- y->val);
         }
         ey = &dumy;
         ey->d[0] = (unsigned short)t;
         ey->d[1] = (unsigned short)(t >> 16);
         ey->d[2] = ey->d[3] = ey->d[4] = ey->d[5] = 0;
         ey->next = NULL;
      }
      else
      {  sy = y->val;
         xassert(sy == +1 || sy == -1);
         ey = y->ptr;
      }
      /* determine the number of digits of [x] */
      nx = n = 0;
      for (e = ex; e != NULL; e = e->next)
      for (k = 0; k <= 5; k++)
      {  n++;
         if (e->d[k]) nx = n;
      }
      xassert(nx > 0);
      /* determine the number of digits of [y] */
      ny = n = 0;
      for (e = ey; e != NULL; e = e->next)
      for (k = 0; k <= 5; k++)
      {  n++;
         if (e->d[k]) ny = n;
      }
      xassert(ny > 0);
      /* if nx < ny then [q] = 0 and [r] = [x] */
      if (nx < ny)
      {  if (r != NULL) mpz_set(r, x);
         if (q != NULL) mpz_set_si(q, 0);
         goto done;
      }
      /* we need working array containing at least nx+ny+1 places */
      work = gmp_get_work(nx+ny+1);
      /* load digits of [x] */
      wx = &work[0];
      for (n = 0; n < nx; n++) wx[n] = 0;
      for (n = 0, e = ex; e != NULL; e = e->next)
         for (k = 0; k <= 5; k++, n++)
            if (e->d[k]) wx[n] = e->d[k];
      /* load digits of [y] */
      wy = &work[nx+1];
      for (n = 0; n < ny; n++) wy[n] = 0;
      for (n = 0, e = ey; e != NULL; e = e->next)
         for (k = 0; k <= 5; k++, n++)
            if (e->d[k]) wy[n] = e->d[k];
      /* compute quotient and remainder */
      xassert(wy[ny-1] != 0);
      bigdiv(nx-ny, ny, wx, wy);
      /* construct and normalize quotient */
      if (q != NULL)
      {  mpz_set_si(q, 0);
         q->val = sx * sy;
         es = NULL;
         k = 6;
         for (n = ny; n <= nx; n++)
         {  if (k > 5)
            {  e = gmp_get_atom(sizeof(struct mpz_seg));
               e->d[0] = e->d[1] = e->d[2] = 0;
               e->d[3] = e->d[4] = e->d[5] = 0;
               e->next = NULL;
               if (q->ptr == NULL)
                  q->ptr = e;
               else
                  es->next = e;
               es = e;
               k = 0;
            }
            es->d[k++] = wx[n];
         }
         normalize(q);
      }
      /* construct and normalize remainder */
      if (r != NULL)
      {  mpz_set_si(r, 0);
         r->val = sx;
         es = NULL;
         k = 6;
         for (n = 0; n < ny; n++)
         {  if (k > 5)
            {  e = gmp_get_atom(sizeof(struct mpz_seg));
               e->d[0] = e->d[1] = e->d[2] = 0;
               e->d[3] = e->d[4] = e->d[5] = 0;
               e->next = NULL;
               if (r->ptr == NULL)
                  r->ptr = e;
               else
                  es->next = e;
               es = e;
               k = 0;
            }
            es->d[k++] = wx[n];
         }
         normalize(r);
      }
done: return;
}

void mpz_gcd(mpz_t z, mpz_t x, mpz_t y)
{     /* set z to the greatest common divisor of x and y */
      /* in case of arbitrary integers GCD(x, y) = GCD(|x|, |y|), and,
         in particular, GCD(0, 0) = 0 */
      mpz_t u, v, r;
      mpz_init(u);
      mpz_init(v);
      mpz_init(r);
      mpz_abs(u, x);
      mpz_abs(v, y);
      while (mpz_sgn(v))
      {  mpz_div(NULL, r, u, v);
         mpz_set(u, v);
         mpz_set(v, r);
      }
      mpz_set(z, u);
      mpz_clear(u);
      mpz_clear(v);
      mpz_clear(r);
      return;
}

int mpz_cmp(mpz_t x, mpz_t y)
{     /* compare x and y; return a positive value if x > y, zero if
         x = y, or a nefative value if x < y */
      static struct mpz_seg zero = { { 0, 0, 0, 0, 0, 0 }, NULL };
      struct mpz_seg dumx, dumy, *ex, *ey;
      int cc, sx, sy, k;
      unsigned int t;
      if (x == y)
      {  cc = 0;
         goto done;
      }
      /* special case when both [x] and [y] are in short format */
      if (x->ptr == NULL && y->ptr == NULL)
      {  int xval = x->val, yval = y->val;
         xassert(xval != 0x80000000 && yval != 0x80000000);
         cc = (xval > yval ? +1 : xval < yval ? -1 : 0);
         goto done;
      }
      /* special case when [x] and [y] have different signs */
      if (x->val > 0 && y->val <= 0 || x->val == 0 && y->val < 0)
      {  cc = +1;
         goto done;
      }
      if (x->val < 0 && y->val >= 0 || x->val == 0 && y->val > 0)
      {  cc = -1;
         goto done;
      }
      /* convert [x] to long format, if necessary */
      if (x->ptr == NULL)
      {  xassert(x->val != 0x80000000);
         if (x->val >= 0)
         {  sx = +1;
            t = (unsigned int)(+ x->val);
         }
         else
         {  sx = -1;
            t = (unsigned int)(- x->val);
         }
         ex = &dumx;
         ex->d[0] = (unsigned short)t;
         ex->d[1] = (unsigned short)(t >> 16);
         ex->d[2] = ex->d[3] = ex->d[4] = ex->d[5] = 0;
         ex->next = NULL;
      }
      else
      {  sx = x->val;
         xassert(sx == +1 || sx == -1);
         ex = x->ptr;
      }
      /* convert [y] to long format, if necessary */
      if (y->ptr == NULL)
      {  xassert(y->val != 0x80000000);
         if (y->val >= 0)
         {  sy = +1;
            t = (unsigned int)(+ y->val);
         }
         else
         {  sy = -1;
            t = (unsigned int)(- y->val);
         }
         ey = &dumy;
         ey->d[0] = (unsigned short)t;
         ey->d[1] = (unsigned short)(t >> 16);
         ey->d[2] = ey->d[3] = ey->d[4] = ey->d[5] = 0;
         ey->next = NULL;
      }
      else
      {  sy = y->val;
         xassert(sy == +1 || sy == -1);
         ey = y->ptr;
      }
      /* main fragment */
      xassert(sx > 0 && sy > 0 || sx < 0 && sy < 0);
      cc = 0;
      for (; ex || ey; ex = ex->next, ey = ey->next)
      {  if (ex == NULL) ex = &zero;
         if (ey == NULL) ey = &zero;
         for (k = 0; k <= 5; k++)
         {  if (ex->d[k] > ey->d[k]) cc = +1;
            if (ex->d[k] < ey->d[k]) cc = -1;
         }
      }
      if (sx < 0) cc = - cc;
done: return cc;
}

int mpz_sgn(mpz_t x)
{     /* return +1 if x > 0, 0 if x = 0, and -1 if x < 0 */
      int s;
      s = (x->val > 0 ? +1 : x->val < 0 ? -1 : 0);
      return s;
}

int mpz_out_str(void *_fp, int base, mpz_t x)
{     /* output x on stream fp, as a string in given base; the base
         may vary from 2 to 36;
         return the number of bytes written, or if an error occurred,
         return 0 */
      FILE *fp = _fp;
      mpz_t b, y, r;
      int n, j, nwr = 0;
      unsigned char *d;
      static char *set = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ";
      if (!(2 <= base && base <= 36))
         xfault("mpz_out_str: base = %d; invalid base\n", base);
      mpz_init(b);
      mpz_set_si(b, base);
      mpz_init(y);
      mpz_init(r);
      /* determine the number of digits */
      mpz_abs(y, x);
      for (n = 0; mpz_sgn(y) != 0; n++)
         mpz_div(y, NULL, y, b);
      if (n == 0) n = 1;
      /* compute the digits */
      d = xmalloc(n);
      mpz_abs(y, x);
      for (j = 0; j < n; j++)
      {  mpz_div(y, r, y, b);
         xassert(0 <= r->val && r->val < base && r->ptr == NULL);
         d[j] = (unsigned char)r->val;
      }
      /* output the integer to the stream */
      /* if (fp == NULL) fp = stdout; */
      if (mpz_sgn(x) < 0)
         fputc('-', fp), nwr++;
      for (j = n-1; j >= 0; j--)
         fputc(set[d[j]], fp), nwr++;
      if (ferror(fp)) nwr = 0;
      mpz_clear(b);
      mpz_clear(y);
      mpz_clear(r);
      xfree(d);
      return nwr;
}

/*====================================================================*/

mpq_t _mpq_init(void)
{     /* initialize x, and set its value to 0/1 */
      mpq_t x;
      x = gmp_get_atom(sizeof(struct mpq));
      x->p.val = 0;
      x->p.ptr = NULL;
      x->q.val = 1;
      x->q.ptr = NULL;
      return x;
}

void mpq_clear(mpq_t x)
{     /* free the space occupied by x */
      mpz_set_si(&x->p, 0);
      xassert(x->p.ptr == NULL);
      mpz_set_si(&x->q, 0);
      xassert(x->q.ptr == NULL);
      /* free the number descriptor */
      gmp_free_atom(x, sizeof(struct mpq));
      return;
}

void mpq_canonicalize(mpq_t x)
{     /* remove any factors that are common to the numerator and
         denominator of x, and make the denominator positive */
      mpz_t f;
      xassert(x->q.val != 0);
      if (x->q.val < 0)
      {  mpz_neg(&x->p, &x->p);
         mpz_neg(&x->q, &x->q);
      }
      mpz_init(f);
      mpz_gcd(f, &x->p, &x->q);
      if (!(f->val == 1 && f->ptr == NULL))
      {  mpz_div(&x->p, NULL, &x->p, f);
         mpz_div(&x->q, NULL, &x->q, f);
      }
      mpz_clear(f);
      return;
}

void mpq_set(mpq_t z, mpq_t x)
{     /* set the value of z from x */
      if (z != x)
      {  mpz_set(&z->p, &x->p);
         mpz_set(&z->q, &x->q);
      }
      return;
}

void mpq_set_si(mpq_t x, int p, unsigned int q)
{     /* set the value of x to p/q */
      if (q == 0)
         xfault("mpq_set_si: zero denominator not allowed\n");
      mpz_set_si(&x->p, p);
      xassert(q <= 0x7FFFFFFF);
      mpz_set_si(&x->q, q);
      return;
}

double mpq_get_d(mpq_t x)
{     /* convert x to a double, truncating if necessary */
      int np, nq;
      double p, q;
      p = mpz_get_d_2exp(&np, &x->p);
      q = mpz_get_d_2exp(&nq, &x->q);
      return ldexp(p / q, np - nq);
}

void mpq_set_d(mpq_t x, double val)
{     /* set x to val; there is no rounding, the conversion is exact */
      int s, n, d, j;
      double f;
      mpz_t temp;
      xassert(-DBL_MAX <= val && val <= +DBL_MAX);
      mpq_set_si(x, 0, 1);
      if (val > 0.0)
         s = +1;
      else if (val < 0.0)
         s = -1;
      else
         goto done;
      f = frexp(fabs(val), &n);
      /* |val| = f * 2^n, where 0.5 <= f < 1.0 */
      mpz_init(temp);
      while (f != 0.0)
      {  f *= 16.0, n -= 4;
         d = (int)f;
         xassert(0 <= d && d <= 15);
         f -= (double)d;
         /* x := 16 * x + d */
         mpz_set_si(temp, 16);
         mpz_mul(&x->p, &x->p, temp);
         mpz_set_si(temp, d);
         mpz_add(&x->p, &x->p, temp);
      }
      mpz_clear(temp);
      /* x := x * 2^n */
      if (n > 0)
      {  for (j = 1; j <= n; j++)
            mpz_add(&x->p, &x->p, &x->p);
      }
      else if (n < 0)
      {  for (j = 1; j <= -n; j++)
            mpz_add(&x->q, &x->q, &x->q);
         mpq_canonicalize(x);
      }
      if (s < 0) mpq_neg(x, x);
done: return;
}

void mpq_add(mpq_t z, mpq_t x, mpq_t y)
{     /* set z to x + y */
      mpz_t p, q;
      mpz_init(p);
      mpz_init(q);
      mpz_mul(p, &x->p, &y->q);
      mpz_mul(q, &x->q, &y->p);
      mpz_add(p, p, q);
      mpz_mul(q, &x->q, &y->q);
      mpz_set(&z->p, p);
      mpz_set(&z->q, q);
      mpz_clear(p);
      mpz_clear(q);
      mpq_canonicalize(z);
      return;
}

void mpq_sub(mpq_t z, mpq_t x, mpq_t y)
{     /* set z to x - y */
      mpz_t p, q;
      mpz_init(p);
      mpz_init(q);
      mpz_mul(p, &x->p, &y->q);
      mpz_mul(q, &x->q, &y->p);
      mpz_sub(p, p, q);
      mpz_mul(q, &x->q, &y->q);
      mpz_set(&z->p, p);
      mpz_set(&z->q, q);
      mpz_clear(p);
      mpz_clear(q);
      mpq_canonicalize(z);
      return;
}

void mpq_mul(mpq_t z, mpq_t x, mpq_t y)
{     /* set z to x * y */
      mpz_mul(&z->p, &x->p, &y->p);
      mpz_mul(&z->q, &x->q, &y->q);
      mpq_canonicalize(z);
      return;
}

void mpq_div(mpq_t z, mpq_t x, mpq_t y)
{     /* set z to x / y */
      mpz_t p, q;
      if (mpq_sgn(y) == 0)
         xfault("mpq_div: zero divisor not allowed\n");
      mpz_init(p);
      mpz_init(q);
      mpz_mul(p, &x->p, &y->q);
      mpz_mul(q, &x->q, &y->p);
      mpz_set(&z->p, p);
      mpz_set(&z->q, q);
      mpz_clear(p);
      mpz_clear(q);
      mpq_canonicalize(z);
      return;
}

void mpq_neg(mpq_t z, mpq_t x)
{     /* set z to 0 - x */
      mpq_set(z, x);
      mpz_neg(&z->p, &z->p);
      return;
}

void mpq_abs(mpq_t z, mpq_t x)
{     /* set z to the absolute value of x */
      mpq_set(z, x);
      mpz_abs(&z->p, &z->p);
      xassert(mpz_sgn(&x->q) > 0);
      return;
}

int mpq_cmp(mpq_t x, mpq_t y)
{     /* compare x and y; return a positive value if x > y, zero if
         x = y, or a nefative value if x < y */
      mpq_t temp;
      int s;
      mpq_init(temp);
      mpq_sub(temp, x, y);
      s = mpq_sgn(temp);
      mpq_clear(temp);
      return s;
}

int mpq_sgn(mpq_t x)
{     /* return +1 if x > 0, 0 if x = 0, and -1 if x < 0 */
      int s;
      s = mpz_sgn(&x->p);
      xassert(mpz_sgn(&x->q) > 0);
      return s;
}

int mpq_out_str(void *_fp, int base, mpq_t x)
{     /* output x on stream fp, as a string in given base; the base
         may vary from 2 to 36; output is in the form 'num/den' or if
         the denominator is 1 then just 'num';
         if the parameter fp is a null pointer, stdout is assumed;
         return the number of bytes written, or if an error occurred,
         return 0 */
      FILE *fp = _fp;
      int nwr;
      if (!(2 <= base && base <= 36))
         xfault("mpq_out_str: base = %d; invalid base\n", base);
      /* if (fp == NULL) fp = stdout; */
      nwr = mpz_out_str(fp, base, &x->p);
      if (x->q.val == 1 && x->q.ptr == NULL)
         ;
      else
      {  fputc('/', fp), nwr++;
         nwr += mpz_out_str(fp, base, &x->q);
      }
      if (ferror(fp)) nwr = 0;
      return nwr;
}

#endif

/* eof */
