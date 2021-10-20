/* hbm.c (Harwell-Boeing sparse matrix format) */

/***********************************************************************
*  This code is part of GLPK (GNU Linear Programming Kit).
*  Copyright (C) 2004-2018 Free Software Foundation, Inc.
*  Written by Andrew Makhorin <mao@gnu.org>.
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

#include "env.h"
#include "hbm.h"
#include "misc.h"

/***********************************************************************
*  NAME
*
*  hbm_read_mat - read sparse matrix in Harwell-Boeing format
*
*  SYNOPSIS
*
*  #include "glphbm.h"
*  HBM *hbm_read_mat(const char *fname);
*
*  DESCRIPTION
*
*  The routine hbm_read_mat reads a sparse matrix in the Harwell-Boeing
*  format from a text file whose name is the character string fname.
*
*  Detailed description of the Harwell-Boeing format recognised by this
*  routine is given in the following report:
*
*  I.S.Duff, R.G.Grimes, J.G.Lewis. User's Guide for the Harwell-Boeing
*  Sparse Matrix Collection (Release I), TR/PA/92/86, October 1992.
*
*  RETURNS
*
*  If no error occured, the routine hbm_read_mat returns a pointer to
*  a data structure containing the matrix. In case of error the routine
*  prints an appropriate error message and returns NULL. */

struct dsa
{     /* working area used by routine hbm_read_mat */
      const char *fname;
      /* name of input text file */
      FILE *fp;
      /* stream assigned to input text file */
      int seqn;
      /* card sequential number */
      char card[80+1];
      /* card image buffer */
      int fmt_p;
      /* scale factor */
      int fmt_k;
      /* iterator */
      int fmt_f;
      /* format code */
      int fmt_w;
      /* field width */
      int fmt_d;
      /* number of decimal places after point */
};

/***********************************************************************
*  read_card - read next data card
*
*  This routine reads the next 80-column card from the input text file
*  and stores its image into the character string card. If the card was
*  read successfully, the routine returns zero, otherwise non-zero. */

#if 1 /* 11/III-2012 */
static int read_card(struct dsa *dsa)
{     int c, len = 0;
      char buf[255+1];
      dsa->seqn++;
      for (;;)
      {  c = fgetc(dsa->fp);
         if (c == EOF)
         {  if (ferror(dsa->fp))
               xprintf("%s:%d: read error\n",
                  dsa->fname, dsa->seqn);
            else
               xprintf("%s:%d: unexpected end-of-file\n",
                  dsa->fname, dsa->seqn);
            return 1;
         }
         else if (c == '\r')
            /* nop */;
         else if (c == '\n')
            break;
         else if (iscntrl(c))
         {  xprintf("%s:%d: invalid control character\n",
               dsa->fname, dsa->seqn, c);
            return 1;
         }
         else
         {  if (len == sizeof(buf)-1)
               goto err;
            buf[len++] = (char)c;
         }
      }
      /* remove trailing spaces */
      while (len > 80 && buf[len-1] == ' ')
         len--;
      buf[len] = '\0';
      /* line should not be longer than 80 chars */
      if (len > 80)
err:  {  xerror("%s:%d: card image too long\n",
            dsa->fname, dsa->seqn);
         return 1;
      }
      /* padd by spaces to 80-column card image */
      strcpy(dsa->card, buf);
      memset(&dsa->card[len], ' ', 80 - len);
      dsa->card[80] = '\0';
      return 0;
}
#endif

/***********************************************************************
*  scan_int - scan integer value from the current card
*
*  This routine scans an integer value from the current card, where fld
*  is the name of the field, pos is the position of the field, width is
*  the width of the field, val points to a location to which the scanned
*  value should be stored. If the value was scanned successfully, the
*  routine returns zero, otherwise non-zero. */

static int scan_int(struct dsa *dsa, char *fld, int pos, int width,
      int *val)
{     char str[80+1];
      xassert(1 <= width && width <= 80);
      memcpy(str, dsa->card + pos, width), str[width] = '\0';
      if (str2int(strspx(str), val))
      {  xprintf("%s:%d: field '%s' contains invalid value '%s'\n",
            dsa->fname, dsa->seqn, fld, str);
         return 1;
      }
      return 0;
}

/***********************************************************************
*  parse_fmt - parse Fortran format specification
*
*  This routine parses the Fortran format specification represented as
*  character string which fmt points to and stores format elements into
*  appropriate static locations. Should note that not all valid Fortran
*  format specifications may be recognised. If the format specification
*  was recognised, the routine returns zero, otherwise non-zero. */

static int parse_fmt(struct dsa *dsa, char *fmt)
{     int k, s, val;
      char str[80+1];
      /* first character should be left parenthesis */
      if (fmt[0] != '(')
fail: {  xprintf("hbm_read_mat: format '%s' not recognised\n", fmt);
         return 1;
      }
      k = 1;
      /* optional scale factor */
      dsa->fmt_p = 0;
      if (isdigit((unsigned char)fmt[k]))
      {  s = 0;
         while (isdigit((unsigned char)fmt[k]))
         {  if (s == 80) goto fail;
            str[s++] = fmt[k++];
         }
         str[s] = '\0';
         if (str2int(str, &val)) goto fail;
         if (toupper((unsigned char)fmt[k]) != 'P') goto iter;
         dsa->fmt_p = val, k++;
         if (!(0 <= dsa->fmt_p && dsa->fmt_p <= 255)) goto fail;
         /* optional comma may follow scale factor */
         if (fmt[k] == ',') k++;
      }
      /* optional iterator */
      dsa->fmt_k = 1;
      if (isdigit((unsigned char)fmt[k]))
      {  s = 0;
         while (isdigit((unsigned char)fmt[k]))
         {  if (s == 80) goto fail;
            str[s++] = fmt[k++];
         }
         str[s] = '\0';
         if (str2int(str, &val)) goto fail;
iter:    dsa->fmt_k = val;
         if (!(1 <= dsa->fmt_k && dsa->fmt_k <= 255)) goto fail;
      }
      /* format code */
      dsa->fmt_f = toupper((unsigned char)fmt[k++]);
      if (!(dsa->fmt_f == 'D' || dsa->fmt_f == 'E' ||
            dsa->fmt_f == 'F' || dsa->fmt_f == 'G' ||
            dsa->fmt_f == 'I')) goto fail;
      /* field width */
      if (!isdigit((unsigned char)fmt[k])) goto fail;
      s = 0;
      while (isdigit((unsigned char)fmt[k]))
      {  if (s == 80) goto fail;
         str[s++] = fmt[k++];
      }
      str[s] = '\0';
      if (str2int(str, &dsa->fmt_w)) goto fail;
      if (!(1 <= dsa->fmt_w && dsa->fmt_w <= 255)) goto fail;
      /* optional number of decimal places after point */
      dsa->fmt_d = 0;
      if (fmt[k] == '.')
      {  k++;
         if (!isdigit((unsigned char)fmt[k])) goto fail;
         s = 0;
         while (isdigit((unsigned char)fmt[k]))
         {  if (s == 80) goto fail;
            str[s++] = fmt[k++];
         }
         str[s] = '\0';
         if (str2int(str, &dsa->fmt_d)) goto fail;
         if (!(0 <= dsa->fmt_d && dsa->fmt_d <= 255)) goto fail;
      }
      /* last character should be right parenthesis */
      if (!(fmt[k] == ')' && fmt[k+1] == '\0')) goto fail;
      return 0;
}

/***********************************************************************
*  read_int_array - read array of integer type
*
*  This routine reads an integer array from the input text file, where
*  name is array name, fmt is Fortran format specification that controls
*  reading, n is number of array elements, val is array of integer type.
*  If the array was read successful, the routine returns zero, otherwise
*  non-zero. */

static int read_int_array(struct dsa *dsa, char *name, char *fmt,
      int n, int val[])
{     int k, pos;
      char str[80+1];
      if (parse_fmt(dsa, fmt)) return 1;
      if (!(dsa->fmt_f == 'I' && dsa->fmt_w <= 80 &&
            dsa->fmt_k * dsa->fmt_w <= 80))
      {  xprintf(
            "%s:%d: can't read array '%s' - invalid format '%s'\n",
            dsa->fname, dsa->seqn, name, fmt);
         return 1;
      }
      for (k = 1, pos = INT_MAX; k <= n; k++, pos++)
      {  if (pos >= dsa->fmt_k)
         {  if (read_card(dsa)) return 1;
            pos = 0;
         }
         memcpy(str, dsa->card + dsa->fmt_w * pos, dsa->fmt_w);
         str[dsa->fmt_w] = '\0';
         strspx(str);
         if (str2int(str, &val[k]))
         {  xprintf(
               "%s:%d: can't read array '%s' - invalid value '%s'\n",
               dsa->fname, dsa->seqn, name, str);
            return 1;
         }
      }
      return 0;
}

/***********************************************************************
*  read_real_array - read array of real type
*
*  This routine reads a real array from the input text file, where name
*  is array name, fmt is Fortran format specification that controls
*  reading, n is number of array elements, val is array of real type.
*  If the array was read successful, the routine returns zero, otherwise
*  non-zero. */

static int read_real_array(struct dsa *dsa, char *name, char *fmt,
      int n, double val[])
{     int k, pos;
      char str[80+1], *ptr;
      if (parse_fmt(dsa, fmt)) return 1;
      if (!(dsa->fmt_f != 'I' && dsa->fmt_w <= 80 &&
            dsa->fmt_k * dsa->fmt_w <= 80))
      {  xprintf(
            "%s:%d: can't read array '%s' - invalid format '%s'\n",
            dsa->fname, dsa->seqn, name, fmt);
         return 1;
      }
      for (k = 1, pos = INT_MAX; k <= n; k++, pos++)
      {  if (pos >= dsa->fmt_k)
         {  if (read_card(dsa)) return 1;
            pos = 0;
         }
         memcpy(str, dsa->card + dsa->fmt_w * pos, dsa->fmt_w);
         str[dsa->fmt_w] = '\0';
         strspx(str);
         if (strchr(str, '.') == NULL && strcmp(str, "0"))
         {  xprintf("%s(%d): can't read array '%s' - value '%s' has no "
               "decimal point\n", dsa->fname, dsa->seqn, name, str);
            return 1;
         }
         /* sometimes lower case letters appear */
         for (ptr = str; *ptr; ptr++)
            *ptr = (char)toupper((unsigned char)*ptr);
         ptr = strchr(str, 'D');
         if (ptr != NULL) *ptr = 'E';
         /* value may appear with decimal exponent but without letters
            E or D (for example, -123.456-012), so missing letter should
            be inserted */
         ptr = strchr(str+1, '+');
         if (ptr == NULL) ptr = strchr(str+1, '-');
         if (ptr != NULL && *(ptr-1) != 'E')
         {  xassert(strlen(str) < 80);
            memmove(ptr+1, ptr, strlen(ptr)+1);
            *ptr = 'E';
         }
         if (str2num(str, &val[k]))
         {  xprintf(
               "%s:%d: can't read array '%s' - invalid value '%s'\n",
               dsa->fname, dsa->seqn, name, str);
            return 1;
         }
      }
      return 0;
}

HBM *hbm_read_mat(const char *fname)
{     struct dsa _dsa, *dsa = &_dsa;
      HBM *hbm = NULL;
      dsa->fname = fname;
      xprintf("hbm_read_mat: reading matrix from '%s'...\n",
         dsa->fname);
      dsa->fp = fopen(dsa->fname, "r");
      if (dsa->fp == NULL)
      {  xprintf("hbm_read_mat: unable to open '%s' - %s\n",
#if 0 /* 29/I-2017 */
            dsa->fname, strerror(errno));
#else
            dsa->fname, xstrerr(errno));
#endif
         goto fail;
      }
      dsa->seqn = 0;
      hbm = xmalloc(sizeof(HBM));
      memset(hbm, 0, sizeof(HBM));
      /* read the first heading card */
      if (read_card(dsa)) goto fail;
      memcpy(hbm->title, dsa->card, 72), hbm->title[72] = '\0';
      strtrim(hbm->title);
      xprintf("%s\n", hbm->title);
      memcpy(hbm->key, dsa->card+72, 8), hbm->key[8] = '\0';
      strspx(hbm->key);
      xprintf("key = %s\n", hbm->key);
      /* read the second heading card */
      if (read_card(dsa)) goto fail;
      if (scan_int(dsa, "totcrd",  0, 14, &hbm->totcrd)) goto fail;
      if (scan_int(dsa, "ptrcrd", 14, 14, &hbm->ptrcrd)) goto fail;
      if (scan_int(dsa, "indcrd", 28, 14, &hbm->indcrd)) goto fail;
      if (scan_int(dsa, "valcrd", 42, 14, &hbm->valcrd)) goto fail;
      if (scan_int(dsa, "rhscrd", 56, 14, &hbm->rhscrd)) goto fail;
      xprintf("totcrd = %d; ptrcrd = %d; indcrd = %d; valcrd = %d; rhsc"
         "rd = %d\n", hbm->totcrd, hbm->ptrcrd, hbm->indcrd,
         hbm->valcrd, hbm->rhscrd);
      /* read the third heading card */
      if (read_card(dsa)) goto fail;
      memcpy(hbm->mxtype, dsa->card, 3), hbm->mxtype[3] = '\0';
      if (strchr("RCP",   hbm->mxtype[0]) == NULL ||
          strchr("SUHZR", hbm->mxtype[1]) == NULL ||
          strchr("AE",    hbm->mxtype[2]) == NULL)
      {  xprintf("%s:%d: matrix type '%s' not recognised\n",
            dsa->fname, dsa->seqn, hbm->mxtype);
         goto fail;
      }
      if (scan_int(dsa, "nrow", 14, 14, &hbm->nrow)) goto fail;
      if (scan_int(dsa, "ncol", 28, 14, &hbm->ncol)) goto fail;
      if (scan_int(dsa, "nnzero", 42, 14, &hbm->nnzero)) goto fail;
      if (scan_int(dsa, "neltvl", 56, 14, &hbm->neltvl)) goto fail;
      xprintf("mxtype = %s; nrow = %d; ncol = %d; nnzero = %d; neltvl ="
         " %d\n", hbm->mxtype, hbm->nrow, hbm->ncol, hbm->nnzero,
         hbm->neltvl);
      /* read the fourth heading card */
      if (read_card(dsa)) goto fail;
      memcpy(hbm->ptrfmt, dsa->card, 16), hbm->ptrfmt[16] = '\0';
      strspx(hbm->ptrfmt);
      memcpy(hbm->indfmt, dsa->card+16, 16), hbm->indfmt[16] = '\0';
      strspx(hbm->indfmt);
      memcpy(hbm->valfmt, dsa->card+32, 20), hbm->valfmt[20] = '\0';
      strspx(hbm->valfmt);
      memcpy(hbm->rhsfmt, dsa->card+52, 20), hbm->rhsfmt[20] = '\0';
      strspx(hbm->rhsfmt);
      xprintf("ptrfmt = %s; indfmt = %s; valfmt = %s; rhsfmt = %s\n",
         hbm->ptrfmt, hbm->indfmt, hbm->valfmt, hbm->rhsfmt);
      /* read the fifth heading card (optional) */
      if (hbm->rhscrd <= 0)
      {  strcpy(hbm->rhstyp, "???");
         hbm->nrhs = 0;
         hbm->nrhsix = 0;
      }
      else
      {  if (read_card(dsa)) goto fail;
         memcpy(hbm->rhstyp, dsa->card, 3), hbm->rhstyp[3] = '\0';
         if (scan_int(dsa, "nrhs", 14, 14, &hbm->nrhs)) goto fail;
         if (scan_int(dsa, "nrhsix", 28, 14, &hbm->nrhsix)) goto fail;
         xprintf("rhstyp = '%s'; nrhs = %d; nrhsix = %d\n",
            hbm->rhstyp, hbm->nrhs, hbm->nrhsix);
      }
      /* read matrix structure */
      hbm->colptr = xcalloc(1+hbm->ncol+1, sizeof(int));
      if (read_int_array(dsa, "colptr", hbm->ptrfmt, hbm->ncol+1,
         hbm->colptr)) goto fail;
      hbm->rowind = xcalloc(1+hbm->nnzero, sizeof(int));
      if (read_int_array(dsa, "rowind", hbm->indfmt, hbm->nnzero,
         hbm->rowind)) goto fail;
      /* read matrix values */
      if (hbm->valcrd <= 0) goto done;
      if (hbm->mxtype[2] == 'A')
      {  /* assembled matrix */
         hbm->values = xcalloc(1+hbm->nnzero, sizeof(double));
         if (read_real_array(dsa, "values", hbm->valfmt, hbm->nnzero,
            hbm->values)) goto fail;
      }
      else
      {  /* elemental (unassembled) matrix */
         hbm->values = xcalloc(1+hbm->neltvl, sizeof(double));
         if (read_real_array(dsa, "values", hbm->valfmt, hbm->neltvl,
            hbm->values)) goto fail;
      }
      /* read right-hand sides */
      if (hbm->nrhs <= 0) goto done;
      if (hbm->rhstyp[0] == 'F')
      {  /* dense format */
         hbm->nrhsvl = hbm->nrow * hbm->nrhs;
         hbm->rhsval = xcalloc(1+hbm->nrhsvl, sizeof(double));
         if (read_real_array(dsa, "rhsval", hbm->rhsfmt, hbm->nrhsvl,
            hbm->rhsval)) goto fail;
      }
      else if (hbm->rhstyp[0] == 'M' && hbm->mxtype[2] == 'A')
      {  /* sparse format */
         /* read pointers */
         hbm->rhsptr = xcalloc(1+hbm->nrhs+1, sizeof(int));
         if (read_int_array(dsa, "rhsptr", hbm->ptrfmt, hbm->nrhs+1,
            hbm->rhsptr)) goto fail;
         /* read sparsity pattern */
         hbm->rhsind = xcalloc(1+hbm->nrhsix, sizeof(int));
         if (read_int_array(dsa, "rhsind", hbm->indfmt, hbm->nrhsix,
            hbm->rhsind)) goto fail;
         /* read values */
         hbm->rhsval = xcalloc(1+hbm->nrhsix, sizeof(double));
         if (read_real_array(dsa, "rhsval", hbm->rhsfmt, hbm->nrhsix,
            hbm->rhsval)) goto fail;
      }
      else if (hbm->rhstyp[0] == 'M' && hbm->mxtype[2] == 'E')
      {  /* elemental format */
         hbm->rhsval = xcalloc(1+hbm->nrhsvl, sizeof(double));
         if (read_real_array(dsa, "rhsval", hbm->rhsfmt, hbm->nrhsvl,
            hbm->rhsval)) goto fail;
      }
      else
      {  xprintf("%s:%d: right-hand side type '%c' not recognised\n",
            dsa->fname, dsa->seqn, hbm->rhstyp[0]);
         goto fail;
      }
      /* read starting guesses */
      if (hbm->rhstyp[1] == 'G')
      {  hbm->nguess = hbm->nrow * hbm->nrhs;
         hbm->sguess = xcalloc(1+hbm->nguess, sizeof(double));
         if (read_real_array(dsa, "sguess", hbm->rhsfmt, hbm->nguess,
            hbm->sguess)) goto fail;
      }
      /* read solution vectors */
      if (hbm->rhstyp[2] == 'X')
      {  hbm->nexact = hbm->nrow * hbm->nrhs;
         hbm->xexact = xcalloc(1+hbm->nexact, sizeof(double));
         if (read_real_array(dsa, "xexact", hbm->rhsfmt, hbm->nexact,
            hbm->xexact)) goto fail;
      }
done: /* reading has been completed */
      xprintf("hbm_read_mat: %d cards were read\n", dsa->seqn);
      fclose(dsa->fp);
      return hbm;
fail: /* something wrong in Danish kingdom */
      if (hbm != NULL)
      {  if (hbm->colptr != NULL) xfree(hbm->colptr);
         if (hbm->rowind != NULL) xfree(hbm->rowind);
         if (hbm->rhsptr != NULL) xfree(hbm->rhsptr);
         if (hbm->rhsind != NULL) xfree(hbm->rhsind);
         if (hbm->values != NULL) xfree(hbm->values);
         if (hbm->rhsval != NULL) xfree(hbm->rhsval);
         if (hbm->sguess != NULL) xfree(hbm->sguess);
         if (hbm->xexact != NULL) xfree(hbm->xexact);
         xfree(hbm);
      }
      if (dsa->fp != NULL) fclose(dsa->fp);
      return NULL;
}

/***********************************************************************
*  NAME
*
*  hbm_free_mat - free sparse matrix in Harwell-Boeing format
*
*  SYNOPSIS
*
*  #include "glphbm.h"
*  void hbm_free_mat(HBM *hbm);
*
*  DESCRIPTION
*
*  The hbm_free_mat routine frees all the memory allocated to the data
*  structure containing a sparse matrix in the Harwell-Boeing format. */

void hbm_free_mat(HBM *hbm)
{     if (hbm->colptr != NULL) xfree(hbm->colptr);
      if (hbm->rowind != NULL) xfree(hbm->rowind);
      if (hbm->rhsptr != NULL) xfree(hbm->rhsptr);
      if (hbm->rhsind != NULL) xfree(hbm->rhsind);
      if (hbm->values != NULL) xfree(hbm->values);
      if (hbm->rhsval != NULL) xfree(hbm->rhsval);
      if (hbm->sguess != NULL) xfree(hbm->sguess);
      if (hbm->xexact != NULL) xfree(hbm->xexact);
      xfree(hbm);
      return;
}

/* eof */
