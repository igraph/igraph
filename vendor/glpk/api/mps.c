/* mps.c (MPS format routines) */

/***********************************************************************
*  This code is part of GLPK (GNU Linear Programming Kit).
*  Copyright (C) 2008-2016 Free Software Foundation, Inc.
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
#include "misc.h"
#include "prob.h"

#define xfprintf glp_format

/***********************************************************************
*  NAME
*
*  glp_init_mpscp - initialize MPS format control parameters
*
*  SYNOPSIS
*
*  void glp_init_mpscp(glp_mpscp *parm);
*
*  DESCRIPTION
*
*  The routine glp_init_mpscp initializes control parameters, which are
*  used by the MPS input/output routines glp_read_mps and glp_write_mps,
*  with default values.
*
*  Default values of the control parameters are stored in the glp_mpscp
*  structure, which the parameter parm points to. */

void glp_init_mpscp(glp_mpscp *parm)
{     parm->blank = '\0';
      parm->obj_name = NULL;
      parm->tol_mps = 1e-12;
      return;
}

static void check_parm(const char *func, const glp_mpscp *parm)
{     /* check control parameters */
      if (!(0x00 <= parm->blank && parm->blank <= 0xFF) ||
          !(parm->blank == '\0' || isprint(parm->blank)))
         xerror("%s: blank = 0x%02X; invalid parameter\n",
            func, parm->blank);
      if (!(parm->obj_name == NULL || strlen(parm->obj_name) <= 255))
         xerror("%s: obj_name = \"%.12s...\"; parameter too long\n",
            func, parm->obj_name);
      if (!(0.0 <= parm->tol_mps && parm->tol_mps < 1.0))
         xerror("%s: tol_mps = %g; invalid parameter\n",
            func, parm->tol_mps);
      return;
}

/***********************************************************************
*  NAME
*
*  glp_read_mps - read problem data in MPS format
*
*  SYNOPSIS
*
*  int glp_read_mps(glp_prob *P, int fmt, const glp_mpscp *parm,
*     const char *fname);
*
*  DESCRIPTION
*
*  The routine glp_read_mps reads problem data in MPS format from a
*  text file.
*
*  The parameter fmt specifies the version of MPS format:
*
*  GLP_MPS_DECK - fixed (ancient) MPS format;
*  GLP_MPS_FILE - free (modern) MPS format.
*
*  The parameter parm is a pointer to the structure glp_mpscp, which
*  specifies control parameters used by the routine. If parm is NULL,
*  the routine uses default settings.
*
*  The character string fname specifies a name of the text file to be
*  read.
*
*  Note that before reading data the current content of the problem
*  object is completely erased with the routine glp_erase_prob.
*
*  RETURNS
*
*  If the operation was successful, the routine glp_read_mps returns
*  zero. Otherwise, it prints an error message and returns non-zero. */

struct csa
{     /* common storage area */
      glp_prob *P;
      /* pointer to problem object */
      int deck;
      /* MPS format (0 - free, 1 - fixed) */
      const glp_mpscp *parm;
      /* pointer to control parameters */
      const char *fname;
      /* name of input MPS file */
      glp_file *fp;
      /* stream assigned to input MPS file */
      jmp_buf jump;
      /* label for go to in case of error */
      int recno;
      /* current record (card) number */
      int recpos;
      /* current record (card) position */
      int c;
      /* current character */
      int fldno;
      /* current field number */
      char field[255+1];
      /* current field content */
      int w80;
      /* warning 'record must not be longer than 80 chars' issued */
      int wef;
      /* warning 'extra fields detected beyond field 6' issued */
      int obj_row;
      /* objective row number */
      void *work1, *work2, *work3;
      /* working arrays */
};

static void error(struct csa *csa, const char *fmt, ...)
{     /* print error message and terminate processing */
      va_list arg;
      xprintf("%s:%d: ", csa->fname, csa->recno);
      va_start(arg, fmt);
      xvprintf(fmt, arg);
      va_end(arg);
      longjmp(csa->jump, 1);
      /* no return */
}

static void warning(struct csa *csa, const char *fmt, ...)
{     /* print warning message and continue processing */
      va_list arg;
      xprintf("%s:%d: warning: ", csa->fname, csa->recno);
      va_start(arg, fmt);
      xvprintf(fmt, arg);
      va_end(arg);
      return;
}

static void read_char(struct csa *csa)
{     /* read next character */
      int c;
      if (csa->c == '\n')
         csa->recno++, csa->recpos = 0;
      csa->recpos++;
read: c = glp_getc(csa->fp);
      if (c < 0)
      {  if (glp_ioerr(csa->fp))
            error(csa, "read error - %s\n", get_err_msg());
         else if (csa->c == '\n')
            error(csa, "unexpected end of file\n");
         else
         {  warning(csa, "missing final end of line\n");
            c = '\n';
         }
      }
      else if (c == '\n')
         ;
      else if (csa->c == '\r')
      {  c = '\r';
         goto badc;
      }
      else if (csa->deck && c == '\r')
      {  csa->c = '\r';
         goto read;
      }
      else if (c == ' ')
         ;
      else if (isspace(c))
      {  if (csa->deck)
badc:       error(csa, "in fixed MPS format white-space character 0x%02"
               "X is not allowed\n", c);
         c = ' ';
      }
      else if (iscntrl(c))
         error(csa, "invalid control character 0x%02X\n", c);
      if (csa->deck && csa->recpos == 81 && c != '\n' && csa->w80 < 1)
      {  warning(csa, "in fixed MPS format record must not be longer th"
            "an 80 characters\n");
         csa->w80++;
      }
      csa->c = c;
      return;
}

static int indicator(struct csa *csa, int name)
{     /* skip comment records and read possible indicator record */
      int ret;
      /* reset current field number */
      csa->fldno = 0;
loop: /* read the very first character of the next record */
      xassert(csa->c == '\n');
      read_char(csa);
      if (csa->c == ' ' || csa->c == '\n')
      {  /* data record */
         ret = 0;
      }
      else if (csa->c == '*')
      {  /* comment record */
         while (csa->c != '\n')
            read_char(csa);
         goto loop;
      }
      else
      {  /* indicator record */
         int len = 0;
         while (csa->c != ' ' && csa->c != '\n' && len < 12)
         {  csa->field[len++] = (char)csa->c;
            read_char(csa);
         }
         csa->field[len] = '\0';
         if (!(strcmp(csa->field, "NAME")    == 0 ||
               strcmp(csa->field, "ROWS")    == 0 ||
               strcmp(csa->field, "COLUMNS") == 0 ||
               strcmp(csa->field, "RHS")     == 0 ||
               strcmp(csa->field, "RANGES")  == 0 ||
               strcmp(csa->field, "BOUNDS")  == 0 ||
               strcmp(csa->field, "ENDATA")  == 0))
            error(csa, "invalid indicator record\n");
         if (!name)
         {  while (csa->c != '\n')
               read_char(csa);
         }
         ret = 1;
      }
      return ret;
}

static void read_field(struct csa *csa)
{     /* read next field of the current data record */
      csa->fldno++;
      if (csa->deck)
      {  /* fixed MPS format */
         int beg, end, pos;
         /* determine predefined field positions */
         if (csa->fldno == 1)
            beg = 2, end = 3;
         else if (csa->fldno == 2)
            beg = 5, end = 12;
         else if (csa->fldno == 3)
            beg = 15, end = 22;
         else if (csa->fldno == 4)
            beg = 25, end = 36;
         else if (csa->fldno == 5)
            beg = 40, end = 47;
         else if (csa->fldno == 6)
            beg = 50, end = 61;
         else
            xassert(csa != csa);
         /* skip blanks preceding the current field */
         if (csa->c != '\n')
         {  pos = csa->recpos;
            while (csa->recpos < beg)
            {  if (csa->c == ' ')
                  ;
               else if (csa->c == '\n')
                  break;
               else
                  error(csa, "in fixed MPS format positions %d-%d must "
                     "be blank\n", pos, beg-1);
               read_char(csa);
            }
         }
         /* skip possible comment beginning in the field 3 or 5 */
         if ((csa->fldno == 3 || csa->fldno == 5) && csa->c == '$')
         {  while (csa->c != '\n')
               read_char(csa);
         }
         /* read the current field */
         for (pos = beg; pos <= end; pos++)
         {  if (csa->c == '\n') break;
            csa->field[pos-beg] = (char)csa->c;
            read_char(csa);
         }
         csa->field[pos-beg] = '\0';
         strtrim(csa->field);
         /* skip blanks following the last field */
         if (csa->fldno == 6 && csa->c != '\n')
         {  while (csa->recpos <= 72)
            {  if (csa->c == ' ')
                  ;
               else if (csa->c == '\n')
                  break;
               else
                  error(csa, "in fixed MPS format positions 62-72 must "
                     "be blank\n");
               read_char(csa);
            }
            while (csa->c != '\n')
               read_char(csa);
         }
      }
      else
      {  /* free MPS format */
         int len;
         /* skip blanks preceding the current field */
         while (csa->c == ' ')
            read_char(csa);
         /* skip possible comment */
         if (csa->c == '$')
         {  while (csa->c != '\n')
               read_char(csa);
         }
         /* read the current field */
         len = 0;
         while (!(csa->c == ' ' || csa->c == '\n'))
         {  if (len == 255)
               error(csa, "length of field %d exceeds 255 characters\n",
                  csa->fldno++);
            csa->field[len++] = (char)csa->c;
            read_char(csa);
         }
         csa->field[len] = '\0';
         /* skip anything following the last field (any extra fields
            are considered to be comments) */
         if (csa->fldno == 6)
         {  while (csa->c == ' ')
               read_char(csa);
            if (csa->c != '$' && csa->c != '\n' && csa->wef < 1)
            {  warning(csa, "some extra field(s) detected beyond field "
                  "6; field(s) ignored\n");
               csa->wef++;
            }
            while (csa->c != '\n')
               read_char(csa);
         }
      }
      return;
}

static void patch_name(struct csa *csa, char *name)
{     /* process embedded blanks in symbolic name */
      int blank = csa->parm->blank;
      if (blank == '\0')
      {  /* remove emedded blanks */
         strspx(name);
      }
      else
      {  /* replace embedded blanks by specified character */
         for (; *name != '\0'; name++)
            if (*name == ' ') *name = (char)blank;
      }
      return;
}

static double read_number(struct csa *csa)
{     /* read next field and convert it to floating-point number */
      double x;
      char *s;
      /* read next field */
      read_field(csa);
      xassert(csa->fldno == 4 || csa->fldno == 6);
      if (csa->field[0] == '\0')
         error(csa, "missing numeric value in field %d\n", csa->fldno);
      /* skip initial spaces of the field */
      for (s = csa->field; *s == ' '; s++);
      /* perform conversion */
      if (str2num(s, &x) != 0)
         error(csa, "cannot convert '%s' to floating-point number\n",
            s);
      return x;
}

static void skip_field(struct csa *csa)
{     /* read and skip next field (assumed to be blank) */
      read_field(csa);
      if (csa->field[0] != '\0')
         error(csa, "field %d must be blank\n", csa->fldno);
      return;
}

static void read_name(struct csa *csa)
{     /* read NAME indicator record */
      if (!(indicator(csa, 1) && strcmp(csa->field, "NAME") == 0))
         error(csa, "missing NAME indicator record\n");
      /* this indicator record looks like a data record; simulate that
         fields 1 and 2 were read */
      csa->fldno = 2;
      /* field 3: model name */
      read_field(csa), patch_name(csa, csa->field);
      if (csa->field[0] == '\0')
         warning(csa, "missing model name in field 3\n");
      else
         glp_set_prob_name(csa->P, csa->field);
      /* skip anything following field 3 */
      while (csa->c != '\n')
         read_char(csa);
      return;
}

static void read_rows(struct csa *csa)
{     /* read ROWS section */
      int i, type;
loop: if (indicator(csa, 0)) goto done;
      /* field 1: row type */
      read_field(csa), strspx(csa->field);
      if (strcmp(csa->field, "N") == 0)
         type = GLP_FR;
      else if (strcmp(csa->field, "G") == 0)
         type = GLP_LO;
      else if (strcmp(csa->field, "L") == 0)
         type = GLP_UP;
      else if (strcmp(csa->field, "E") == 0)
         type = GLP_FX;
      else if (csa->field[0] == '\0')
         error(csa, "missing row type in field 1\n");
      else
         error(csa, "invalid row type in field 1\n");
      /* field 2: row name */
      read_field(csa), patch_name(csa, csa->field);
      if (csa->field[0] == '\0')
         error(csa, "missing row name in field 2\n");
      if (glp_find_row(csa->P, csa->field) != 0)
         error(csa, "row '%s' multiply specified\n", csa->field);
      i = glp_add_rows(csa->P, 1);
      glp_set_row_name(csa->P, i, csa->field);
      glp_set_row_bnds(csa->P, i, type, 0.0, 0.0);
      /* fields 3, 4, 5, and 6 must be blank */
      skip_field(csa);
      skip_field(csa);
      skip_field(csa);
      skip_field(csa);
      goto loop;
done: return;
}

static void read_columns(struct csa *csa)
{     /* read COLUMNS section */
      int i, j, f, len, kind = GLP_CV, *ind;
      double aij, *val;
      char name[255+1], *flag;
      /* allocate working arrays */
      csa->work1 = ind = xcalloc(1+csa->P->m, sizeof(int));
      csa->work2 = val = xcalloc(1+csa->P->m, sizeof(double));
      csa->work3 = flag = xcalloc(1+csa->P->m, sizeof(char));
      memset(&flag[1], 0, csa->P->m);
      /* no current column exists */
      j = 0, len = 0;
loop: if (indicator(csa, 0)) goto done;
      /* field 1 must be blank */
      if (csa->deck)
      {  read_field(csa);
         if (csa->field[0] != '\0')
            error(csa, "field 1 must be blank\n");
      }
      else
         csa->fldno++;
      /* field 2: column or kind name */
      read_field(csa), patch_name(csa, csa->field);
      strcpy(name, csa->field);
      /* field 3: row name or keyword 'MARKER' */
      read_field(csa), patch_name(csa, csa->field);
      if (strcmp(csa->field, "'MARKER'") == 0)
      {  /* process kind data record */
         /* field 4 must be blank */
         if (csa->deck)
         {  read_field(csa);
            if (csa->field[0] != '\0')
               error(csa, "field 4 must be blank\n");
         }
         else
            csa->fldno++;
         /* field 5: keyword 'INTORG' or 'INTEND' */
         read_field(csa), patch_name(csa, csa->field);
         if (strcmp(csa->field, "'INTORG'") == 0)
            kind = GLP_IV;
         else if (strcmp(csa->field, "'INTEND'") == 0)
            kind = GLP_CV;
         else if (csa->field[0] == '\0')
            error(csa, "missing keyword in field 5\n");
         else
            error(csa, "invalid keyword in field 5\n");
         /* field 6 must be blank */
         skip_field(csa);
         goto loop;
      }
      /* process column name specified in field 2 */
      if (name[0] == '\0')
      {  /* the same column as in previous data record */
         if (j == 0)
            error(csa, "missing column name in field 2\n");
      }
      else if (j != 0 && strcmp(name, csa->P->col[j]->name) == 0)
      {  /* the same column as in previous data record */
         xassert(j != 0);
      }
      else
      {  /* store the current column */
         if (j != 0)
         {  glp_set_mat_col(csa->P, j, len, ind, val);
            while (len > 0) flag[ind[len--]] = 0;
         }
         /* create new column */
         if (glp_find_col(csa->P, name) != 0)
            error(csa, "column '%s' multiply specified\n", name);
         j = glp_add_cols(csa->P, 1);
         glp_set_col_name(csa->P, j, name);
         glp_set_col_kind(csa->P, j, kind);
         if (kind == GLP_CV)
            glp_set_col_bnds(csa->P, j, GLP_LO, 0.0, 0.0);
         else if (kind == GLP_IV)
            glp_set_col_bnds(csa->P, j, GLP_DB, 0.0, 1.0);
         else
            xassert(kind != kind);
      }
      /* process fields 3-4 and 5-6 */
      for (f = 3; f <= 5; f += 2)
      {  /* field 3 or 5: row name */
         if (f == 3)
         {  if (csa->field[0] == '\0')
               error(csa, "missing row name in field 3\n");
         }
         else
         {  read_field(csa), patch_name(csa, csa->field);
            if (csa->field[0] == '\0')
            {  /* if field 5 is blank, field 6 also must be blank */
               skip_field(csa);
               continue;
            }
         }
         i = glp_find_row(csa->P, csa->field);
         if (i == 0)
            error(csa, "row '%s' not found\n", csa->field);
         if (flag[i])
            error(csa, "duplicate coefficient in row '%s'\n",
               csa->field);
         /* field 4 or 6: coefficient value */
         aij = read_number(csa);
         if (fabs(aij) < csa->parm->tol_mps) aij = 0.0;
         len++, ind[len] = i, val[len] = aij, flag[i] = 1;
      }
      goto loop;
done: /* store the last column */
      if (j != 0)
         glp_set_mat_col(csa->P, j, len, ind, val);
      /* free working arrays */
      xfree(ind);
      xfree(val);
      xfree(flag);
      csa->work1 = csa->work2 = csa->work3 = NULL;
      return;
}

static void read_rhs(struct csa *csa)
{     /* read RHS section */
      int i, f, v, type;
      double rhs;
      char name[255+1], *flag;
      /* allocate working array */
      csa->work3 = flag = xcalloc(1+csa->P->m, sizeof(char));
      memset(&flag[1], 0, csa->P->m);
      /* no current RHS vector exists */
      v = 0;
loop: if (indicator(csa, 0)) goto done;
      /* field 1 must be blank */
      if (csa->deck)
      {  read_field(csa);
         if (csa->field[0] != '\0')
            error(csa, "field 1 must be blank\n");
      }
      else
         csa->fldno++;
      /* field 2: RHS vector name */
      read_field(csa), patch_name(csa, csa->field);
      if (csa->field[0] == '\0')
      {  /* the same RHS vector as in previous data record */
         if (v == 0)
         {  warning(csa, "missing RHS vector name in field 2\n");
            goto blnk;
         }
      }
      else if (v != 0 && strcmp(csa->field, name) == 0)
      {  /* the same RHS vector as in previous data record */
         xassert(v != 0);
      }
      else
blnk: {  /* new RHS vector */
         if (v != 0)
            error(csa, "multiple RHS vectors not supported\n");
         v++;
         strcpy(name, csa->field);
      }
      /* process fields 3-4 and 5-6 */
      for (f = 3; f <= 5; f += 2)
      {  /* field 3 or 5: row name */
         read_field(csa), patch_name(csa, csa->field);
         if (csa->field[0] == '\0')
         {  if (f == 3)
               error(csa, "missing row name in field 3\n");
            else
            {  /* if field 5 is blank, field 6 also must be blank */
               skip_field(csa);
               continue;
            }
         }
         i = glp_find_row(csa->P, csa->field);
         if (i == 0)
            error(csa, "row '%s' not found\n", csa->field);
         if (flag[i])
            error(csa, "duplicate right-hand side for row '%s'\n",
               csa->field);
         /* field 4 or 6: right-hand side value */
         rhs = read_number(csa);
         if (fabs(rhs) < csa->parm->tol_mps) rhs = 0.0;
         type = csa->P->row[i]->type;
         if (type == GLP_FR)
         {  if (i == csa->obj_row)
               glp_set_obj_coef(csa->P, 0, rhs);
            else if (rhs != 0.0)
               warning(csa, "non-zero right-hand side for free row '%s'"
                  " ignored\n", csa->P->row[i]->name);
         }
         else
            glp_set_row_bnds(csa->P, i, type, rhs, rhs);
         flag[i] = 1;
      }
      goto loop;
done: /* free working array */
      xfree(flag);
      csa->work3 = NULL;
      return;
}

static void read_ranges(struct csa *csa)
{     /* read RANGES section */
      int i, f, v, type;
      double rhs, rng;
      char name[255+1], *flag;
      /* allocate working array */
      csa->work3 = flag = xcalloc(1+csa->P->m, sizeof(char));
      memset(&flag[1], 0, csa->P->m);
      /* no current RANGES vector exists */
      v = 0;
loop: if (indicator(csa, 0)) goto done;
      /* field 1 must be blank */
      if (csa->deck)
      {  read_field(csa);
         if (csa->field[0] != '\0')
            error(csa, "field 1 must be blank\n");
      }
      else
         csa->fldno++;
      /* field 2: RANGES vector name */
      read_field(csa), patch_name(csa, csa->field);
      if (csa->field[0] == '\0')
      {  /* the same RANGES vector as in previous data record */
         if (v == 0)
         {  warning(csa, "missing RANGES vector name in field 2\n");
            goto blnk;
         }
      }
      else if (v != 0 && strcmp(csa->field, name) == 0)
      {  /* the same RANGES vector as in previous data record */
         xassert(v != 0);
      }
      else
blnk: {  /* new RANGES vector */
         if (v != 0)
            error(csa, "multiple RANGES vectors not supported\n");
         v++;
         strcpy(name, csa->field);
      }
      /* process fields 3-4 and 5-6 */
      for (f = 3; f <= 5; f += 2)
      {  /* field 3 or 5: row name */
         read_field(csa), patch_name(csa, csa->field);
         if (csa->field[0] == '\0')
         {  if (f == 3)
               error(csa, "missing row name in field 3\n");
            else
            {  /* if field 5 is blank, field 6 also must be blank */
               skip_field(csa);
               continue;
            }
         }
         i = glp_find_row(csa->P, csa->field);
         if (i == 0)
            error(csa, "row '%s' not found\n", csa->field);
         if (flag[i])
            error(csa, "duplicate range for row '%s'\n", csa->field);
         /* field 4 or 6: range value */
         rng = read_number(csa);
         if (fabs(rng) < csa->parm->tol_mps) rng = 0.0;
         type = csa->P->row[i]->type;
         if (type == GLP_FR)
            warning(csa, "range for free row '%s' ignored\n",
               csa->P->row[i]->name);
         else if (type == GLP_LO)
         {  rhs = csa->P->row[i]->lb;
#if 0 /* 26/V-2017 by cmatraki */
            glp_set_row_bnds(csa->P, i, rhs == 0.0 ? GLP_FX : GLP_DB,
#else
            glp_set_row_bnds(csa->P, i, rng == 0.0 ? GLP_FX : GLP_DB,
#endif
               rhs, rhs + fabs(rng));
         }
         else if (type == GLP_UP)
         {  rhs = csa->P->row[i]->ub;
#if 0 /* 26/V-2017 by cmatraki */
            glp_set_row_bnds(csa->P, i, rhs == 0.0 ? GLP_FX : GLP_DB,
#else
            glp_set_row_bnds(csa->P, i, rng == 0.0 ? GLP_FX : GLP_DB,
#endif
               rhs - fabs(rng), rhs);
         }
         else if (type == GLP_FX)
         {  rhs = csa->P->row[i]->lb;
            if (rng > 0.0)
               glp_set_row_bnds(csa->P, i, GLP_DB, rhs, rhs + rng);
            else if (rng < 0.0)
               glp_set_row_bnds(csa->P, i, GLP_DB, rhs + rng, rhs);
         }
         else
            xassert(type != type);
         flag[i] = 1;
      }
      goto loop;
done: /* free working array */
      xfree(flag);
      csa->work3 = NULL;
      return;
}

static void read_bounds(struct csa *csa)
{     /* read BOUNDS section */
      GLPCOL *col;
      int j, v, mask, data;
      double bnd, lb, ub;
      char type[2+1], name[255+1], *flag;
      /* allocate working array */
      csa->work3 = flag = xcalloc(1+csa->P->n, sizeof(char));
      memset(&flag[1], 0, csa->P->n);
      /* no current BOUNDS vector exists */
      v = 0;
loop: if (indicator(csa, 0)) goto done;
      /* field 1: bound type */
      read_field(csa);
      if (strcmp(csa->field, "LO") == 0)
         mask = 0x01, data = 1;
      else if (strcmp(csa->field, "UP") == 0)
         mask = 0x10, data = 1;
      else if (strcmp(csa->field, "FX") == 0)
         mask = 0x11, data = 1;
      else if (strcmp(csa->field, "FR") == 0)
         mask = 0x11, data = 0;
      else if (strcmp(csa->field, "MI") == 0)
         mask = 0x01, data = 0;
      else if (strcmp(csa->field, "PL") == 0)
         mask = 0x10, data = 0;
      else if (strcmp(csa->field, "LI") == 0)
         mask = 0x01, data = 1;
      else if (strcmp(csa->field, "UI") == 0)
         mask = 0x10, data = 1;
      else if (strcmp(csa->field, "BV") == 0)
         mask = 0x11, data = 0;
      else if (csa->field[0] == '\0')
         error(csa, "missing bound type in field 1\n");
      else
         error(csa, "invalid bound type in field 1\n");
      strcpy(type, csa->field);
      /* field 2: BOUNDS vector name */
      read_field(csa), patch_name(csa, csa->field);
      if (csa->field[0] == '\0')
      {  /* the same BOUNDS vector as in previous data record */
         if (v == 0)
         {  warning(csa, "missing BOUNDS vector name in field 2\n");
            goto blnk;
         }
      }
      else if (v != 0 && strcmp(csa->field, name) == 0)
      {  /* the same BOUNDS vector as in previous data record */
         xassert(v != 0);
      }
      else
blnk: {  /* new BOUNDS vector */
         if (v != 0)
            error(csa, "multiple BOUNDS vectors not supported\n");
         v++;
         strcpy(name, csa->field);
      }
      /* field 3: column name */
      read_field(csa), patch_name(csa, csa->field);
      if (csa->field[0] == '\0')
         error(csa, "missing column name in field 3\n");
      j = glp_find_col(csa->P, csa->field);
      if (j == 0)
         error(csa, "column '%s' not found\n", csa->field);
      if ((flag[j] & mask) == 0x01)
         error(csa, "duplicate lower bound for column '%s'\n",
            csa->field);
      if ((flag[j] & mask) == 0x10)
         error(csa, "duplicate upper bound for column '%s'\n",
            csa->field);
      xassert((flag[j] & mask) == 0x00);
      /* field 4: bound value */
      if (data)
      {  bnd = read_number(csa);
         if (fabs(bnd) < csa->parm->tol_mps) bnd = 0.0;
      }
      else
         read_field(csa), bnd = 0.0;
      /* get current column bounds */
      col = csa->P->col[j];
      if (col->type == GLP_FR)
         lb = -DBL_MAX, ub = +DBL_MAX;
      else if (col->type == GLP_LO)
         lb = col->lb, ub = +DBL_MAX;
      else if (col->type == GLP_UP)
         lb = -DBL_MAX, ub = col->ub;
      else if (col->type == GLP_DB)
         lb = col->lb, ub = col->ub;
      else if (col->type == GLP_FX)
         lb = ub = col->lb;
      else
         xassert(col != col);
      /* change column bounds */
      if (strcmp(type, "LO") == 0)
         lb = bnd;
      else if (strcmp(type, "UP") == 0)
         ub = bnd;
      else if (strcmp(type, "FX") == 0)
         lb = ub = bnd;
      else if (strcmp(type, "FR") == 0)
         lb = -DBL_MAX, ub = +DBL_MAX;
      else if (strcmp(type, "MI") == 0)
         lb = -DBL_MAX;
      else if (strcmp(type, "PL") == 0)
         ub = +DBL_MAX;
      else if (strcmp(type, "LI") == 0)
      {  glp_set_col_kind(csa->P, j, GLP_IV);
         lb = ceil(bnd);
#if 1 /* 16/VII-2013 */
         /* if column upper bound has not been explicitly specified,
            take it as +inf */
         if (!(flag[j] & 0x10))
            ub = +DBL_MAX;
#endif
      }
      else if (strcmp(type, "UI") == 0)
      {  glp_set_col_kind(csa->P, j, GLP_IV);
         ub = floor(bnd);
      }
      else if (strcmp(type, "BV") == 0)
      {  glp_set_col_kind(csa->P, j, GLP_IV);
         lb = 0.0, ub = 1.0;
      }
      else
         xassert(type != type);
      /* set new column bounds */
      if (lb == -DBL_MAX && ub == +DBL_MAX)
         glp_set_col_bnds(csa->P, j, GLP_FR, lb, ub);
      else if (ub == +DBL_MAX)
         glp_set_col_bnds(csa->P, j, GLP_LO, lb, ub);
      else if (lb == -DBL_MAX)
         glp_set_col_bnds(csa->P, j, GLP_UP, lb, ub);
      else if (lb != ub)
         glp_set_col_bnds(csa->P, j, GLP_DB, lb, ub);
      else
         glp_set_col_bnds(csa->P, j, GLP_FX, lb, ub);
      flag[j] |= (char)mask;
      /* fields 5 and 6 must be blank */
      skip_field(csa);
      skip_field(csa);
      goto loop;
done: /* free working array */
      xfree(flag);
      csa->work3 = NULL;
      return;
}

int glp_read_mps(glp_prob *P, int fmt, const glp_mpscp *parm,
      const char *fname)
{     /* read problem data in MPS format */
      glp_mpscp _parm;
      struct csa _csa, *csa = &_csa;
      int ret;
      xprintf("Reading problem data from '%s'...\n", fname);
      if (!(fmt == GLP_MPS_DECK || fmt == GLP_MPS_FILE))
         xerror("glp_read_mps: fmt = %d; invalid parameter\n", fmt);
      if (parm == NULL)
         glp_init_mpscp(&_parm), parm = &_parm;
      /* check control parameters */
      check_parm("glp_read_mps", parm);
      /* initialize common storage area */
      csa->P = P;
      csa->deck = (fmt == GLP_MPS_DECK);
      csa->parm = parm;
      csa->fname = fname;
      csa->fp = NULL;
      if (setjmp(csa->jump))
      {  ret = 1;
         goto done;
      }
      csa->recno = csa->recpos = 0;
      csa->c = '\n';
      csa->fldno = 0;
      csa->field[0] = '\0';
      csa->w80 = csa->wef = 0;
      csa->obj_row = 0;
      csa->work1 = csa->work2 = csa->work3 = NULL;
      /* erase problem object */
      glp_erase_prob(P);
      glp_create_index(P);
      /* open input MPS file */
      csa->fp = glp_open(fname, "r");
      if (csa->fp == NULL)
      {  xprintf("Unable to open '%s' - %s\n", fname, get_err_msg());
         ret = 1;
         goto done;
      }
      /* read NAME indicator record */
      read_name(csa);
      if (P->name != NULL)
         xprintf("Problem: %s\n", P->name);
      /* read ROWS section */
      if (!(indicator(csa, 0) && strcmp(csa->field, "ROWS") == 0))
         error(csa, "missing ROWS indicator record\n");
      read_rows(csa);
      /* determine objective row */
      if (parm->obj_name == NULL || parm->obj_name[0] == '\0')
      {  /* use the first row of N type */
         int i;
         for (i = 1; i <= P->m; i++)
         {  if (P->row[i]->type == GLP_FR)
            {  csa->obj_row = i;
               break;
            }
         }
         if (csa->obj_row == 0)
            warning(csa, "unable to determine objective row\n");
      }
      else
      {  /* use a row with specified name */
         int i;
         for (i = 1; i <= P->m; i++)
         {  xassert(P->row[i]->name != NULL);
            if (strcmp(parm->obj_name, P->row[i]->name) == 0)
            {  csa->obj_row = i;
               break;
            }
         }
         if (csa->obj_row == 0)
            error(csa, "objective row '%s' not found\n",
               parm->obj_name);
      }
      if (csa->obj_row != 0)
      {  glp_set_obj_name(P, P->row[csa->obj_row]->name);
         xprintf("Objective: %s\n", P->obj);
      }
      /* read COLUMNS section */
      if (strcmp(csa->field, "COLUMNS") != 0)
         error(csa, "missing COLUMNS indicator record\n");
      read_columns(csa);
      /* set objective coefficients */
      if (csa->obj_row != 0)
      {  GLPAIJ *aij;
         for (aij = P->row[csa->obj_row]->ptr; aij != NULL; aij =
            aij->r_next) glp_set_obj_coef(P, aij->col->j, aij->val);
      }
      /* read optional RHS section */
      if (strcmp(csa->field, "RHS") == 0)
         read_rhs(csa);
      /* read optional RANGES section */
      if (strcmp(csa->field, "RANGES") == 0)
         read_ranges(csa);
      /* read optional BOUNDS section */
      if (strcmp(csa->field, "BOUNDS") == 0)
         read_bounds(csa);
      /* read ENDATA indicator record */
      if (strcmp(csa->field, "ENDATA") != 0)
         error(csa, "invalid use of %s indicator record\n",
            csa->field);
      /* print some statistics */
      xprintf("%d row%s, %d column%s, %d non-zero%s\n",
         P->m, P->m == 1 ? "" : "s", P->n, P->n == 1 ? "" : "s",
         P->nnz, P->nnz == 1 ? "" : "s");
      if (glp_get_num_int(P) > 0)
      {  int ni = glp_get_num_int(P);
         int nb = glp_get_num_bin(P);
         if (ni == 1)
         {  if (nb == 0)
               xprintf("One variable is integer\n");
            else
               xprintf("One variable is binary\n");
         }
         else
         {  xprintf("%d integer variables, ", ni);
            if (nb == 0)
               xprintf("none");
            else if (nb == 1)
               xprintf("one");
            else if (nb == ni)
               xprintf("all");
            else
               xprintf("%d", nb);
            xprintf(" of which %s binary\n", nb == 1 ? "is" : "are");
         }
      }
      xprintf("%d records were read\n", csa->recno);
#if 1 /* 31/III-2016 */
      /* free (unbounded) row(s) in MPS file are intended to specify
       * objective function(s), so all such rows can be removed */
#if 1 /* 08/VIII-2013 */
      /* remove free rows */
      {  int i, nrs, *num;
         num = talloc(1+P->m, int);
         nrs = 0;
         for (i = 1; i <= P->m; i++)
         {  if (P->row[i]->type == GLP_FR)
               num[++nrs] = i;
         }
         if (nrs > 0)
         {  glp_del_rows(P, nrs, num);
            if (nrs == 1)
               xprintf("One free row was removed\n");
            else
               xprintf("%d free rows were removed\n", nrs);
         }
         tfree(num);
      }
#endif
#else
      /* if objective function row is free, remove it */
      if (csa->obj_row != 0 && P->row[csa->obj_row]->type == GLP_FR)
      {  int num[1+1];
         num[1] = csa->obj_row;
         glp_del_rows(P, 1, num);
         xprintf("Free objective row was removed\n");
      }
#endif
      /* problem data has been successfully read */
      glp_delete_index(P);
      glp_sort_matrix(P);
      ret = 0;
done: if (csa->fp != NULL) glp_close(csa->fp);
      if (csa->work1 != NULL) xfree(csa->work1);
      if (csa->work2 != NULL) xfree(csa->work2);
      if (csa->work3 != NULL) xfree(csa->work3);
      if (ret != 0) glp_erase_prob(P);
      return ret;
}

/***********************************************************************
*  NAME
*
*  glp_write_mps - write problem data in MPS format
*
*  SYNOPSIS
*
*  int glp_write_mps(glp_prob *P, int fmt, const glp_mpscp *parm,
*     const char *fname);
*
*  DESCRIPTION
*
*  The routine glp_write_mps writes problem data in MPS format to a
*  text file.
*
*  The parameter fmt specifies the version of MPS format:
*
*  GLP_MPS_DECK - fixed (ancient) MPS format;
*  GLP_MPS_FILE - free (modern) MPS format.
*
*  The parameter parm is a pointer to the structure glp_mpscp, which
*  specifies control parameters used by the routine. If parm is NULL,
*  the routine uses default settings.
*
*  The character string fname specifies a name of the text file to be
*  written.
*
*  RETURNS
*
*  If the operation was successful, the routine glp_read_mps returns
*  zero. Otherwise, it prints an error message and returns non-zero. */

#define csa csa1

struct csa
{     /* common storage area */
      glp_prob *P;
      /* pointer to problem object */
      int deck;
      /* MPS format (0 - free, 1 - fixed) */
      const glp_mpscp *parm;
      /* pointer to control parameters */
      char field[255+1];
      /* field buffer */
};

static char *mps_name(struct csa *csa)
{     /* make problem name */
      char *f;
      if (csa->P->name == NULL)
         csa->field[0] = '\0';
      else if (csa->deck)
      {  strncpy(csa->field, csa->P->name, 8);
         csa->field[8] = '\0';
      }
      else
         strcpy(csa->field, csa->P->name);
      for (f = csa->field; *f != '\0'; f++)
         if (*f == ' ') *f = '_';
      return csa->field;
}

static char *row_name(struct csa *csa, int i)
{     /* make i-th row name */
      char *f;
      xassert(0 <= i && i <= csa->P->m);
      if (i == 0 || csa->P->row[i]->name == NULL ||
          csa->deck && strlen(csa->P->row[i]->name) > 8)
         sprintf(csa->field, "R%07d", i);
      else
      {  strcpy(csa->field, csa->P->row[i]->name);
         for (f = csa->field; *f != '\0'; f++)
            if (*f == ' ') *f = '_';
      }
      return csa->field;
}

static char *col_name(struct csa *csa, int j)
{     /* make j-th column name */
      char *f;
      xassert(1 <= j && j <= csa->P->n);
      if (csa->P->col[j]->name == NULL ||
          csa->deck && strlen(csa->P->col[j]->name) > 8)
         sprintf(csa->field, "C%07d", j);
      else
      {  strcpy(csa->field, csa->P->col[j]->name);
         for (f = csa->field; *f != '\0'; f++)
            if (*f == ' ') *f = '_';
      }
      return csa->field;
}

static char *mps_numb(struct csa *csa, double val)
{     /* format floating-point number */
      int dig;
      char *exp;
      for (dig = 12; dig >= 6; dig--)
      {  if (val != 0.0 && fabs(val) < 0.002)
            sprintf(csa->field, "%.*E", dig-1, val);
         else
            sprintf(csa->field, "%.*G", dig, val);
         exp = strchr(csa->field, 'E');
         if (exp != NULL)
            sprintf(exp+1, "%d", atoi(exp+1));
         if (strlen(csa->field) <= 12) break;
      }
      xassert(strlen(csa->field) <= 12);
      return csa->field;
}

int glp_write_mps(glp_prob *P, int fmt, const glp_mpscp *parm,
      const char *fname)
{     /* write problem data in MPS format */
      glp_mpscp _parm;
      struct csa _csa, *csa = &_csa;
      glp_file *fp;
      int out_obj, one_col = 0, empty = 0;
      int i, j, recno, marker, count, gap, ret;
      xprintf("Writing problem data to '%s'...\n", fname);
      if (!(fmt == GLP_MPS_DECK || fmt == GLP_MPS_FILE))
         xerror("glp_write_mps: fmt = %d; invalid parameter\n", fmt);
      if (parm == NULL)
         glp_init_mpscp(&_parm), parm = &_parm;
      /* check control parameters */
      check_parm("glp_write_mps", parm);
      /* initialize common storage area */
      csa->P = P;
      csa->deck = (fmt == GLP_MPS_DECK);
      csa->parm = parm;
      /* create output MPS file */
      fp = glp_open(fname, "w"), recno = 0;
      if (fp == NULL)
      {  xprintf("Unable to create '%s' - %s\n", fname, get_err_msg());
         ret = 1;
         goto done;
      }
      /* write comment records */
      xfprintf(fp, "* %-*s%s\n", P->name == NULL ? 1 : 12, "Problem:",
         P->name == NULL ? "" : P->name), recno++;
      xfprintf(fp, "* %-12s%s\n", "Class:", glp_get_num_int(P) == 0 ?
         "LP" : "MIP"), recno++;
      xfprintf(fp, "* %-12s%d\n", "Rows:", P->m), recno++;
      if (glp_get_num_int(P) == 0)
         xfprintf(fp, "* %-12s%d\n", "Columns:", P->n), recno++;
      else
         xfprintf(fp, "* %-12s%d (%d integer, %d binary)\n",
            "Columns:", P->n, glp_get_num_int(P), glp_get_num_bin(P)),
            recno++;
      xfprintf(fp, "* %-12s%d\n", "Non-zeros:", P->nnz), recno++;
      xfprintf(fp, "* %-12s%s\n", "Format:", csa->deck ? "Fixed MPS" :
         "Free MPS"), recno++;
      xfprintf(fp, "*\n", recno++);
      /* write NAME indicator record */
      xfprintf(fp, "NAME%*s%s\n",
         P->name == NULL ? 0 : csa->deck ? 10 : 1, "", mps_name(csa)),
         recno++;
#if 1
      /* determine whether to write the objective row */
      out_obj = 1;
      for (i = 1; i <= P->m; i++)
      {  if (P->row[i]->type == GLP_FR)
         {  out_obj = 0;
            break;
         }
      }
#endif
      /* write ROWS section */
      xfprintf(fp, "ROWS\n"), recno++;
      for (i = (out_obj ? 0 : 1); i <= P->m; i++)
      {  int type;
         type = (i == 0 ? GLP_FR : P->row[i]->type);
         if (type == GLP_FR)
            type = 'N';
         else if (type == GLP_LO)
            type = 'G';
         else if (type == GLP_UP)
            type = 'L';
         else if (type == GLP_DB || type == GLP_FX)
            type = 'E';
         else
            xassert(type != type);
         xfprintf(fp, " %c%*s%s\n", type, csa->deck ? 2 : 1, "",
            row_name(csa, i)), recno++;
      }
      /* write COLUMNS section */
      xfprintf(fp, "COLUMNS\n"), recno++;
      marker = 0;
      for (j = 1; j <= P->n; j++)
      {  GLPAIJ cj, *aij;
         int kind;
         kind = P->col[j]->kind;
         if (kind == GLP_CV)
         {  if (marker % 2 == 1)
            {  /* close current integer block */
               marker++;
               xfprintf(fp, "%*sM%07d%*s'MARKER'%*s'INTEND'\n",
                  csa->deck ? 4 : 1, "", marker, csa->deck ? 2 : 1, "",
                  csa->deck ? 17 : 1, ""), recno++;
            }
         }
         else if (kind == GLP_IV)
         {  if (marker % 2 == 0)
            {  /* open new integer block */
               marker++;
               xfprintf(fp, "%*sM%07d%*s'MARKER'%*s'INTORG'\n",
                  csa->deck ? 4 : 1, "", marker, csa->deck ? 2 : 1, "",
                  csa->deck ? 17 : 1, ""), recno++;
            }
         }
         else
            xassert(kind != kind);
         if (out_obj && P->col[j]->coef != 0.0)
         {  /* make fake objective coefficient */
            aij = &cj;
            aij->row = NULL;
            aij->val = P->col[j]->coef;
            aij->c_next = P->col[j]->ptr;
         }
         else
            aij = P->col[j]->ptr;
#if 1 /* FIXME */
         if (aij == NULL)
         {  /* empty column */
            empty++;
            xfprintf(fp, "%*s%-*s", csa->deck ? 4 : 1, "",
               csa->deck ? 8 : 1, col_name(csa, j));
            /* we need a row */
            xassert(P->m > 0);
            xfprintf(fp, "%*s%-*s",
               csa->deck ? 2 : 1, "", csa->deck ? 8 : 1,
               row_name(csa, 1));
            xfprintf(fp, "%*s0%*s$ empty column\n",
               csa->deck ? 13 : 1, "", csa->deck ? 3 : 1, ""), recno++;
         }
#endif
         count = 0;
         for (aij = aij; aij != NULL; aij = aij->c_next)
         {  if (one_col || count % 2 == 0)
               xfprintf(fp, "%*s%-*s", csa->deck ? 4 : 1, "",
                  csa->deck ? 8 : 1, col_name(csa, j));
            gap = (one_col || count % 2 == 0 ? 2 : 3);
            xfprintf(fp, "%*s%-*s",
               csa->deck ? gap : 1, "", csa->deck ? 8 : 1,
               row_name(csa, aij->row == NULL ? 0 : aij->row->i));
            xfprintf(fp, "%*s%*s",
               csa->deck ? 2 : 1, "", csa->deck ? 12 : 1,
               mps_numb(csa, aij->val)), count++;
            if (one_col || count % 2 == 0)
               xfprintf(fp, "\n"), recno++;
         }
         if (!(one_col || count % 2 == 0))
            xfprintf(fp, "\n"), recno++;
      }
      if (marker % 2 == 1)
      {  /* close last integer block */
         marker++;
         xfprintf(fp, "%*sM%07d%*s'MARKER'%*s'INTEND'\n",
            csa->deck ? 4 : 1, "", marker, csa->deck ? 2 : 1, "",
            csa->deck ? 17 : 1, ""), recno++;
      }
#if 1
      if (empty > 0)
         xprintf("Warning: problem has %d empty column(s)\n", empty);
#endif
      /* write RHS section */
      xfprintf(fp, "RHS\n"), recno++;
      count = 0;
      for (i = (out_obj ? 0 : 1); i <= P->m; i++)
      {  int type;
         double rhs;
         if (i == 0)
            rhs = P->c0;
         else
         {  type = P->row[i]->type;
            if (type == GLP_FR)
               rhs = 0.0;
            else if (type == GLP_LO)
               rhs = P->row[i]->lb;
            else if (type == GLP_UP)
               rhs = P->row[i]->ub;
            else if (type == GLP_DB || type == GLP_FX)
               rhs = P->row[i]->lb;
            else
               xassert(type != type);
         }
         if (rhs != 0.0)
         {  if (one_col || count % 2 == 0)
               xfprintf(fp, "%*s%-*s", csa->deck ? 4 : 1, "",
                  csa->deck ? 8 : 1, "RHS1");
            gap = (one_col || count % 2 == 0 ? 2 : 3);
            xfprintf(fp, "%*s%-*s",
               csa->deck ? gap : 1, "", csa->deck ? 8 : 1,
               row_name(csa, i));
            xfprintf(fp, "%*s%*s",
               csa->deck ? 2 : 1, "", csa->deck ? 12 : 1,
               mps_numb(csa, rhs)), count++;
            if (one_col || count % 2 == 0)
               xfprintf(fp, "\n"), recno++;
         }
      }
      if (!(one_col || count % 2 == 0))
         xfprintf(fp, "\n"), recno++;
      /* write RANGES section */
      for (i = P->m; i >= 1; i--)
         if (P->row[i]->type == GLP_DB) break;
      if (i == 0) goto bnds;
      xfprintf(fp, "RANGES\n"), recno++;
      count = 0;
      for (i = 1; i <= P->m; i++)
      {  if (P->row[i]->type == GLP_DB)
         {  if (one_col || count % 2 == 0)
               xfprintf(fp, "%*s%-*s", csa->deck ? 4 : 1, "",
                  csa->deck ? 8 : 1, "RNG1");
            gap = (one_col || count % 2 == 0 ? 2 : 3);
            xfprintf(fp, "%*s%-*s",
               csa->deck ? gap : 1, "", csa->deck ? 8 : 1,
               row_name(csa, i));
            xfprintf(fp, "%*s%*s",
               csa->deck ? 2 : 1, "", csa->deck ? 12 : 1,
               mps_numb(csa, P->row[i]->ub - P->row[i]->lb)), count++;
            if (one_col || count % 2 == 0)
               xfprintf(fp, "\n"), recno++;
         }
      }
      if (!(one_col || count % 2 == 0))
         xfprintf(fp, "\n"), recno++;
bnds: /* write BOUNDS section */
      for (j = P->n; j >= 1; j--)
         if (!(P->col[j]->kind == GLP_CV &&
               P->col[j]->type == GLP_LO && P->col[j]->lb == 0.0))
            break;
      if (j == 0) goto endt;
      xfprintf(fp, "BOUNDS\n"), recno++;
      for (j = 1; j <= P->n; j++)
      {  int type, data[2];
         double bnd[2];
         char *spec[2];
         spec[0] = spec[1] = NULL;
         type = P->col[j]->type;
         if (type == GLP_FR)
            spec[0] = "FR", data[0] = 0;
         else if (type == GLP_LO)
         {  if (P->col[j]->lb != 0.0)
               spec[0] = "LO", data[0] = 1, bnd[0] = P->col[j]->lb;
            if (P->col[j]->kind == GLP_IV)
               spec[1] = "PL", data[1] = 0;
         }
         else if (type == GLP_UP)
         {  spec[0] = "MI", data[0] = 0;
            spec[1] = "UP", data[1] = 1, bnd[1] = P->col[j]->ub;
         }
         else if (type == GLP_DB)
         {  if (P->col[j]->lb != 0.0)
               spec[0] = "LO", data[0] = 1, bnd[0] = P->col[j]->lb;
            spec[1] = "UP", data[1] = 1, bnd[1] = P->col[j]->ub;
         }
         else if (type == GLP_FX)
            spec[0] = "FX", data[0] = 1, bnd[0] = P->col[j]->lb;
         else
            xassert(type != type);
         for (i = 0; i <= 1; i++)
         {  if (spec[i] != NULL)
            {  xfprintf(fp, " %s %-*s%*s%-*s", spec[i],
                  csa->deck ? 8 : 1, "BND1", csa->deck ? 2 : 1, "",
                  csa->deck ? 8 : 1, col_name(csa, j));
               if (data[i])
                  xfprintf(fp, "%*s%*s", csa->deck ? 2 : 1, "",
                     csa->deck ? 12 : 1, mps_numb(csa, bnd[i]));
               xfprintf(fp, "\n"), recno++;
            }
         }
      }
endt: /* write ENDATA indicator record */
      xfprintf(fp, "ENDATA\n"), recno++;
#if 0 /* FIXME */
      xfflush(fp);
#endif
      if (glp_ioerr(fp))
      {  xprintf("Write error on '%s' - %s\n", fname, get_err_msg());
         ret = 1;
         goto done;
      }
      /* problem data has been successfully written */
      xprintf("%d records were written\n", recno);
      ret = 0;
done: if (fp != NULL) glp_close(fp);
      return ret;
}

/* eof */
