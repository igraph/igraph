/* cplex.c (CPLEX LP format routines) */

/***********************************************************************
*  This code is part of GLPK (GNU Linear Programming Kit).
*  Copyright (C) 2009-2018 Free Software Foundation, Inc.
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
*  glp_init_cpxcp - initialize CPLEX LP format control parameters
*
*  SYNOPSIS
*
*  void glp_init_cpxcp(glp_cpxcp *parm):
*
*  The routine glp_init_cpxcp initializes control parameters used by
*  the CPLEX LP input/output routines glp_read_lp and glp_write_lp with
*  default values.
*
*  Default values of the control parameters are stored in the glp_cpxcp
*  structure, which the parameter parm points to. */

void glp_init_cpxcp(glp_cpxcp *parm)
{     xassert(parm != NULL);
      return;
}

static void check_parm(const char *func, const glp_cpxcp *parm)
{     /* check control parameters */
      xassert(func != NULL);
      xassert(parm != NULL);
      return;
}

/***********************************************************************
*  NAME
*
*  glp_read_lp - read problem data in CPLEX LP format
*
*  SYNOPSIS
*
*  int glp_read_lp(glp_prob *P, const glp_cpxcp *parm, const char
*     *fname);
*
*  DESCRIPTION
*
*  The routine glp_read_lp reads problem data in CPLEX LP format from
*  a text file.
*
*  The parameter parm is a pointer to the structure glp_cpxcp, which
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
*  If the operation was successful, the routine glp_read_lp returns
*  zero. Otherwise, it prints an error message and returns non-zero. */

struct csa
{     /* common storage area */
      glp_prob *P;
      /* LP/MIP problem object */
      const glp_cpxcp *parm;
      /* pointer to control parameters */
      const char *fname;
      /* name of input CPLEX LP file */
      glp_file *fp;
      /* stream assigned to input CPLEX LP file */
      jmp_buf jump;
      /* label for go to in case of error */
      int count;
      /* line count */
      int c;
      /* current character or EOF */
      int token;
      /* current token: */
#define T_EOF        0x00  /* end of file */
#define T_MINIMIZE   0x01  /* keyword 'minimize' */
#define T_MAXIMIZE   0x02  /* keyword 'maximize' */
#define T_SUBJECT_TO 0x03  /* keyword 'subject to' */
#define T_BOUNDS     0x04  /* keyword 'bounds' */
#define T_GENERAL    0x05  /* keyword 'general' */
#define T_INTEGER    0x06  /* keyword 'integer' */
#define T_BINARY     0x07  /* keyword 'binary' */
#define T_END        0x08  /* keyword 'end' */
#define T_NAME       0x09  /* symbolic name */
#define T_NUMBER     0x0A  /* numeric constant */
#define T_PLUS       0x0B  /* delimiter '+' */
#define T_MINUS      0x0C  /* delimiter '-' */
#define T_COLON      0x0D  /* delimiter ':' */
#define T_LE         0x0E  /* delimiter '<=' */
#define T_GE         0x0F  /* delimiter '>=' */
#define T_EQ         0x10  /* delimiter '=' */
      char image[255+1];
      /* image of current token */
      int imlen;
      /* length of token image */
      double value;
      /* value of numeric constant */
      int n_max;
      /* length of the following five arrays (enlarged automatically,
         if necessary) */
      int *ind; /* int ind[1+n_max]; */
      double *val; /* double val[1+n_max]; */
      char *flag; /* char flag[1+n_max]; */
      /* working arrays used to construct linear forms */
      double *lb; /* double lb[1+n_max]; */
      double *ub; /* double ub[1+n_max]; */
      /* lower and upper bounds of variables (columns) */
#if 1 /* 27/VII-2013 */
      int lb_warn, ub_warn;
      /* warning 'lower/upper bound redefined' already issued */
#endif
};

#define CHAR_SET "!\"#$%&()/,.;?@_`'{}|~"
/* characters that may appear in symbolic names */

static void error(struct csa *csa, const char *fmt, ...)
{     /* print error message and terminate processing */
      va_list arg;
      xprintf("%s:%d: ", csa->fname, csa->count);
      va_start(arg, fmt);
      xvprintf(fmt, arg);
      va_end(arg);
      longjmp(csa->jump, 1);
      /* no return */
}

static void warning(struct csa *csa, const char *fmt, ...)
{     /* print warning message and continue processing */
      va_list arg;
      xprintf("%s:%d: warning: ", csa->fname, csa->count);
      va_start(arg, fmt);
      xvprintf(fmt, arg);
      va_end(arg);
      return;
}

static void read_char(struct csa *csa)
{     /* read next character from input file */
      int c;
      xassert(csa->c != EOF);
      if (csa->c == '\n') csa->count++;
      c = glp_getc(csa->fp);
      if (c < 0)
      {  if (glp_ioerr(csa->fp))
            error(csa, "read error - %s\n", get_err_msg());
         else if (csa->c == '\n')
         {  csa->count--;
            c = EOF;
         }
         else
         {  warning(csa, "missing final end of line\n");
            c = '\n';
         }
      }
      else if (c == '\n')
         ;
      else if (isspace(c))
         c = ' ';
      else if (iscntrl(c))
         error(csa, "invalid control character 0x%02X\n", c);
      csa->c = c;
      return;
}

static void add_char(struct csa *csa)
{     /* append current character to current token */
      if (csa->imlen == sizeof(csa->image)-1)
         error(csa, "token '%.15s...' too long\n", csa->image);
      csa->image[csa->imlen++] = (char)csa->c;
      csa->image[csa->imlen] = '\0';
      read_char(csa);
      return;
}

static int the_same(char *s1, char *s2)
{     /* compare two character strings ignoring case sensitivity */
      for (; *s1 != '\0'; s1++, s2++)
      {  if (tolower((unsigned char)*s1) != tolower((unsigned char)*s2))
            return 0;
      }
      return 1;
}

static void scan_token(struct csa *csa)
{     /* scan next token */
      int flag;
      csa->token = -1;
      csa->image[0] = '\0';
      csa->imlen = 0;
      csa->value = 0.0;
loop: flag = 0;
      /* skip non-significant characters */
      while (csa->c == ' ') read_char(csa);
      /* recognize and scan current token */
      if (csa->c == EOF)
         csa->token = T_EOF;
      else if (csa->c == '\n')
      {  read_char(csa);
         /* if the next character is letter, it may begin a keyword */
         if (isalpha(csa->c))
         {  flag = 1;
            goto name;
         }
         goto loop;
      }
      else if (csa->c == '\\')
      {  /* comment; ignore everything until end-of-line */
         while (csa->c != '\n') read_char(csa);
         goto loop;
      }
      else if (isalpha(csa->c) || csa->c != '.' && strchr(CHAR_SET,
         csa->c) != NULL)
name: {  /* symbolic name */
         csa->token = T_NAME;
         while (isalnum(csa->c) || strchr(CHAR_SET, csa->c) != NULL)
            add_char(csa);
         if (flag)
         {  /* check for keyword */
            if (the_same(csa->image, "minimize"))
               csa->token = T_MINIMIZE;
            else if (the_same(csa->image, "minimum"))
               csa->token = T_MINIMIZE;
            else if (the_same(csa->image, "min"))
               csa->token = T_MINIMIZE;
            else if (the_same(csa->image, "maximize"))
               csa->token = T_MAXIMIZE;
            else if (the_same(csa->image, "maximum"))
               csa->token = T_MAXIMIZE;
            else if (the_same(csa->image, "max"))
               csa->token = T_MAXIMIZE;
            else if (the_same(csa->image, "subject"))
            {  if (csa->c == ' ')
               {  read_char(csa);
                  if (tolower(csa->c) == 't')
                  {  csa->token = T_SUBJECT_TO;
                     csa->image[csa->imlen++] = ' ';
                     csa->image[csa->imlen] = '\0';
                     add_char(csa);
                     if (tolower(csa->c) != 'o')
                        error(csa, "keyword 'subject to' incomplete\n");
                     add_char(csa);
                     if (isalpha(csa->c))
                        error(csa, "keyword '%s%c...' not recognized\n",
                           csa->image, csa->c);
                  }
               }
            }
            else if (the_same(csa->image, "such"))
            {  if (csa->c == ' ')
               {  read_char(csa);
                  if (tolower(csa->c) == 't')
                  {  csa->token = T_SUBJECT_TO;
                     csa->image[csa->imlen++] = ' ';
                     csa->image[csa->imlen] = '\0';
                     add_char(csa);
                     if (tolower(csa->c) != 'h')
err:                    error(csa, "keyword 'such that' incomplete\n");
                     add_char(csa);
                     if (tolower(csa->c) != 'a') goto err;
                     add_char(csa);
                     if (tolower(csa->c) != 't') goto err;
                     add_char(csa);
                     if (isalpha(csa->c))
                        error(csa, "keyword '%s%c...' not recognized\n",
                           csa->image, csa->c);
                  }
               }
            }
            else if (the_same(csa->image, "st"))
               csa->token = T_SUBJECT_TO;
            else if (the_same(csa->image, "s.t."))
               csa->token = T_SUBJECT_TO;
            else if (the_same(csa->image, "st."))
               csa->token = T_SUBJECT_TO;
            else if (the_same(csa->image, "bounds"))
               csa->token = T_BOUNDS;
            else if (the_same(csa->image, "bound"))
               csa->token = T_BOUNDS;
            else if (the_same(csa->image, "general"))
               csa->token = T_GENERAL;
            else if (the_same(csa->image, "generals"))
               csa->token = T_GENERAL;
            else if (the_same(csa->image, "gen"))
               csa->token = T_GENERAL;
            else if (the_same(csa->image, "integer"))
               csa->token = T_INTEGER;
            else if (the_same(csa->image, "integers"))
               csa->token = T_INTEGER;
            else if (the_same(csa->image, "int"))
              csa->token = T_INTEGER;
            else if (the_same(csa->image, "binary"))
               csa->token = T_BINARY;
            else if (the_same(csa->image, "binaries"))
               csa->token = T_BINARY;
            else if (the_same(csa->image, "bin"))
               csa->token = T_BINARY;
            else if (the_same(csa->image, "end"))
               csa->token = T_END;
         }
      }
      else if (isdigit(csa->c) || csa->c == '.')
      {  /* numeric constant */
         csa->token = T_NUMBER;
         /* scan integer part */
         while (isdigit(csa->c)) add_char(csa);
         /* scan optional fractional part (it is mandatory, if there is
            no integer part) */
         if (csa->c == '.')
         {  add_char(csa);
            if (csa->imlen == 1 && !isdigit(csa->c))
               error(csa, "invalid use of decimal point\n");
            while (isdigit(csa->c)) add_char(csa);
         }
         /* scan optional decimal exponent */
         if (csa->c == 'e' || csa->c == 'E')
         {  add_char(csa);
            if (csa->c == '+' || csa->c == '-') add_char(csa);
            if (!isdigit(csa->c))
               error(csa, "numeric constant '%s' incomplete\n",
                  csa->image);
            while (isdigit(csa->c)) add_char(csa);
         }
         /* convert the numeric constant to floating-point */
         if (str2num(csa->image, &csa->value))
            error(csa, "numeric constant '%s' out of range\n",
               csa->image);
      }
      else if (csa->c == '+')
         csa->token = T_PLUS, add_char(csa);
      else if (csa->c == '-')
         csa->token = T_MINUS, add_char(csa);
      else if (csa->c == ':')
         csa->token = T_COLON, add_char(csa);
      else if (csa->c == '<')
      {  csa->token = T_LE, add_char(csa);
         if (csa->c == '=') add_char(csa);
      }
      else if (csa->c == '>')
      {  csa->token = T_GE, add_char(csa);
         if (csa->c == '=') add_char(csa);
      }
      else if (csa->c == '=')
      {  csa->token = T_EQ, add_char(csa);
         if (csa->c == '<')
            csa->token = T_LE, add_char(csa);
         else if (csa->c == '>')
            csa->token = T_GE, add_char(csa);
      }
      else
         error(csa, "character '%c' not recognized\n", csa->c);
      /* skip non-significant characters */
      while (csa->c == ' ') read_char(csa);
      return;
}

static int find_col(struct csa *csa, char *name)
{     /* find column by its symbolic name */
      int j;
      j = glp_find_col(csa->P, name);
      if (j == 0)
      {  /* not found; create new column */
         j = glp_add_cols(csa->P, 1);
         glp_set_col_name(csa->P, j, name);
         /* enlarge working arrays, if necessary */
         if (csa->n_max < j)
         {  int n_max = csa->n_max;
            int *ind = csa->ind;
            double *val = csa->val;
            char *flag = csa->flag;
            double *lb = csa->lb;
            double *ub = csa->ub;
            csa->n_max += csa->n_max;
            csa->ind = xcalloc(1+csa->n_max, sizeof(int));
            memcpy(&csa->ind[1], &ind[1], n_max * sizeof(int));
            xfree(ind);
            csa->val = xcalloc(1+csa->n_max, sizeof(double));
            memcpy(&csa->val[1], &val[1], n_max * sizeof(double));
            xfree(val);
            csa->flag = xcalloc(1+csa->n_max, sizeof(char));
            memset(&csa->flag[1], 0, csa->n_max * sizeof(char));
            memcpy(&csa->flag[1], &flag[1], n_max * sizeof(char));
            xfree(flag);
            csa->lb = xcalloc(1+csa->n_max, sizeof(double));
            memcpy(&csa->lb[1], &lb[1], n_max * sizeof(double));
            xfree(lb);
            csa->ub = xcalloc(1+csa->n_max, sizeof(double));
            memcpy(&csa->ub[1], &ub[1], n_max * sizeof(double));
            xfree(ub);
         }
         csa->lb[j] = +DBL_MAX, csa->ub[j] = -DBL_MAX;
      }
      return j;
}

/***********************************************************************
*  parse_linear_form - parse linear form
*
*  This routine parses the linear form using the following syntax:
*
*  <variable> ::= <symbolic name>
*  <coefficient> ::= <numeric constant>
*  <term> ::= <variable> | <numeric constant> <variable>
*  <linear form> ::= <term> | + <term> | - <term> |
*     <linear form> + <term> | <linear form> - <term>
*
*  The routine returns the number of terms in the linear form. */

static int parse_linear_form(struct csa *csa)
{     int j, k, len = 0, newlen;
      double s, coef;
loop: /* parse an optional sign */
      if (csa->token == T_PLUS)
         s = +1.0, scan_token(csa);
      else if (csa->token == T_MINUS)
         s = -1.0, scan_token(csa);
      else
         s = +1.0;
      /* parse an optional coefficient */
      if (csa->token == T_NUMBER)
         coef = csa->value, scan_token(csa);
      else
         coef = 1.0;
      /* parse a variable name */
      if (csa->token != T_NAME)
         error(csa, "missing variable name\n");
      /* find the corresponding column */
      j = find_col(csa, csa->image);
      /* check if the variable is already used in the linear form */
      if (csa->flag[j])
         error(csa, "multiple use of variable '%s' not allowed\n",
            csa->image);
      /* add new term to the linear form */
      len++, csa->ind[len] = j, csa->val[len] = s * coef;
      /* and mark that the variable is used in the linear form */
      csa->flag[j] = 1;
      scan_token(csa);
      /* if the next token is a sign, there is another term */
      if (csa->token == T_PLUS || csa->token == T_MINUS) goto loop;
      /* clear marks of the variables used in the linear form */
      for (k = 1; k <= len; k++) csa->flag[csa->ind[k]] = 0;
      /* remove zero coefficients */
      newlen = 0;
      for (k = 1; k <= len; k++)
      {  if (csa->val[k] != 0.0)
         {  newlen++;
            csa->ind[newlen] = csa->ind[k];
            csa->val[newlen] = csa->val[k];
         }
      }
      return newlen;
}

/***********************************************************************
*  parse_objective - parse objective function
*
*  This routine parses definition of the objective function using the
*  following syntax:
*
*  <obj sense> ::= minimize | minimum | min | maximize | maximum | max
*  <obj name> ::= <empty> | <symbolic name> :
*  <obj function> ::= <obj sense> <obj name> <linear form> */

static void parse_objective(struct csa *csa)
{     /* parse objective sense */
      int k, len;
      /* parse the keyword 'minimize' or 'maximize' */
      if (csa->token == T_MINIMIZE)
         glp_set_obj_dir(csa->P, GLP_MIN);
      else if (csa->token == T_MAXIMIZE)
         glp_set_obj_dir(csa->P, GLP_MAX);
      else
         xassert(csa != csa);
      scan_token(csa);
      /* parse objective name */
      if (csa->token == T_NAME && csa->c == ':')
      {  /* objective name is followed by a colon */
         glp_set_obj_name(csa->P, csa->image);
         scan_token(csa);
         xassert(csa->token == T_COLON);
         scan_token(csa);
      }
      else
      {  /* objective name is not specified; use default */
         glp_set_obj_name(csa->P, "obj");
      }
      /* parse linear form */
      len = parse_linear_form(csa);
      for (k = 1; k <= len; k++)
         glp_set_obj_coef(csa->P, csa->ind[k], csa->val[k]);
      return;
}

/***********************************************************************
*  parse_constraints - parse constraints section
*
*  This routine parses the constraints section using the following
*  syntax:
*
*  <row name> ::= <empty> | <symbolic name> :
*  <row sense> ::= < | <= | =< | > | >= | => | =
*  <right-hand side> ::= <numeric constant> | + <numeric constant> |
*     - <numeric constant>
*  <constraint> ::= <row name> <linear form> <row sense>
*     <right-hand side>
*  <subject to> ::= subject to | such that | st | s.t. | st.
*  <constraints section> ::= <subject to> <constraint> |
*     <constraints section> <constraint> */

static void parse_constraints(struct csa *csa)
{     int i, len, type;
      double s;
      /* parse the keyword 'subject to' */
      xassert(csa->token == T_SUBJECT_TO);
      scan_token(csa);
loop: /* create new row (constraint) */
      i = glp_add_rows(csa->P, 1);
      /* parse row name */
      if (csa->token == T_NAME && csa->c == ':')
      {  /* row name is followed by a colon */
         if (glp_find_row(csa->P, csa->image) != 0)
            error(csa, "constraint '%s' multiply defined\n",
               csa->image);
         glp_set_row_name(csa->P, i, csa->image);
         scan_token(csa);
         xassert(csa->token == T_COLON);
         scan_token(csa);
      }
      else
      {  /* row name is not specified; use default */
         char name[50];
         sprintf(name, "r.%d", csa->count);
         glp_set_row_name(csa->P, i, name);
      }
      /* parse linear form */
      len = parse_linear_form(csa);
      glp_set_mat_row(csa->P, i, len, csa->ind, csa->val);
      /* parse constraint sense */
      if (csa->token == T_LE)
         type = GLP_UP, scan_token(csa);
      else if (csa->token == T_GE)
         type = GLP_LO, scan_token(csa);
      else if (csa->token == T_EQ)
         type = GLP_FX, scan_token(csa);
      else
         error(csa, "missing constraint sense\n");
      /* parse right-hand side */
      if (csa->token == T_PLUS)
         s = +1.0, scan_token(csa);
      else if (csa->token == T_MINUS)
         s = -1.0, scan_token(csa);
      else
         s = +1.0;
      if (csa->token != T_NUMBER)
         error(csa, "missing right-hand side\n");
      glp_set_row_bnds(csa->P, i, type, s * csa->value, s * csa->value);
      /* the rest of the current line must be empty */
      if (!(csa->c == '\n' || csa->c == EOF))
         error(csa, "invalid symbol(s) beyond right-hand side\n");
      scan_token(csa);
      /* if the next token is a sign, numeric constant, or a symbolic
         name, here is another constraint */
      if (csa->token == T_PLUS || csa->token == T_MINUS ||
          csa->token == T_NUMBER || csa->token == T_NAME) goto loop;
      return;
}

static void set_lower_bound(struct csa *csa, int j, double lb)
{     /* set lower bound of j-th variable */
      if (csa->lb[j] != +DBL_MAX && !csa->lb_warn)
      {  warning(csa, "lower bound of variable '%s' redefined\n",
            glp_get_col_name(csa->P, j));
         csa->lb_warn = 1;
      }
      csa->lb[j] = lb;
      return;
}

static void set_upper_bound(struct csa *csa, int j, double ub)
{     /* set upper bound of j-th variable */
      if (csa->ub[j] != -DBL_MAX && !csa->ub_warn)
      {  warning(csa, "upper bound of variable '%s' redefined\n",
            glp_get_col_name(csa->P, j));
         csa->ub_warn = 1;
      }
      csa->ub[j] = ub;
      return;
}

/***********************************************************************
*  parse_bounds - parse bounds section
*
*  This routine parses the bounds section using the following syntax:
*
*  <variable> ::= <symbolic name>
*  <infinity> ::= infinity | inf
*  <bound> ::= <numeric constant> | + <numeric constant> |
*     - <numeric constant> | + <infinity> | - <infinity>
*  <lt> ::= < | <= | =<
*  <gt> ::= > | >= | =>
*  <bound definition> ::= <bound> <lt> <variable> <lt> <bound> |
*     <bound> <lt> <variable> | <variable> <lt> <bound> |
*     <variable> <gt> <bound> | <variable> = <bound> | <variable> free
*  <bounds> ::= bounds | bound
*  <bounds section> ::= <bounds> |
*     <bounds section> <bound definition> */

static void parse_bounds(struct csa *csa)
{     int j, lb_flag;
      double lb, s;
      /* parse the keyword 'bounds' */
      xassert(csa->token == T_BOUNDS);
      scan_token(csa);
loop: /* bound definition can start with a sign, numeric constant, or
         a symbolic name */
      if (!(csa->token == T_PLUS || csa->token == T_MINUS ||
            csa->token == T_NUMBER || csa->token == T_NAME)) goto done;
      /* parse bound definition */
      if (csa->token == T_PLUS || csa->token == T_MINUS)
      {  /* parse signed lower bound */
         lb_flag = 1;
         s = (csa->token == T_PLUS ? +1.0 : -1.0);
         scan_token(csa);
         if (csa->token == T_NUMBER)
            lb = s * csa->value, scan_token(csa);
         else if (the_same(csa->image, "infinity") ||
                  the_same(csa->image, "inf"))
         {  if (s > 0.0)
               error(csa, "invalid use of '+inf' as lower bound\n");
            lb = -DBL_MAX, scan_token(csa);
         }
         else
            error(csa, "missing lower bound\n");
      }
      else if (csa->token == T_NUMBER)
      {  /* parse unsigned lower bound */
         lb_flag = 1;
         lb = csa->value, scan_token(csa);
      }
      else
      {  /* lower bound is not specified */
         lb_flag = 0;
      }
      /* parse the token that should follow the lower bound */
      if (lb_flag)
      {  if (csa->token != T_LE)
            error(csa, "missing '<', '<=', or '=<' after lower bound\n")
               ;
         scan_token(csa);
      }
      /* parse variable name */
      if (csa->token != T_NAME)
         error(csa, "missing variable name\n");
      j = find_col(csa, csa->image);
      /* set lower bound */
      if (lb_flag) set_lower_bound(csa, j, lb);
      scan_token(csa);
      /* parse the context that follows the variable name */
      if (csa->token == T_LE)
      {  /* parse upper bound */
         scan_token(csa);
         if (csa->token == T_PLUS || csa->token == T_MINUS)
         {  /* parse signed upper bound */
            s = (csa->token == T_PLUS ? +1.0 : -1.0);
            scan_token(csa);
            if (csa->token == T_NUMBER)
            {  set_upper_bound(csa, j, s * csa->value);
               scan_token(csa);
            }
            else if (the_same(csa->image, "infinity") ||
                     the_same(csa->image, "inf"))
            {  if (s < 0.0)
                  error(csa, "invalid use of '-inf' as upper bound\n");
               set_upper_bound(csa, j, +DBL_MAX);
               scan_token(csa);
            }
            else
               error(csa, "missing upper bound\n");
         }
         else if (csa->token == T_NUMBER)
         {  /* parse unsigned upper bound */
            set_upper_bound(csa, j, csa->value);
            scan_token(csa);
         }
         else
            error(csa, "missing upper bound\n");
      }
      else if (csa->token == T_GE)
      {  /* parse lower bound */
         if (lb_flag)
         {  /* the context '... <= x >= ...' is invalid */
            error(csa, "invalid bound definition\n");
         }
         scan_token(csa);
         if (csa->token == T_PLUS || csa->token == T_MINUS)
         {  /* parse signed lower bound */
            s = (csa->token == T_PLUS ? +1.0 : -1.0);
            scan_token(csa);
            if (csa->token == T_NUMBER)
            {  set_lower_bound(csa, j, s * csa->value);
               scan_token(csa);
            }
            else if (the_same(csa->image, "infinity") ||
                     the_same(csa->image, "inf") == 0)
            {  if (s > 0.0)
                  error(csa, "invalid use of '+inf' as lower bound\n");
               set_lower_bound(csa, j, -DBL_MAX);
               scan_token(csa);
            }
            else
               error(csa, "missing lower bound\n");
         }
         else if (csa->token == T_NUMBER)
         {  /* parse unsigned lower bound */
            set_lower_bound(csa, j, csa->value);
            scan_token(csa);
         }
         else
            error(csa, "missing lower bound\n");
      }
      else if (csa->token == T_EQ)
      {  /* parse fixed value */
         if (lb_flag)
         {  /* the context '... <= x = ...' is invalid */
            error(csa, "invalid bound definition\n");
         }
         scan_token(csa);
         if (csa->token == T_PLUS || csa->token == T_MINUS)
         {  /* parse signed fixed value */
            s = (csa->token == T_PLUS ? +1.0 : -1.0);
            scan_token(csa);
            if (csa->token == T_NUMBER)
            {  set_lower_bound(csa, j, s * csa->value);
               set_upper_bound(csa, j, s * csa->value);
               scan_token(csa);
            }
            else
               error(csa, "missing fixed value\n");
         }
         else if (csa->token == T_NUMBER)
         {  /* parse unsigned fixed value */
            set_lower_bound(csa, j, csa->value);
            set_upper_bound(csa, j, csa->value);
            scan_token(csa);
         }
         else
            error(csa, "missing fixed value\n");
      }
      else if (the_same(csa->image, "free"))
      {  /* parse the keyword 'free' */
         if (lb_flag)
         {  /* the context '... <= x free ...' is invalid */
            error(csa, "invalid bound definition\n");
         }
         set_lower_bound(csa, j, -DBL_MAX);
         set_upper_bound(csa, j, +DBL_MAX);
         scan_token(csa);
      }
      else if (!lb_flag)
      {  /* neither lower nor upper bounds are specified */
         error(csa, "invalid bound definition\n");
      }
      goto loop;
done: return;
}

/***********************************************************************
*  parse_integer - parse general, integer, or binary section
*
*  <variable> ::= <symbolic name>
*  <general> ::= general | generals | gen
*  <integer> ::= integer | integers | int
*  <binary> ::= binary | binaries | bin
*  <section head> ::= <general> <integer> <binary>
*  <additional section> ::= <section head> |
*     <additional section> <variable> */

static void parse_integer(struct csa *csa)
{     int j, binary;
      /* parse the keyword 'general', 'integer', or 'binary' */
      if (csa->token == T_GENERAL)
         binary = 0, scan_token(csa);
      else if (csa->token == T_INTEGER)
         binary = 0, scan_token(csa);
      else if (csa->token == T_BINARY)
         binary = 1, scan_token(csa);
      else
         xassert(csa != csa);
      /* parse list of variables (may be empty) */
      while (csa->token == T_NAME)
      {  /* find the corresponding column */
         j = find_col(csa, csa->image);
         /* change kind of the variable */
         glp_set_col_kind(csa->P, j, GLP_IV);
         /* set bounds for the binary variable */
         if (binary)
#if 0 /* 07/VIII-2013 */
         {  set_lower_bound(csa, j, 0.0);
            set_upper_bound(csa, j, 1.0);
         }
#else
         {  set_lower_bound(csa, j,
               csa->lb[j] == +DBL_MAX ? 0.0 : csa->lb[j]);
            set_upper_bound(csa, j,
               csa->ub[j] == -DBL_MAX ? 1.0 : csa->ub[j]);
         }
#endif
         scan_token(csa);
      }
      return;
}

int glp_read_lp(glp_prob *P, const glp_cpxcp *parm, const char *fname)
{     /* read problem data in CPLEX LP format */
      glp_cpxcp _parm;
      struct csa _csa, *csa = &_csa;
      int ret;
      xprintf("Reading problem data from '%s'...\n", fname);
      if (parm == NULL)
         glp_init_cpxcp(&_parm), parm = &_parm;
      /* check control parameters */
      check_parm("glp_read_lp", parm);
      /* initialize common storage area */
      csa->P = P;
      csa->parm = parm;
      csa->fname = fname;
      csa->fp = NULL;
      if (setjmp(csa->jump))
      {  ret = 1;
         goto done;
      }
      csa->count = 0;
      csa->c = '\n';
      csa->token = T_EOF;
      csa->image[0] = '\0';
      csa->imlen = 0;
      csa->value = 0.0;
      csa->n_max = 100;
      csa->ind = xcalloc(1+csa->n_max, sizeof(int));
      csa->val = xcalloc(1+csa->n_max, sizeof(double));
      csa->flag = xcalloc(1+csa->n_max, sizeof(char));
      memset(&csa->flag[1], 0, csa->n_max * sizeof(char));
      csa->lb = xcalloc(1+csa->n_max, sizeof(double));
      csa->ub = xcalloc(1+csa->n_max, sizeof(double));
#if 1 /* 27/VII-2013 */
      csa->lb_warn = csa->ub_warn = 0;
#endif
      /* erase problem object */
      glp_erase_prob(P);
      glp_create_index(P);
      /* open input CPLEX LP file */
      csa->fp = glp_open(fname, "r");
      if (csa->fp == NULL)
      {  xprintf("Unable to open '%s' - %s\n", fname, get_err_msg());
         ret = 1;
         goto done;
      }
      /* scan very first token */
      scan_token(csa);
      /* parse definition of the objective function */
      if (!(csa->token == T_MINIMIZE || csa->token == T_MAXIMIZE))
         error(csa, "'minimize' or 'maximize' keyword missing\n");
      parse_objective(csa);
      /* parse constraints section */
      if (csa->token != T_SUBJECT_TO)
         error(csa, "constraints section missing\n");
      parse_constraints(csa);
      /* parse optional bounds section */
      if (csa->token == T_BOUNDS) parse_bounds(csa);
      /* parse optional general, integer, and binary sections */
      while (csa->token == T_GENERAL ||
             csa->token == T_INTEGER ||
             csa->token == T_BINARY) parse_integer(csa);
      /* check for the keyword 'end' */
      if (csa->token == T_END)
         scan_token(csa);
      else if (csa->token == T_EOF)
         warning(csa, "keyword 'end' missing\n");
      else
         error(csa, "symbol '%s' in wrong position\n", csa->image);
      /* nothing must follow the keyword 'end' (except comments) */
      if (csa->token != T_EOF)
         error(csa, "extra symbol(s) detected beyond 'end'\n");
      /* set bounds of variables */
      {  int j, type;
         double lb, ub;
         for (j = 1; j <= P->n; j++)
         {  lb = csa->lb[j];
            ub = csa->ub[j];
            if (lb == +DBL_MAX) lb = 0.0;      /* default lb */
            if (ub == -DBL_MAX) ub = +DBL_MAX; /* default ub */
            if (lb == -DBL_MAX && ub == +DBL_MAX)
               type = GLP_FR;
            else if (ub == +DBL_MAX)
               type = GLP_LO;
            else if (lb == -DBL_MAX)
               type = GLP_UP;
            else if (lb != ub)
               type = GLP_DB;
            else
               type = GLP_FX;
            glp_set_col_bnds(csa->P, j, type, lb, ub);
         }
      }
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
      xprintf("%d lines were read\n", csa->count);
      /* problem data has been successfully read */
      glp_delete_index(P);
      glp_sort_matrix(P);
      ret = 0;
done: if (csa->fp != NULL) glp_close(csa->fp);
      xfree(csa->ind);
      xfree(csa->val);
      xfree(csa->flag);
      xfree(csa->lb);
      xfree(csa->ub);
      if (ret != 0) glp_erase_prob(P);
      return ret;
}

/***********************************************************************
*  NAME
*
*  glp_write_lp - write problem data in CPLEX LP format
*
*  SYNOPSIS
*
*  int glp_write_lp(glp_prob *P, const glp_cpxcp *parm, const char
*     *fname);
*
*  DESCRIPTION
*
*  The routine glp_write_lp writes problem data in CPLEX LP format to
*  a text file.
*
*  The parameter parm is a pointer to the structure glp_cpxcp, which
*  specifies control parameters used by the routine. If parm is NULL,
*  the routine uses default settings.
*
*  The character string fname specifies a name of the text file to be
*  written.
*
*  RETURNS
*
*  If the operation was successful, the routine glp_write_lp returns
*  zero. Otherwise, it prints an error message and returns non-zero. */

#define csa csa1

struct csa
{     /* common storage area */
      glp_prob *P;
      /* pointer to problem object */
      const glp_cpxcp *parm;
      /* pointer to control parameters */
};

static int check_name(char *name)
{     /* check if specified name is valid for CPLEX LP format */
      if (*name == '.') return 1;
      if (isdigit((unsigned char)*name)) return 1;
      for (; *name; name++)
      {  if (!isalnum((unsigned char)*name) &&
             strchr(CHAR_SET, (unsigned char)*name) == NULL) return 1;
      }
      return 0; /* name is ok */
}

static void adjust_name(char *name)
{     /* attempt to adjust specified name to make it valid for CPLEX LP
         format */
      for (; *name; name++)
      {  if (*name == ' ')
            *name = '_';
         else if (*name == '-')
            *name = '~';
         else if (*name == '[')
            *name = '(';
         else if (*name == ']')
            *name = ')';
      }
      return;
}

static char *row_name(struct csa *csa, int i, char rname[255+1])
{     /* construct symbolic name of i-th row (constraint) */
      const char *name;
      if (i == 0)
         name = glp_get_obj_name(csa->P);
      else
         name = glp_get_row_name(csa->P, i);
      if (name == NULL) goto fake;
      strcpy(rname, name);
      adjust_name(rname);
      if (check_name(rname)) goto fake;
      return rname;
fake: if (i == 0)
         strcpy(rname, "obj");
      else
         sprintf(rname, "r_%d", i);
      return rname;
}

static char *col_name(struct csa *csa, int j, char cname[255+1])
{     /* construct symbolic name of j-th column (variable) */
      const char *name;
      name = glp_get_col_name(csa->P, j);
      if (name == NULL) goto fake;
      strcpy(cname, name);
      adjust_name(cname);
      if (check_name(cname)) goto fake;
      return cname;
#if 0 /* 18/I-2018 */
fake: sprintf(cname, "x_%d", j);
#else
fake: /* construct fake name depending on column's attributes */
      {  GLPCOL *col = csa->P->col[j];
         if (col->type == GLP_FX)
         {  /* fixed column */
            sprintf(cname, "s_%d", j);
         }
         else if (col->kind == GLP_CV)
         {  /* continuous variable */
            sprintf(cname, "x_%d", j);
         }
         else if (!(col->lb == 0 && col->ub == 1))
         {  /* general (non-binary) integer variable */
            sprintf(cname, "y_%d", j);
         }
         else
         {  /* binary variable */
            sprintf(cname, "z_%d", j);
         }
      }
#endif
      return cname;
}

int glp_write_lp(glp_prob *P, const glp_cpxcp *parm, const char *fname)
{     /* write problem data in CPLEX LP format */
      glp_cpxcp _parm;
      struct csa _csa, *csa = &_csa;
      glp_file *fp;
      GLPROW *row;
      GLPCOL *col;
      GLPAIJ *aij;
      int i, j, len, flag, count, ret;
      char line[1000+1], term[500+1], name[255+1];
      xprintf("Writing problem data to '%s'...\n", fname);
      if (parm == NULL)
         glp_init_cpxcp(&_parm), parm = &_parm;
      /* check control parameters */
      check_parm("glp_write_lp", parm);
      /* initialize common storage area */
      csa->P = P;
      csa->parm = parm;
      /* create output CPLEX LP file */
      fp = glp_open(fname, "w"), count = 0;
      if (fp == NULL)
      {  xprintf("Unable to create '%s' - %s\n", fname, get_err_msg());
         ret = 1;
         goto done;
      }
      /* write problem name */
      xfprintf(fp, "\\* Problem: %s *\\\n",
         P->name == NULL ? "Unknown" : P->name), count++;
      xfprintf(fp, "\n"), count++;
      /* the problem should contain at least one row and one column */
      if (!(P->m > 0 && P->n > 0))
      {  xprintf("Warning: problem has no rows/columns\n");
         xfprintf(fp, "\\* WARNING: PROBLEM HAS NO ROWS/COLUMNS *\\\n"),
            count++;
         xfprintf(fp, "\n"), count++;
         goto skip;
      }
      /* write the objective function definition */
      if (P->dir == GLP_MIN)
         xfprintf(fp, "Minimize\n"), count++;
      else if (P->dir == GLP_MAX)
         xfprintf(fp, "Maximize\n"), count++;
      else
         xassert(P != P);
      row_name(csa, 0, name);
      sprintf(line, " %s:", name);
      len = 0;
      for (j = 1; j <= P->n; j++)
      {  col = P->col[j];
         if (col->coef != 0.0 || col->ptr == NULL)
         {  len++;
            col_name(csa, j, name);
            if (col->coef == 0.0)
               sprintf(term, " + 0 %s", name); /* empty column */
            else if (col->coef == +1.0)
               sprintf(term, " + %s", name);
            else if (col->coef == -1.0)
               sprintf(term, " - %s", name);
            else if (col->coef > 0.0)
               sprintf(term, " + %.*g %s", DBL_DIG, +col->coef, name);
            else
               sprintf(term, " - %.*g %s", DBL_DIG, -col->coef, name);
            if (strlen(line) + strlen(term) > 72)
               xfprintf(fp, "%s\n", line), line[0] = '\0', count++;
            strcat(line, term);
         }
      }
      if (len == 0)
      {  /* empty objective */
         sprintf(term, " 0 %s", col_name(csa, 1, name));
         strcat(line, term);
      }
      xfprintf(fp, "%s\n", line), count++;
      if (P->c0 != 0.0)
         xfprintf(fp, "\\* constant term = %.*g *\\\n", DBL_DIG, P->c0),
            count++;
      xfprintf(fp, "\n"), count++;
      /* write the constraints section */
      xfprintf(fp, "Subject To\n"), count++;
      for (i = 1; i <= P->m; i++)
      {  row = P->row[i];
         if (row->type == GLP_FR) continue; /* skip free row */
         row_name(csa, i, name);
         sprintf(line, " %s:", name);
         /* linear form */
         for (aij = row->ptr; aij != NULL; aij = aij->r_next)
         {  col_name(csa, aij->col->j, name);
            if (aij->val == +1.0)
               sprintf(term, " + %s", name);
            else if (aij->val == -1.0)
               sprintf(term, " - %s", name);
            else if (aij->val > 0.0)
               sprintf(term, " + %.*g %s", DBL_DIG, +aij->val, name);
            else
               sprintf(term, " - %.*g %s", DBL_DIG, -aij->val, name);
            if (strlen(line) + strlen(term) > 72)
               xfprintf(fp, "%s\n", line), line[0] = '\0', count++;
            strcat(line, term);
         }
         if (row->type == GLP_DB)
         {  /* double-bounded (ranged) constraint */
            sprintf(term, " - ~r_%d", i);
            if (strlen(line) + strlen(term) > 72)
               xfprintf(fp, "%s\n", line), line[0] = '\0', count++;
            strcat(line, term);
         }
         else if (row->ptr == NULL)
         {  /* empty constraint */
            sprintf(term, " 0 %s", col_name(csa, 1, name));
            strcat(line, term);
         }
         /* right hand-side */
         if (row->type == GLP_LO)
            sprintf(term, " >= %.*g", DBL_DIG, row->lb);
         else if (row->type == GLP_UP)
            sprintf(term, " <= %.*g", DBL_DIG, row->ub);
         else if (row->type == GLP_DB || row->type == GLP_FX)
            sprintf(term, " = %.*g", DBL_DIG, row->lb);
         else
            xassert(row != row);
         if (strlen(line) + strlen(term) > 72)
            xfprintf(fp, "%s\n", line), line[0] = '\0', count++;
         strcat(line, term);
         xfprintf(fp, "%s\n", line), count++;
      }
      xfprintf(fp, "\n"), count++;
      /* write the bounds section */
      flag = 0;
      for (i = 1; i <= P->m; i++)
      {  row = P->row[i];
         if (row->type != GLP_DB) continue;
         if (!flag)
            xfprintf(fp, "Bounds\n"), flag = 1, count++;
         xfprintf(fp, " 0 <= ~r_%d <= %.*g\n",
            i, DBL_DIG, row->ub - row->lb), count++;
      }
      for (j = 1; j <= P->n; j++)
      {  col = P->col[j];
         if (col->type == GLP_LO && col->lb == 0.0) continue;
         if (!flag)
            xfprintf(fp, "Bounds\n"), flag = 1, count++;
         col_name(csa, j, name);
         if (col->type == GLP_FR)
            xfprintf(fp, " %s free\n", name), count++;
         else if (col->type == GLP_LO)
            xfprintf(fp, " %s >= %.*g\n",
               name, DBL_DIG, col->lb), count++;
         else if (col->type == GLP_UP)
            xfprintf(fp, " -Inf <= %s <= %.*g\n",
               name, DBL_DIG, col->ub), count++;
         else if (col->type == GLP_DB)
            xfprintf(fp, " %.*g <= %s <= %.*g\n",
               DBL_DIG, col->lb, name, DBL_DIG, col->ub), count++;
         else if (col->type == GLP_FX)
            xfprintf(fp, " %s = %.*g\n",
               name, DBL_DIG, col->lb), count++;
         else
            xassert(col != col);
      }
      if (flag) xfprintf(fp, "\n"), count++;
      /* write the integer section */
      flag = 0;
      for (j = 1; j <= P->n; j++)
      {  col = P->col[j];
         if (col->kind == GLP_CV) continue;
         xassert(col->kind == GLP_IV);
         if (!flag)
            xfprintf(fp, "Generals\n"), flag = 1, count++;
         xfprintf(fp, " %s\n", col_name(csa, j, name)), count++;
      }
      if (flag) xfprintf(fp, "\n"), count++;
skip: /* write the end keyword */
      xfprintf(fp, "End\n"), count++;
#if 0 /* FIXME */
      xfflush(fp);
#endif
      if (glp_ioerr(fp))
      {  xprintf("Write error on '%s' - %s\n", fname, get_err_msg());
         ret = 1;
         goto done;
      }
      /* problem data has been successfully written */
      xprintf("%d lines were written\n", count);
      ret = 0;
done: if (fp != NULL) glp_close(fp);
      return ret;
}

/* eof */
