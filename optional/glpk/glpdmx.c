/* glpdmx.c (reading/writing data in DIMACS format) */

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

#define _GLPSTD_STDIO
#include "glpapi.h"

struct csa
{     /* common storage area */
      jmp_buf jump;
      /* label for go to in case of error */
      const char *fname;
      /* name of input text file */
      XFILE *fp;
      /* stream assigned to input text file */
      int count;
      /* line count */
      int c;
      /* current character */
      char field[255+1];
      /* data field */
      int empty;
      /* warning 'empty line ignored' was printed */
      int nonint;
      /* warning 'non-integer data detected' was printed */
};

static void error(struct csa *csa, const char *fmt, ...)
{     /* print error message and terminate processing */
      va_list arg;
      xprintf("%s:%d: error: ", csa->fname, csa->count);
      va_start(arg, fmt);
      xvprintf(fmt, arg);
      va_end(arg);
      xprintf("\n");
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
      xprintf("\n");
      return;
}

static void read_char(struct csa *csa)
{     /* read character from input text file */
      int c;
      if (csa->c == '\n') csa->count++;
      c = xfgetc(csa->fp);
      if (c < 0)
      {  if (xferror(csa->fp))
            error(csa, "read error - %s", xerrmsg());
         else if (csa->c == '\n')
            error(csa, "unexpected end of file");
         else
         {  warning(csa, "missing final end of line");
            c = '\n';
         }
      }
      else if (c == '\n')
         ;
      else if (isspace(c))
         c = ' ';
      else if (iscntrl(c))
         error(csa, "invalid control character 0x%02X", c);
      csa->c = c;
      return;
}

static void read_designator(struct csa *csa)
{     /* read one-character line designator */
      xassert(csa->c == '\n');
      read_char(csa);
      for (;;)
      {  /* skip preceding white-space characters */
         while (csa->c == ' ')
            read_char(csa);
         if (csa->c == '\n')
         {  /* ignore empty line */
            if (!csa->empty)
            {  warning(csa, "empty line ignored");
               csa->empty = 1;
            }
            read_char(csa);
         }
         else if (csa->c == 'c')
         {  /* skip comment line */
            while (csa->c != '\n')
               read_char(csa);
            read_char(csa);
         }
         else
         {  /* hmm... looks like a line designator */
            csa->field[0] = (char)csa->c, csa->field[1] = '\0';
            /* check that it is followed by a white-space character */
            read_char(csa);
            if (!(csa->c == ' ' || csa->c == '\n'))
               error(csa, "line designator missing or invalid");
            break;
         }
      }
      return;
}

static void read_field(struct csa *csa)
{     /* read data field */
      int len = 0;
      /* skip preceding white-space characters */
      while (csa->c == ' ')
         read_char(csa);
      /* scan data field */
      if (csa->c == '\n')
         error(csa, "unexpected end of line");
      while (!(csa->c == ' ' || csa->c == '\n'))
      {  if (len == sizeof(csa->field)-1)
            error(csa, "data field `%.15s...' too long", csa->field);
         csa->field[len++] = (char)csa->c;
         read_char(csa);
      }
      csa->field[len] = '\0';
      return;
}

static void end_of_line(struct csa *csa)
{     /* skip white-space characters until end of line */
      while (csa->c == ' ')
         read_char(csa);
      if (csa->c != '\n')
         error(csa, "too many data fields specified");
      return;
}

static void check_int(struct csa *csa, double num)
{     /* print a warning if non-integer data are detected */
      if (!csa->nonint && num != floor(num))
      {  warning(csa, "non-integer data detected");
         csa->nonint = 1;
      }
      return;
}

/***********************************************************************
*  NAME
*
*  glp_read_mincost - read min-cost flow problem data in DIMACS format
*
*  SYNOPSIS
*
*  int glp_read_mincost(glp_graph *G, int v_rhs, int a_low, int a_cap,
*     int a_cost, const char *fname);
*
*  DESCRIPTION
*
*  The routine glp_read_mincost reads minimum cost flow problem data in
*  DIMACS format from a text file.
*
*  RETURNS
*
*  If the operation was successful, the routine returns zero. Otherwise
*  it prints an error message and returns non-zero. */

int glp_read_mincost(glp_graph *G, int v_rhs, int a_low, int a_cap,
      int a_cost, const char *fname)
{     struct csa _csa, *csa = &_csa;
      glp_vertex *v;
      glp_arc *a;
      int i, j, k, nv, na, ret = 0;
      double rhs, low, cap, cost;
      char *flag = NULL;
      if (v_rhs >= 0 && v_rhs > G->v_size - (int)sizeof(double))
         xerror("glp_read_mincost: v_rhs = %d; invalid offset\n",
            v_rhs);
      if (a_low >= 0 && a_low > G->a_size - (int)sizeof(double))
         xerror("glp_read_mincost: a_low = %d; invalid offset\n",
            a_low);
      if (a_cap >= 0 && a_cap > G->a_size - (int)sizeof(double))
         xerror("glp_read_mincost: a_cap = %d; invalid offset\n",
            a_cap);
      if (a_cost >= 0 && a_cost > G->a_size - (int)sizeof(double))
         xerror("glp_read_mincost: a_cost = %d; invalid offset\n",
            a_cost);
      glp_erase_graph(G, G->v_size, G->a_size);
      if (setjmp(csa->jump))
      {  ret = 1;
         goto done;
      }
      csa->fname = fname;
      csa->fp = NULL;
      csa->count = 0;
      csa->c = '\n';
      csa->field[0] = '\0';
      csa->empty = csa->nonint = 0;
      xprintf("Reading min-cost flow problem data from `%s'...\n",
         fname);
      csa->fp = xfopen(fname, "r");
      if (csa->fp == NULL)
      {  xprintf("Unable to open `%s' - %s\n", fname, xerrmsg());
         longjmp(csa->jump, 1);
      }
      /* read problem line */
      read_designator(csa);
      if (strcmp(csa->field, "p") != 0)
         error(csa, "problem line missing or invalid");
      read_field(csa);
      if (strcmp(csa->field, "min") != 0)
         error(csa, "wrong problem designator; `min' expected");
      read_field(csa);
      if (!(str2int(csa->field, &nv) == 0 && nv >= 0))
         error(csa, "number of nodes missing or invalid");
      read_field(csa);
      if (!(str2int(csa->field, &na) == 0 && na >= 0))
         error(csa, "number of arcs missing or invalid");
      xprintf("Flow network has %d node%s and %d arc%s\n",
         nv, nv == 1 ? "" : "s", na, na == 1 ? "" : "s");
      if (nv > 0) glp_add_vertices(G, nv);
      end_of_line(csa);
      /* read node descriptor lines */
      flag = xcalloc(1+nv, sizeof(char));
      memset(&flag[1], 0, nv * sizeof(char));
      if (v_rhs >= 0)
      {  rhs = 0.0;
         for (i = 1; i <= nv; i++)
         {  v = G->v[i];
            memcpy((char *)v->data + v_rhs, &rhs, sizeof(double));
         }
      }
      for (;;)
      {  read_designator(csa);
         if (strcmp(csa->field, "n") != 0) break;
         read_field(csa);
         if (str2int(csa->field, &i) != 0)
            error(csa, "node number missing or invalid");
         if (!(1 <= i && i <= nv))
            error(csa, "node number %d out of range", i);
         if (flag[i])
            error(csa, "duplicate descriptor of node %d", i);
         read_field(csa);
         if (str2num(csa->field, &rhs) != 0)
            error(csa, "node supply/demand missing or invalid");
         check_int(csa, rhs);
         if (v_rhs >= 0)
         {  v = G->v[i];
            memcpy((char *)v->data + v_rhs, &rhs, sizeof(double));
         }
         flag[i] = 1;
         end_of_line(csa);
      }
      xfree(flag), flag = NULL;
      /* read arc descriptor lines */
      for (k = 1; k <= na; k++)
      {  if (k > 1) read_designator(csa);
         if (strcmp(csa->field, "a") != 0)
            error(csa, "wrong line designator; `a' expected");
         read_field(csa);
         if (str2int(csa->field, &i) != 0)
            error(csa, "starting node number missing or invalid");
         if (!(1 <= i && i <= nv))
            error(csa, "starting node number %d out of range", i);
         read_field(csa);
         if (str2int(csa->field, &j) != 0)
            error(csa, "ending node number missing or invalid");
         if (!(1 <= j && j <= nv))
            error(csa, "ending node number %d out of range", j);
         read_field(csa);
         if (!(str2num(csa->field, &low) == 0 && low >= 0.0))
            error(csa, "lower bound of arc flow missing or invalid");
         check_int(csa, low);
         read_field(csa);
         if (!(str2num(csa->field, &cap) == 0 && cap >= low))
            error(csa, "upper bound of arc flow missing or invalid");
         check_int(csa, cap);
         read_field(csa);
         if (str2num(csa->field, &cost) != 0)
            error(csa, "per-unit cost of arc flow missing or invalid");
         check_int(csa, cost);
         a = glp_add_arc(G, i, j);
         if (a_low >= 0)
            memcpy((char *)a->data + a_low, &low, sizeof(double));
         if (a_cap >= 0)
            memcpy((char *)a->data + a_cap, &cap, sizeof(double));
         if (a_cost >= 0)
            memcpy((char *)a->data + a_cost, &cost, sizeof(double));
         end_of_line(csa);
      }
      xprintf("%d lines were read\n", csa->count);
done: if (ret) glp_erase_graph(G, G->v_size, G->a_size);
      if (csa->fp != NULL) xfclose(csa->fp);
      if (flag != NULL) xfree(flag);
      return ret;
}

/***********************************************************************
*  NAME
*
*  glp_write_mincost - write min-cost flow problem data in DIMACS format
*
*  SYNOPSIS
*
*  int glp_write_mincost(glp_graph *G, int v_rhs, int a_low, int a_cap,
*     int a_cost, const char *fname);
*
*  DESCRIPTION
*
*  The routine glp_write_mincost writes minimum cost flow problem data
*  in DIMACS format to a text file.
*
*  RETURNS
*
*  If the operation was successful, the routine returns zero. Otherwise
*  it prints an error message and returns non-zero. */

int glp_write_mincost(glp_graph *G, int v_rhs, int a_low, int a_cap,
      int a_cost, const char *fname)
{     XFILE *fp;
      glp_vertex *v;
      glp_arc *a;
      int i, count = 0, ret;
      double rhs, low, cap, cost;
      if (v_rhs >= 0 && v_rhs > G->v_size - (int)sizeof(double))
         xerror("glp_write_mincost: v_rhs = %d; invalid offset\n",
            v_rhs);
      if (a_low >= 0 && a_low > G->a_size - (int)sizeof(double))
         xerror("glp_write_mincost: a_low = %d; invalid offset\n",
            a_low);
      if (a_cap >= 0 && a_cap > G->a_size - (int)sizeof(double))
         xerror("glp_write_mincost: a_cap = %d; invalid offset\n",
            a_cap);
      if (a_cost >= 0 && a_cost > G->a_size - (int)sizeof(double))
         xerror("glp_write_mincost: a_cost = %d; invalid offset\n",
            a_cost);
      xprintf("Writing min-cost flow problem data to `%s'...\n",
         fname);
      fp = xfopen(fname, "w");
      if (fp == NULL)
      {  xprintf("Unable to create `%s' - %s\n", fname, xerrmsg());
         ret = 1;
         goto done;
      }
      xfprintf(fp, "c %s\n",
         G->name == NULL ? "unknown" : G->name), count++;
      xfprintf(fp, "p min %d %d\n", G->nv, G->na), count++;
      if (v_rhs >= 0)
      {  for (i = 1; i <= G->nv; i++)
         {  v = G->v[i];
            memcpy(&rhs, (char *)v->data + v_rhs, sizeof(double));
            if (rhs != 0.0)
               xfprintf(fp, "n %d %.*g\n", i, DBL_DIG, rhs), count++;
         }
      }
      for (i = 1; i <= G->nv; i++)
      {  v = G->v[i];
         for (a = v->out; a != NULL; a = a->t_next)
         {  if (a_low >= 0)
               memcpy(&low, (char *)a->data + a_low, sizeof(double));
            else
               low = 0.0;
            if (a_cap >= 0)
               memcpy(&cap, (char *)a->data + a_cap, sizeof(double));
            else
               cap = 1.0;
            if (a_cost >= 0)
               memcpy(&cost, (char *)a->data + a_cost, sizeof(double));
            else
               cost = 0.0;
            xfprintf(fp, "a %d %d %.*g %.*g %.*g\n",
               a->tail->i, a->head->i, DBL_DIG, low, DBL_DIG, cap,
               DBL_DIG, cost), count++;
         }
      }
      xfprintf(fp, "c eof\n"), count++;
      xfflush(fp);
      if (xferror(fp))
      {  xprintf("Write error on `%s' - %s\n", fname, xerrmsg());
         ret = 1;
         goto done;
      }
      xprintf("%d lines were written\n", count);
      ret = 0;
done: if (fp != NULL) xfclose(fp);
      return ret;
}

/***********************************************************************
*  NAME
*
*  glp_read_maxflow - read maximum flow problem data in DIMACS format
*
*  SYNOPSIS
*
*  int glp_read_maxflow(glp_graph *G, int *s, int *t, int a_cap,
*     const char *fname);
*
*  DESCRIPTION
*
*  The routine glp_read_maxflow reads maximum flow problem data in
*  DIMACS format from a text file.
*
*  RETURNS
*
*  If the operation was successful, the routine returns zero. Otherwise
*  it prints an error message and returns non-zero. */

int glp_read_maxflow(glp_graph *G, int *_s, int *_t, int a_cap,
      const char *fname)
{     struct csa _csa, *csa = &_csa;
      glp_arc *a;
      int i, j, k, s, t, nv, na, ret = 0;
      double cap;
      if (a_cap >= 0 && a_cap > G->a_size - (int)sizeof(double))
         xerror("glp_read_maxflow: a_cap = %d; invalid offset\n",
            a_cap);
      glp_erase_graph(G, G->v_size, G->a_size);
      if (setjmp(csa->jump))
      {  ret = 1;
         goto done;
      }
      csa->fname = fname;
      csa->fp = NULL;
      csa->count = 0;
      csa->c = '\n';
      csa->field[0] = '\0';
      csa->empty = csa->nonint = 0;
      xprintf("Reading maximum flow problem data from `%s'...\n",
         fname);
      csa->fp = xfopen(fname, "r");
      if (csa->fp == NULL)
      {  xprintf("Unable to open `%s' - %s\n", fname, xerrmsg());
         longjmp(csa->jump, 1);
      }
      /* read problem line */
      read_designator(csa);
      if (strcmp(csa->field, "p") != 0)
         error(csa, "problem line missing or invalid");
      read_field(csa);
      if (strcmp(csa->field, "max") != 0)
         error(csa, "wrong problem designator; `max' expected");
      read_field(csa);
      if (!(str2int(csa->field, &nv) == 0 && nv >= 2))
         error(csa, "number of nodes missing or invalid");
      read_field(csa);
      if (!(str2int(csa->field, &na) == 0 && na >= 0))
         error(csa, "number of arcs missing or invalid");
      xprintf("Flow network has %d node%s and %d arc%s\n",
         nv, nv == 1 ? "" : "s", na, na == 1 ? "" : "s");
      if (nv > 0) glp_add_vertices(G, nv);
      end_of_line(csa);
      /* read node descriptor lines */
      s = t = 0;
      for (;;)
      {  read_designator(csa);
         if (strcmp(csa->field, "n") != 0) break;
         read_field(csa);
         if (str2int(csa->field, &i) != 0)
            error(csa, "node number missing or invalid");
         if (!(1 <= i && i <= nv))
            error(csa, "node number %d out of range", i);
         read_field(csa);
         if (strcmp(csa->field, "s") == 0)
         {  if (s > 0)
               error(csa, "only one source node allowed");
            s = i;
         }
         else if (strcmp(csa->field, "t") == 0)
         {  if (t > 0)
               error(csa, "only one sink node allowed");
            t = i;
         }
         else
            error(csa, "wrong node designator; `s' or `t' expected");
         if (s > 0 && s == t)
            error(csa, "source and sink nodes must be distinct");
         end_of_line(csa);
      }
      if (s == 0)
         error(csa, "source node descriptor missing\n");
      if (t == 0)
         error(csa, "sink node descriptor missing\n");
      if (_s != NULL) *_s = s;
      if (_t != NULL) *_t = t;
      /* read arc descriptor lines */
      for (k = 1; k <= na; k++)
      {  if (k > 1) read_designator(csa);
         if (strcmp(csa->field, "a") != 0)
            error(csa, "wrong line designator; `a' expected");
         read_field(csa);
         if (str2int(csa->field, &i) != 0)
            error(csa, "starting node number missing or invalid");
         if (!(1 <= i && i <= nv))
            error(csa, "starting node number %d out of range", i);
         read_field(csa);
         if (str2int(csa->field, &j) != 0)
            error(csa, "ending node number missing or invalid");
         if (!(1 <= j && j <= nv))
            error(csa, "ending node number %d out of range", j);
         read_field(csa);
         if (!(str2num(csa->field, &cap) == 0 && cap >= 0.0))
            error(csa, "arc capacity missing or invalid");
         check_int(csa, cap);
         a = glp_add_arc(G, i, j);
         if (a_cap >= 0)
            memcpy((char *)a->data + a_cap, &cap, sizeof(double));
         end_of_line(csa);
      }
      xprintf("%d lines were read\n", csa->count);
done: if (ret) glp_erase_graph(G, G->v_size, G->a_size);
      if (csa->fp != NULL) xfclose(csa->fp);
      return ret;
}

/***********************************************************************
*  NAME
*
*  glp_write_maxflow - write maximum flow problem data in DIMACS format
*
*  SYNOPSIS
*
*  int glp_write_maxflow(glp_graph *G, int s, int t, int a_cap,
*     const char *fname);
*
*  DESCRIPTION
*
*  The routine glp_write_maxflow writes maximum flow problem data in
*  DIMACS format to a text file.
*
*  RETURNS
*
*  If the operation was successful, the routine returns zero. Otherwise
*  it prints an error message and returns non-zero. */

int glp_write_maxflow(glp_graph *G, int s, int t, int a_cap,
      const char *fname)
{     XFILE *fp;
      glp_vertex *v;
      glp_arc *a;
      int i, count = 0, ret;
      double cap;
      if (!(1 <= s && s <= G->nv))
         xerror("glp_write_maxflow: s = %d; source node number out of r"
            "ange\n", s);
      if (!(1 <= t && t <= G->nv))
         xerror("glp_write_maxflow: t = %d: sink node number out of ran"
            "ge\n", t);
      if (a_cap >= 0 && a_cap > G->a_size - (int)sizeof(double))
         xerror("glp_write_mincost: a_cap = %d; invalid offset\n",
            a_cap);
      xprintf("Writing maximum flow problem data to `%s'...\n",
         fname);
      fp = xfopen(fname, "w");
      if (fp == NULL)
      {  xprintf("Unable to create `%s' - %s\n", fname, xerrmsg());
         ret = 1;
         goto done;
      }
      xfprintf(fp, "c %s\n",
         G->name == NULL ? "unknown" : G->name), count++;
      xfprintf(fp, "p max %d %d\n", G->nv, G->na), count++;
      xfprintf(fp, "n %d s\n", s), count++;
      xfprintf(fp, "n %d t\n", t), count++;
      for (i = 1; i <= G->nv; i++)
      {  v = G->v[i];
         for (a = v->out; a != NULL; a = a->t_next)
         {  if (a_cap >= 0)
               memcpy(&cap, (char *)a->data + a_cap, sizeof(double));
            else
               cap = 1.0;
            xfprintf(fp, "a %d %d %.*g\n",
               a->tail->i, a->head->i, DBL_DIG, cap), count++;
         }
      }
      xfprintf(fp, "c eof\n"), count++;
      xfflush(fp);
      if (xferror(fp))
      {  xprintf("Write error on `%s' - %s\n", fname, xerrmsg());
         ret = 1;
         goto done;
      }
      xprintf("%d lines were written\n", count);
      ret = 0;
done: if (fp != NULL) xfclose(fp);
      return ret;
}

/***********************************************************************
*  NAME
*
*  glp_read_asnprob - read assignment problem data in DIMACS format
*
*  SYNOPSIS
*
*  int glp_read_asnprob(glp_graph *G, int v_set, int a_cost,
*     const char *fname);
*
*  DESCRIPTION
*
*  The routine glp_read_asnprob reads assignment problem data in DIMACS
*  format from a text file.
*
*  RETURNS
*
*  If the operation was successful, the routine returns zero. Otherwise
*  it prints an error message and returns non-zero. */

int glp_read_asnprob(glp_graph *G, int v_set, int a_cost, const char
      *fname)
{     struct csa _csa, *csa = &_csa;
      glp_vertex *v;
      glp_arc *a;
      int nv, na, n1, i, j, k, ret = 0;
      double cost;
      char *flag = NULL;
      if (v_set >= 0 && v_set > G->v_size - (int)sizeof(int))
         xerror("glp_read_asnprob: v_set = %d; invalid offset\n",
            v_set);
      if (a_cost >= 0 && a_cost > G->a_size - (int)sizeof(double))
         xerror("glp_read_asnprob: a_cost = %d; invalid offset\n",
            a_cost);
      glp_erase_graph(G, G->v_size, G->a_size);
      if (setjmp(csa->jump))
      {  ret = 1;
         goto done;
      }
      csa->fname = fname;
      csa->fp = NULL;
      csa->count = 0;
      csa->c = '\n';
      csa->field[0] = '\0';
      csa->empty = csa->nonint = 0;
      xprintf("Reading assignment problem data from `%s'...\n", fname);
      csa->fp = xfopen(fname, "r");
      if (csa->fp == NULL)
      {  xprintf("Unable to open `%s' - %s\n", fname, xerrmsg());
         longjmp(csa->jump, 1);
      }
      /* read problem line */
      read_designator(csa);
      if (strcmp(csa->field, "p") != 0)
         error(csa, "problem line missing or invalid");
      read_field(csa);
      if (strcmp(csa->field, "asn") != 0)
         error(csa, "wrong problem designator; `asn' expected");
      read_field(csa);
      if (!(str2int(csa->field, &nv) == 0 && nv >= 0))
         error(csa, "number of nodes missing or invalid");
      read_field(csa);
      if (!(str2int(csa->field, &na) == 0 && na >= 0))
         error(csa, "number of arcs missing or invalid");
      if (nv > 0) glp_add_vertices(G, nv);
      end_of_line(csa);
      /* read node descriptor lines */
      flag = xcalloc(1+nv, sizeof(char));
      memset(&flag[1], 0, nv * sizeof(char));
      n1 = 0;
      for (;;)
      {  read_designator(csa);
         if (strcmp(csa->field, "n") != 0) break;
         read_field(csa);
         if (str2int(csa->field, &i) != 0)
            error(csa, "node number missing or invalid");
         if (!(1 <= i && i <= nv))
            error(csa, "node number %d out of range", i);
         if (flag[i])
            error(csa, "duplicate descriptor of node %d", i);
         flag[i] = 1, n1++;
         end_of_line(csa);
      }
      xprintf(
         "Assignment problem has %d + %d = %d node%s and %d arc%s\n",
         n1, nv - n1, nv, nv == 1 ? "" : "s", na, na == 1 ? "" : "s");
      if (v_set >= 0)
      {  for (i = 1; i <= nv; i++)
         {  v = G->v[i];
            k = (flag[i] ? 0 : 1);
            memcpy((char *)v->data + v_set, &k, sizeof(int));
         }
      }
      /* read arc descriptor lines */
      for (k = 1; k <= na; k++)
      {  if (k > 1) read_designator(csa);
         if (strcmp(csa->field, "a") != 0)
            error(csa, "wrong line designator; `a' expected");
         read_field(csa);
         if (str2int(csa->field, &i) != 0)
            error(csa, "starting node number missing or invalid");
         if (!(1 <= i && i <= nv))
            error(csa, "starting node number %d out of range", i);
         if (!flag[i])
            error(csa, "node %d cannot be a starting node", i);
         read_field(csa);
         if (str2int(csa->field, &j) != 0)
            error(csa, "ending node number missing or invalid");
         if (!(1 <= j && j <= nv))
            error(csa, "ending node number %d out of range", j);
         if (flag[j])
            error(csa, "node %d cannot be an ending node", j);
         read_field(csa);
         if (str2num(csa->field, &cost) != 0)
            error(csa, "arc cost missing or invalid");
         check_int(csa, cost);
         a = glp_add_arc(G, i, j);
         if (a_cost >= 0)
            memcpy((char *)a->data + a_cost, &cost, sizeof(double));
         end_of_line(csa);
      }
      xprintf("%d lines were read\n", csa->count);
done: if (ret) glp_erase_graph(G, G->v_size, G->a_size);
      if (csa->fp != NULL) xfclose(csa->fp);
      if (flag != NULL) xfree(flag);
      return ret;
}

/***********************************************************************
*  NAME
*
*  glp_write_asnprob - write assignment problem data in DIMACS format
*
*  SYNOPSIS
*
*  int glp_write_asnprob(glp_graph *G, int v_set, int a_cost,
*     const char *fname);
*
*  DESCRIPTION
*
*  The routine glp_write_asnprob writes assignment problem data in
*  DIMACS format to a text file.
*
*  RETURNS
*
*  If the operation was successful, the routine returns zero. Otherwise
*  it prints an error message and returns non-zero. */

int glp_write_asnprob(glp_graph *G, int v_set, int a_cost, const char
      *fname)
{     XFILE *fp;
      glp_vertex *v;
      glp_arc *a;
      int i, k, count = 0, ret;
      double cost;
      if (v_set >= 0 && v_set > G->v_size - (int)sizeof(int))
         xerror("glp_write_asnprob: v_set = %d; invalid offset\n",
            v_set);
      if (a_cost >= 0 && a_cost > G->a_size - (int)sizeof(double))
         xerror("glp_write_asnprob: a_cost = %d; invalid offset\n",
            a_cost);
      xprintf("Writing assignment problem data to `%s'...\n", fname);
      fp = xfopen(fname, "w");
      if (fp == NULL)
      {  xprintf("Unable to create `%s' - %s\n", fname, xerrmsg());
         ret = 1;
         goto done;
      }
      xfprintf(fp, "c %s\n",
         G->name == NULL ? "unknown" : G->name), count++;
      xfprintf(fp, "p asn %d %d\n", G->nv, G->na), count++;
      for (i = 1; i <= G->nv; i++)
      {  v = G->v[i];
         if (v_set >= 0)
            memcpy(&k, (char *)v->data + v_set, sizeof(int));
         else
            k = (v->out != NULL ? 0 : 1);
         if (k == 0)
            xfprintf(fp, "n %d\n", i), count++;
      }
      for (i = 1; i <= G->nv; i++)
      {  v = G->v[i];
         for (a = v->out; a != NULL; a = a->t_next)
         {  if (a_cost >= 0)
               memcpy(&cost, (char *)a->data + a_cost, sizeof(double));
            else
               cost = 1.0;
            xfprintf(fp, "a %d %d %.*g\n",
               a->tail->i, a->head->i, DBL_DIG, cost), count++;
         }
      }
      xfprintf(fp, "c eof\n"), count++;
      xfflush(fp);
      if (xferror(fp))
      {  xprintf("Write error on `%s' - %s\n", fname, xerrmsg());
         ret = 1;
         goto done;
      }
      xprintf("%d lines were written\n", count);
      ret = 0;
done: if (fp != NULL) xfclose(fp);
      return ret;
}

/***********************************************************************
*  NAME
*
*  glp_read_ccdata - read graph in DIMACS clique/coloring format
*
*  SYNOPSIS
*
*  int glp_read_ccdata(glp_graph *G, int v_wgt, const char *fname);
*
*  DESCRIPTION
*
*  The routine glp_read_ccdata reads an (undirected) graph in DIMACS
*  clique/coloring format from a text file.
*
*  RETURNS
*
*  If the operation was successful, the routine returns zero. Otherwise
*  it prints an error message and returns non-zero. */

int glp_read_ccdata(glp_graph *G, int v_wgt, const char *fname)
{     struct csa _csa, *csa = &_csa;
      glp_vertex *v;
      int i, j, k, nv, ne, ret = 0;
      double w;
      char *flag = NULL;
      if (v_wgt >= 0 && v_wgt > G->v_size - (int)sizeof(double))
         xerror("glp_read_ccdata: v_wgt = %d; invalid offset\n",
            v_wgt);
      glp_erase_graph(G, G->v_size, G->a_size);
      if (setjmp(csa->jump))
      {  ret = 1;
         goto done;
      }
      csa->fname = fname;
      csa->fp = NULL;
      csa->count = 0;
      csa->c = '\n';
      csa->field[0] = '\0';
      csa->empty = csa->nonint = 0;
      xprintf("Reading graph from `%s'...\n", fname);
      csa->fp = xfopen(fname, "r");
      if (csa->fp == NULL)
      {  xprintf("Unable to open `%s' - %s\n", fname, xerrmsg());
         longjmp(csa->jump, 1);
      }
      /* read problem line */
      read_designator(csa);
      if (strcmp(csa->field, "p") != 0)
         error(csa, "problem line missing or invalid");
      read_field(csa);
      if (strcmp(csa->field, "edge") != 0)
         error(csa, "wrong problem designator; `edge' expected");
      read_field(csa);
      if (!(str2int(csa->field, &nv) == 0 && nv >= 0))
         error(csa, "number of vertices missing or invalid");
      read_field(csa);
      if (!(str2int(csa->field, &ne) == 0 && ne >= 0))
         error(csa, "number of edges missing or invalid");
      xprintf("Graph has %d vert%s and %d edge%s\n",
         nv, nv == 1 ? "ex" : "ices", ne, ne == 1 ? "" : "s");
      if (nv > 0) glp_add_vertices(G, nv);
      end_of_line(csa);
      /* read node descriptor lines */
      flag = xcalloc(1+nv, sizeof(char));
      memset(&flag[1], 0, nv * sizeof(char));
      if (v_wgt >= 0)
      {  w = 1.0;
         for (i = 1; i <= nv; i++)
         {  v = G->v[i];
            memcpy((char *)v->data + v_wgt, &w, sizeof(double));
         }
      }
      for (;;)
      {  read_designator(csa);
         if (strcmp(csa->field, "n") != 0) break;
         read_field(csa);
         if (str2int(csa->field, &i) != 0)
            error(csa, "vertex number missing or invalid");
         if (!(1 <= i && i <= nv))
            error(csa, "vertex number %d out of range", i);
         if (flag[i])
            error(csa, "duplicate descriptor of vertex %d", i);
         read_field(csa);
         if (str2num(csa->field, &w) != 0)
            error(csa, "vertex weight missing or invalid");
         check_int(csa, w);
         if (v_wgt >= 0)
         {  v = G->v[i];
            memcpy((char *)v->data + v_wgt, &w, sizeof(double));
         }
         flag[i] = 1;
         end_of_line(csa);
      }
      xfree(flag), flag = NULL;
      /* read edge descriptor lines */
      for (k = 1; k <= ne; k++)
      {  if (k > 1) read_designator(csa);
         if (strcmp(csa->field, "e") != 0)
            error(csa, "wrong line designator; `e' expected");
         read_field(csa);
         if (str2int(csa->field, &i) != 0)
            error(csa, "first vertex number missing or invalid");
         if (!(1 <= i && i <= nv))
            error(csa, "first vertex number %d out of range", i);
         read_field(csa);
         if (str2int(csa->field, &j) != 0)
            error(csa, "second vertex number missing or invalid");
         if (!(1 <= j && j <= nv))
            error(csa, "second vertex number %d out of range", j);
         glp_add_arc(G, i, j);
         end_of_line(csa);
      }
      xprintf("%d lines were read\n", csa->count);
done: if (ret) glp_erase_graph(G, G->v_size, G->a_size);
      if (csa->fp != NULL) xfclose(csa->fp);
      if (flag != NULL) xfree(flag);
      return ret;
}

/***********************************************************************
*  NAME
*
*  glp_write_ccdata - write graph in DIMACS clique/coloring format
*
*  SYNOPSIS
*
*  int glp_write_ccdata(glp_graph *G, int v_wgt, const char *fname);
*
*  DESCRIPTION
*
*  The routine glp_write_ccdata writes the specified graph in DIMACS
*  clique/coloring format to a text file.
*
*  RETURNS
*
*  If the operation was successful, the routine returns zero. Otherwise
*  it prints an error message and returns non-zero. */

int glp_write_ccdata(glp_graph *G, int v_wgt, const char *fname)
{     XFILE *fp;
      glp_vertex *v;
      glp_arc *e;
      int i, count = 0, ret;
      double w;
      if (v_wgt >= 0 && v_wgt > G->v_size - (int)sizeof(double))
         xerror("glp_write_ccdata: v_wgt = %d; invalid offset\n",
            v_wgt);
      xprintf("Writing graph to `%s'\n", fname);
      fp = xfopen(fname, "w");
      if (fp == NULL)
      {  xprintf("Unable to create `%s' - %s\n", fname, xerrmsg());
         ret = 1;
         goto done;
      }
      xfprintf(fp, "c %s\n",
         G->name == NULL ? "unknown" : G->name), count++;
      xfprintf(fp, "p edge %d %d\n", G->nv, G->na), count++;
      if (v_wgt >= 0)
      {  for (i = 1; i <= G->nv; i++)
         {  v = G->v[i];
            memcpy(&w, (char *)v->data + v_wgt, sizeof(double));
            if (w != 1.0)
               xfprintf(fp, "n %d %.*g\n", i, DBL_DIG, w), count++;
         }
      }
      for (i = 1; i <= G->nv; i++)
      {  v = G->v[i];
         for (e = v->out; e != NULL; e = e->t_next)
            xfprintf(fp, "e %d %d\n", e->tail->i, e->head->i), count++;
      }
      xfprintf(fp, "c eof\n"), count++;
      xfflush(fp);
      if (xferror(fp))
      {  xprintf("Write error on `%s' - %s\n", fname, xerrmsg());
         ret = 1;
         goto done;
      }
      xprintf("%d lines were written\n", count);
      ret = 0;
done: if (fp != NULL) xfclose(fp);
      return ret;
}

/***********************************************************************
*  NAME
*
*  glp_read_prob - read problem data in GLPK format
*
*  SYNOPSIS
*
*  int glp_read_prob(glp_prob *P, int flags, const char *fname);
*
*  The routine glp_read_prob reads problem data in GLPK LP/MIP format
*  from a text file.
*
*  RETURNS
*
*  If the operation was successful, the routine returns zero. Otherwise
*  it prints an error message and returns non-zero. */

int glp_read_prob(glp_prob *P, int flags, const char *fname)
{     struct csa _csa, *csa = &_csa;
      int mip, m, n, nnz, ne, i, j, k, type, kind, ret, *ln = NULL,
         *ia = NULL, *ja = NULL;
      double lb, ub, temp, *ar = NULL;
      char *rf = NULL, *cf = NULL;
      if (P == NULL || P->magic != GLP_PROB_MAGIC)
         xerror("glp_read_prob: P = %p; invalid problem object\n",
            P);
      if (flags != 0)
         xerror("glp_read_prob: flags = %d; invalid parameter\n",
            flags);
      if (fname == NULL)
         xerror("glp_read_prob: fname = %d; invalid parameter\n",
            fname);
      glp_erase_prob(P);
      if (setjmp(csa->jump))
      {  ret = 1;
         goto done;
      }
      csa->fname = fname;
      csa->fp = NULL;
      csa->count = 0;
      csa->c = '\n';
      csa->field[0] = '\0';
      csa->empty = csa->nonint = 0;
      xprintf("Reading problem data from `%s'...\n", fname);
      csa->fp = xfopen(fname, "r");
      if (csa->fp == NULL)
      {  xprintf("Unable to open `%s' - %s\n", fname, xerrmsg());
         longjmp(csa->jump, 1);
      }
      /* read problem line */
      read_designator(csa);
      if (strcmp(csa->field, "p") != 0)
         error(csa, "problem line missing or invalid");
      read_field(csa);
      if (strcmp(csa->field, "lp") == 0)
         mip = 0;
      else if (strcmp(csa->field, "mip") == 0)
         mip = 1;
      else
         error(csa, "wrong problem designator; `lp' or `mip' expected\n"
            );
      read_field(csa);
      if (strcmp(csa->field, "min") == 0)
         glp_set_obj_dir(P, GLP_MIN);
      else if (strcmp(csa->field, "max") == 0)
         glp_set_obj_dir(P, GLP_MAX);
      else
         error(csa, "objective sense missing or invalid");
      read_field(csa);
      if (!(str2int(csa->field, &m) == 0 && m >= 0))
         error(csa, "number of rows missing or invalid");
      read_field(csa);
      if (!(str2int(csa->field, &n) == 0 && n >= 0))
         error(csa, "number of columns missing or invalid");
      read_field(csa);
      if (!(str2int(csa->field, &nnz) == 0 && nnz >= 0))
         error(csa, "number of constraint coefficients missing or inval"
            "id");
      if (m > 0)
      {  glp_add_rows(P, m);
         for (i = 1; i <= m; i++)
            glp_set_row_bnds(P, i, GLP_FX, 0.0, 0.0);
      }
      if (n > 0)
      {  glp_add_cols(P, n);
         for (j = 1; j <= n; j++)
         {  if (!mip)
               glp_set_col_bnds(P, j, GLP_LO, 0.0, 0.0);
            else
               glp_set_col_kind(P, j, GLP_BV);
         }
      }
      end_of_line(csa);
      /* allocate working arrays */
      rf = xcalloc(1+m, sizeof(char));
      memset(rf, 0, 1+m);
      cf = xcalloc(1+n, sizeof(char));
      memset(cf, 0, 1+n);
      ln = xcalloc(1+nnz, sizeof(int));
      ia = xcalloc(1+nnz, sizeof(int));
      ja = xcalloc(1+nnz, sizeof(int));
      ar = xcalloc(1+nnz, sizeof(double));
      /* read descriptor lines */
      ne = 0;
      for (;;)
      {  read_designator(csa);
         if (strcmp(csa->field, "i") == 0)
         {  /* row descriptor */
            read_field(csa);
            if (str2int(csa->field, &i) != 0)
               error(csa, "row number missing or invalid");
            if (!(1 <= i && i <= m))
               error(csa, "row number out of range");
            read_field(csa);
            if (strcmp(csa->field, "f") == 0)
               type = GLP_FR;
            else if (strcmp(csa->field, "l") == 0)
               type = GLP_LO;
            else if (strcmp(csa->field, "u") == 0)
               type = GLP_UP;
            else if (strcmp(csa->field, "d") == 0)
               type = GLP_DB;
            else if (strcmp(csa->field, "s") == 0)
               type = GLP_FX;
            else
               error(csa, "row type missing or invalid");
            if (type == GLP_LO || type == GLP_DB || type == GLP_FX)
            {  read_field(csa);
               if (str2num(csa->field, &lb) != 0)
                  error(csa, "row lower bound/fixed value missing or in"
                     "valid");
            }
            else
               lb = 0.0;
            if (type == GLP_UP || type == GLP_DB)
            {  read_field(csa);
               if (str2num(csa->field, &ub) != 0)
                  error(csa, "row upper bound missing or invalid");
            }
            else
               ub = 0.0;
            if (rf[i] & 0x01)
               error(csa, "duplicate row descriptor");
            glp_set_row_bnds(P, i, type, lb, ub), rf[i] |= 0x01;
         }
         else if (strcmp(csa->field, "j") == 0)
         {  /* column descriptor */
            read_field(csa);
            if (str2int(csa->field, &j) != 0)
               error(csa, "column number missing or invalid");
            if (!(1 <= j && j <= n))
               error(csa, "column number out of range");
            if (!mip)
               kind = GLP_CV;
            else
            {  read_field(csa);
               if (strcmp(csa->field, "c") == 0)
                  kind = GLP_CV;
               else if (strcmp(csa->field, "i") == 0)
                  kind = GLP_IV;
               else if (strcmp(csa->field, "b") == 0)
               {  kind = GLP_IV;
                  type = GLP_DB, lb = 0.0, ub = 1.0;
                  goto skip;
               }
               else
                  error(csa, "column kind missing or invalid");
            }
            read_field(csa);
            if (strcmp(csa->field, "f") == 0)
               type = GLP_FR;
            else if (strcmp(csa->field, "l") == 0)
               type = GLP_LO;
            else if (strcmp(csa->field, "u") == 0)
               type = GLP_UP;
            else if (strcmp(csa->field, "d") == 0)
               type = GLP_DB;
            else if (strcmp(csa->field, "s") == 0)
               type = GLP_FX;
            else
               error(csa, "column type missing or invalid");
            if (type == GLP_LO || type == GLP_DB || type == GLP_FX)
            {  read_field(csa);
               if (str2num(csa->field, &lb) != 0)
                  error(csa, "column lower bound/fixed value missing or"
                     " invalid");
            }
            else
               lb = 0.0;
            if (type == GLP_UP || type == GLP_DB)
            {  read_field(csa);
               if (str2num(csa->field, &ub) != 0)
                  error(csa, "column upper bound missing or invalid");
            }
            else
               ub = 0.0;
skip:       if (cf[j] & 0x01)
               error(csa, "duplicate column descriptor");
            glp_set_col_kind(P, j, kind);
            glp_set_col_bnds(P, j, type, lb, ub), cf[j] |= 0x01;
         }
         else if (strcmp(csa->field, "a") == 0)
         {  /* coefficient descriptor */
            read_field(csa);
            if (str2int(csa->field, &i) != 0)
               error(csa, "row number missing or invalid");
            if (!(0 <= i && i <= m))
               error(csa, "row number out of range");
            read_field(csa);
            if (str2int(csa->field, &j) != 0)
               error(csa, "column number missing or invalid");
            if (!((i == 0 ? 0 : 1) <= j && j <= n))
               error(csa, "column number out of range");
            read_field(csa);
            if (i == 0)
            {  if (str2num(csa->field, &temp) != 0)
                  error(csa, "objective %s missing or invalid",
                     j == 0 ? "constant term" : "coefficient");
               if (cf[j] & 0x10)
                  error(csa, "duplicate objective %s",
                     j == 0 ? "constant term" : "coefficient");
               glp_set_obj_coef(P, j, temp), cf[j] |= 0x10;
            }
            else
            {  if (str2num(csa->field, &temp) != 0)
                  error(csa, "constraint coefficient missing or invalid"
                     );
               if (ne == nnz)
                  error(csa, "too many constraint coefficient descripto"
                     "rs");
               ln[++ne] = csa->count;
               ia[ne] = i, ja[ne] = j, ar[ne] = temp;
            }
         }
         else if (strcmp(csa->field, "n") == 0)
         {  /* symbolic name descriptor */
            read_field(csa);
            if (strcmp(csa->field, "p") == 0)
            {  /* problem name */
               read_field(csa);
               if (P->name != NULL)
                  error(csa, "duplicate problem name");
               glp_set_prob_name(P, csa->field);
            }
            else if (strcmp(csa->field, "z") == 0)
            {  /* objective name */
               read_field(csa);
               if (P->obj != NULL)
                  error(csa, "duplicate objective name");
               glp_set_obj_name(P, csa->field);
            }
            else if (strcmp(csa->field, "i") == 0)
            {  /* row name */
               read_field(csa);
               if (str2int(csa->field, &i) != 0)
                  error(csa, "row number missing or invalid");
               if (!(1 <= i && i <= m))
                  error(csa, "row number out of range");
               read_field(csa);
               if (P->row[i]->name != NULL)
                  error(csa, "duplicate row name");
               glp_set_row_name(P, i, csa->field);
            }
            else if (strcmp(csa->field, "j") == 0)
            {  /* column name */
               read_field(csa);
               if (str2int(csa->field, &j) != 0)
                  error(csa, "column number missing or invalid");
               if (!(1 <= j && j <= n))
                  error(csa, "column number out of range");
               read_field(csa);
               if (P->col[j]->name != NULL)
                  error(csa, "duplicate column name");
               glp_set_col_name(P, j, csa->field);
            }
            else
               error(csa, "object designator missing or invalid");
         }
         else if (strcmp(csa->field, "e") == 0)
            break;
         else
            error(csa, "line designator missing or invalid");
         end_of_line(csa);
      }
      if (ne < nnz)
         error(csa, "too few constraint coefficient descriptors");
      xassert(ne == nnz);
      k = glp_check_dup(m, n, ne, ia, ja);
      xassert(0 <= k && k <= nnz);
      if (k > 0)
      {  csa->count = ln[k];
         error(csa, "duplicate constraint coefficient");
      }
      glp_load_matrix(P, ne, ia, ja, ar);
      /* print some statistics */
      if (P->name != NULL)
         xprintf("Problem: %s\n", P->name);
      if (P->obj != NULL)
         xprintf("Objective: %s\n", P->obj);
      xprintf("%d row%s, %d column%s, %d non-zero%s\n",
         m, m == 1 ? "" : "s", n, n == 1 ? "" : "s", nnz, nnz == 1 ?
         "" : "s");
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
      glp_sort_matrix(P);
      ret = 0;
done: if (csa->fp != NULL) xfclose(csa->fp);
      if (rf != NULL) xfree(rf);
      if (cf != NULL) xfree(cf);
      if (ln != NULL) xfree(ln);
      if (ia != NULL) xfree(ia);
      if (ja != NULL) xfree(ja);
      if (ar != NULL) xfree(ar);
      if (ret) glp_erase_prob(P);
      return ret;
}

/***********************************************************************
*  NAME
*
*  glp_write_prob - write problem data in GLPK format
*
*  SYNOPSIS
*
*  int glp_write_prob(glp_prob *P, int flags, const char *fname);
*
*  The routine glp_write_prob writes problem data in GLPK LP/MIP format
*  to a text file.
*
*  RETURNS
*
*  If the operation was successful, the routine returns zero. Otherwise
*  it prints an error message and returns non-zero. */

int glp_write_prob(glp_prob *P, int flags, const char *fname)
{     XFILE *fp;
      GLPROW *row;
      GLPCOL *col;
      GLPAIJ *aij;
      int mip, i, j, count, ret;
      if (P == NULL || P->magic != GLP_PROB_MAGIC)
         xerror("glp_write_prob: P = %p; invalid problem object\n",
            P);
      if (flags != 0)
         xerror("glp_write_prob: flags = %d; invalid parameter\n",
            flags);
      if (fname == NULL)
         xerror("glp_write_prob: fname = %d; invalid parameter\n",
            fname);
      xprintf("Writing problem data to `%s'...\n", fname);
      fp = xfopen(fname, "w"), count = 0;
      if (fp == NULL)
      {  xprintf("Unable to create `%s' - %s\n", fname, xerrmsg());
         ret = 1;
         goto done;
      }
      /* write problem line */
      mip = (glp_get_num_int(P) > 0);
      xfprintf(fp, "p %s %s %d %d %d\n", !mip ? "lp" : "mip",
         P->dir == GLP_MIN ? "min" : P->dir == GLP_MAX ? "max" : "???",
         P->m, P->n, P->nnz), count++;
      if (P->name != NULL)
         xfprintf(fp, "n p %s\n", P->name), count++;
      if (P->obj != NULL)
         xfprintf(fp, "n z %s\n", P->obj), count++;
      /* write row descriptors */
      for (i = 1; i <= P->m; i++)
      {  row = P->row[i];
         if (row->type == GLP_FX && row->lb == 0.0)
            goto skip1;
         xfprintf(fp, "i %d ", i), count++;
         if (row->type == GLP_FR)
            xfprintf(fp, "f\n");
         else if (row->type == GLP_LO)
            xfprintf(fp, "l %.*g\n", DBL_DIG, row->lb);
         else if (row->type == GLP_UP)
            xfprintf(fp, "u %.*g\n", DBL_DIG, row->ub);
         else if (row->type == GLP_DB)
            xfprintf(fp, "d %.*g %.*g\n", DBL_DIG, row->lb, DBL_DIG,
                  row->ub);
         else if (row->type == GLP_FX)
            xfprintf(fp, "s %.*g\n", DBL_DIG, row->lb);
         else
            xassert(row != row);
skip1:   if (row->name != NULL)
            xfprintf(fp, "n i %d %s\n", i, row->name), count++;
      }
      /* write column descriptors */
      for (j = 1; j <= P->n; j++)
      {  col = P->col[j];
         if (!mip && col->type == GLP_LO && col->lb == 0.0)
            goto skip2;
         if (mip && col->kind == GLP_IV && col->type == GLP_DB &&
             col->lb == 0.0 && col->ub == 1.0)
            goto skip2;
         xfprintf(fp, "j %d ", j), count++;
         if (mip)
         {  if (col->kind == GLP_CV)
               xfprintf(fp, "c ");
            else if (col->kind == GLP_IV)
               xfprintf(fp, "i ");
            else
               xassert(col != col);
         }
         if (col->type == GLP_FR)
            xfprintf(fp, "f\n");
         else if (col->type == GLP_LO)
            xfprintf(fp, "l %.*g\n", DBL_DIG, col->lb);
         else if (col->type == GLP_UP)
            xfprintf(fp, "u %.*g\n", DBL_DIG, col->ub);
         else if (col->type == GLP_DB)
            xfprintf(fp, "d %.*g %.*g\n", DBL_DIG, col->lb, DBL_DIG,
                  col->ub);
         else if (col->type == GLP_FX)
            xfprintf(fp, "s %.*g\n", DBL_DIG, col->lb);
         else
            xassert(col != col);
skip2:   if (col->name != NULL)
            xfprintf(fp, "n j %d %s\n", j, col->name), count++;
      }
      /* write objective coefficient descriptors */
      if (P->c0 != 0.0)
         xfprintf(fp, "a 0 0 %.*g\n", DBL_DIG, P->c0), count++;
      for (j = 1; j <= P->n; j++)
      {  col = P->col[j];
         if (col->coef != 0.0)
            xfprintf(fp, "a 0 %d %.*g\n", j, DBL_DIG, col->coef),
               count++;
      }
      /* write constraint coefficient descriptors */
      for (i = 1; i <= P->m; i++)
      {  row = P->row[i];
         for (aij = row->ptr; aij != NULL; aij = aij->r_next)
            xfprintf(fp, "a %d %d %.*g\n", i, aij->col->j, DBL_DIG,
               aij->val), count++;
      }
      /* write end line */
      xfprintf(fp, "e o f\n"), count++;
      xfflush(fp);
      if (xferror(fp))
      {  xprintf("Write error on `%s' - %s\n", fname, xerrmsg());
         ret = 1;
         goto done;
      }
      xprintf("%d lines were written\n", count);
      ret = 0;
done: if (fp != NULL) xfclose(fp);
      return ret;
}

/* eof */
