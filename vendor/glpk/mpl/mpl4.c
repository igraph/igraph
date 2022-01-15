/* mpl4.c */

/***********************************************************************
*  This code is part of GLPK (GNU Linear Programming Kit).
*  Copyright (C) 2003-2016 Free Software Foundation, Inc.
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

#include "mpl.h"

#define xfault xerror
#define xfprintf glp_format
#define dmp_create_poolx(size) dmp_create_pool()

/**********************************************************************/
/* * *              GENERATING AND POSTSOLVING MODEL              * * */
/**********************************************************************/

/*----------------------------------------------------------------------
-- alloc_content - allocate content arrays for all model objects.
--
-- This routine allocates content arrays for all existing model objects
-- and thereby finalizes creating model.
--
-- This routine must be called immediately after reading model section,
-- i.e. before reading data section or generating model. */

void alloc_content(MPL *mpl)
{     STATEMENT *stmt;
      /* walk through all model statements */
      for (stmt = mpl->model; stmt != NULL; stmt = stmt->next)
      {  switch (stmt->type)
         {  case A_SET:
               /* model set */
               xassert(stmt->u.set->array == NULL);
               stmt->u.set->array = create_array(mpl, A_ELEMSET,
                  stmt->u.set->dim);
               break;
            case A_PARAMETER:
               /* model parameter */
               xassert(stmt->u.par->array == NULL);
               switch (stmt->u.par->type)
               {  case A_NUMERIC:
                  case A_INTEGER:
                  case A_BINARY:
                     stmt->u.par->array = create_array(mpl, A_NUMERIC,
                        stmt->u.par->dim);
                     break;
                  case A_SYMBOLIC:
                     stmt->u.par->array = create_array(mpl, A_SYMBOLIC,
                        stmt->u.par->dim);
                     break;
                  default:
                     xassert(stmt != stmt);
               }
               break;
            case A_VARIABLE:
               /* model variable */
               xassert(stmt->u.var->array == NULL);
               stmt->u.var->array = create_array(mpl, A_ELEMVAR,
                  stmt->u.var->dim);
               break;
            case A_CONSTRAINT:
               /* model constraint/objective */
               xassert(stmt->u.con->array == NULL);
               stmt->u.con->array = create_array(mpl, A_ELEMCON,
                  stmt->u.con->dim);
               break;
#if 1 /* 11/II-2008 */
            case A_TABLE:
#endif
            case A_SOLVE:
            case A_CHECK:
            case A_DISPLAY:
            case A_PRINTF:
            case A_FOR:
               /* functional statements have no content array */
               break;
            default:
               xassert(stmt != stmt);
         }
      }
      return;
}

/*----------------------------------------------------------------------
-- generate_model - generate model.
--
-- This routine executes the model statements which precede the solve
-- statement. */

void generate_model(MPL *mpl)
{     STATEMENT *stmt;
      xassert(!mpl->flag_p);
      for (stmt = mpl->model; stmt != NULL; stmt = stmt->next)
      {  execute_statement(mpl, stmt);
         if (mpl->stmt->type == A_SOLVE) break;
      }
      mpl->stmt = stmt;
      return;
}

/*----------------------------------------------------------------------
-- build_problem - build problem instance.
--
-- This routine builds lists of rows and columns for problem instance,
-- which corresponds to the generated model. */

void build_problem(MPL *mpl)
{     STATEMENT *stmt;
      MEMBER *memb;
      VARIABLE *v;
      CONSTRAINT *c;
      FORMULA *t;
      int i, j;
      xassert(mpl->m == 0);
      xassert(mpl->n == 0);
      xassert(mpl->row == NULL);
      xassert(mpl->col == NULL);
      /* check that all elemental variables has zero column numbers */
      for (stmt = mpl->model; stmt != NULL; stmt = stmt->next)
      {  if (stmt->type == A_VARIABLE)
         {  v = stmt->u.var;
            for (memb = v->array->head; memb != NULL; memb = memb->next)
               xassert(memb->value.var->j == 0);
         }
      }
      /* assign row numbers to elemental constraints and objectives */
      for (stmt = mpl->model; stmt != NULL; stmt = stmt->next)
      {  if (stmt->type == A_CONSTRAINT)
         {  c = stmt->u.con;
            for (memb = c->array->head; memb != NULL; memb = memb->next)
            {  xassert(memb->value.con->i == 0);
               memb->value.con->i = ++mpl->m;
               /* walk through linear form and mark elemental variables,
                  which are referenced at least once */
               for (t = memb->value.con->form; t != NULL; t = t->next)
               {  xassert(t->var != NULL);
                  t->var->memb->value.var->j = -1;
               }
            }
         }
      }
      /* assign column numbers to marked elemental variables */
      for (stmt = mpl->model; stmt != NULL; stmt = stmt->next)
      {  if (stmt->type == A_VARIABLE)
         {  v = stmt->u.var;
            for (memb = v->array->head; memb != NULL; memb = memb->next)
               if (memb->value.var->j != 0) memb->value.var->j =
                  ++mpl->n;
         }
      }
      /* build list of rows */
      mpl->row = xcalloc(1+mpl->m, sizeof(ELEMCON *));
      for (i = 1; i <= mpl->m; i++) mpl->row[i] = NULL;
      for (stmt = mpl->model; stmt != NULL; stmt = stmt->next)
      {  if (stmt->type == A_CONSTRAINT)
         {  c = stmt->u.con;
            for (memb = c->array->head; memb != NULL; memb = memb->next)
            {  i = memb->value.con->i;
               xassert(1 <= i && i <= mpl->m);
               xassert(mpl->row[i] == NULL);
               mpl->row[i] = memb->value.con;
            }
         }
      }
      for (i = 1; i <= mpl->m; i++) xassert(mpl->row[i] != NULL);
      /* build list of columns */
      mpl->col = xcalloc(1+mpl->n, sizeof(ELEMVAR *));
      for (j = 1; j <= mpl->n; j++) mpl->col[j] = NULL;
      for (stmt = mpl->model; stmt != NULL; stmt = stmt->next)
      {  if (stmt->type == A_VARIABLE)
         {  v = stmt->u.var;
            for (memb = v->array->head; memb != NULL; memb = memb->next)
            {  j = memb->value.var->j;
               if (j == 0) continue;
               xassert(1 <= j && j <= mpl->n);
               xassert(mpl->col[j] == NULL);
               mpl->col[j] = memb->value.var;
            }
         }
      }
      for (j = 1; j <= mpl->n; j++) xassert(mpl->col[j] != NULL);
      return;
}

/*----------------------------------------------------------------------
-- postsolve_model - postsolve model.
--
-- This routine executes the model statements which follow the solve
-- statement. */

void postsolve_model(MPL *mpl)
{     STATEMENT *stmt;
      xassert(!mpl->flag_p);
      mpl->flag_p = 1;
      for (stmt = mpl->stmt; stmt != NULL; stmt = stmt->next)
         execute_statement(mpl, stmt);
      mpl->stmt = NULL;
      return;
}

/*----------------------------------------------------------------------
-- clean_model - clean model content.
--
-- This routine cleans the model content that assumes deleting all stuff
-- dynamically allocated on generating/postsolving phase.
--
-- Actually cleaning model content is not needed. This function is used
-- mainly to be sure that there were no logical errors on using dynamic
-- memory pools during the generation phase.
--
-- NOTE: This routine must not be called if any errors were detected on
--       the generation phase. */

void clean_model(MPL *mpl)
{     STATEMENT *stmt;
      for (stmt = mpl->model; stmt != NULL; stmt = stmt->next)
         clean_statement(mpl, stmt);
      /* check that all atoms have been returned to their pools */
      if (dmp_in_use(mpl->strings) != 0)
         error(mpl, "internal logic error: %d string segment(s) were lo"
            "st", dmp_in_use(mpl->strings));
      if (dmp_in_use(mpl->symbols) != 0)
         error(mpl, "internal logic error: %d symbol(s) were lost",
            dmp_in_use(mpl->symbols));
      if (dmp_in_use(mpl->tuples) != 0)
         error(mpl, "internal logic error: %d n-tuple component(s) were"
            " lost", dmp_in_use(mpl->tuples));
      if (dmp_in_use(mpl->arrays) != 0)
         error(mpl, "internal logic error: %d array(s) were lost",
            dmp_in_use(mpl->arrays));
      if (dmp_in_use(mpl->members) != 0)
         error(mpl, "internal logic error: %d array member(s) were lost"
            , dmp_in_use(mpl->members));
      if (dmp_in_use(mpl->elemvars) != 0)
         error(mpl, "internal logic error: %d elemental variable(s) wer"
            "e lost", dmp_in_use(mpl->elemvars));
      if (dmp_in_use(mpl->formulae) != 0)
         error(mpl, "internal logic error: %d linear term(s) were lost",
            dmp_in_use(mpl->formulae));
      if (dmp_in_use(mpl->elemcons) != 0)
         error(mpl, "internal logic error: %d elemental constraint(s) w"
            "ere lost", dmp_in_use(mpl->elemcons));
      return;
}

/**********************************************************************/
/* * *                        INPUT/OUTPUT                        * * */
/**********************************************************************/

/*----------------------------------------------------------------------
-- open_input - open input text file.
--
-- This routine opens the input text file for scanning. */

void open_input(MPL *mpl, char *file)
{     mpl->line = 0;
      mpl->c = '\n';
      mpl->token = 0;
      mpl->imlen = 0;
      mpl->image[0] = '\0';
      mpl->value = 0.0;
      mpl->b_token = T_EOF;
      mpl->b_imlen = 0;
      mpl->b_image[0] = '\0';
      mpl->b_value = 0.0;
      mpl->f_dots = 0;
      mpl->f_scan = 0;
      mpl->f_token = 0;
      mpl->f_imlen = 0;
      mpl->f_image[0] = '\0';
      mpl->f_value = 0.0;
      memset(mpl->context, ' ', CONTEXT_SIZE);
      mpl->c_ptr = 0;
      xassert(mpl->in_fp == NULL);
      mpl->in_fp = glp_open(file, "r");
      if (mpl->in_fp == NULL)
         error(mpl, "unable to open %s - %s", file, get_err_msg());
      mpl->in_file = file;
      /* scan the very first character */
      get_char(mpl);
      /* scan the very first token */
      get_token(mpl);
      return;
}

/*----------------------------------------------------------------------
-- read_char - read next character from input text file.
--
-- This routine returns a next ASCII character read from the input text
-- file. If the end of file has been reached, EOF is returned. */

int read_char(MPL *mpl)
{     int c;
      xassert(mpl->in_fp != NULL);
      c = glp_getc(mpl->in_fp);
      if (c < 0)
      {  if (glp_ioerr(mpl->in_fp))
            error(mpl, "read error on %s - %s", mpl->in_file,
               get_err_msg());
         c = EOF;
      }
      return c;
}

/*----------------------------------------------------------------------
-- close_input - close input text file.
--
-- This routine closes the input text file. */

void close_input(MPL *mpl)
{     xassert(mpl->in_fp != NULL);
      glp_close(mpl->in_fp);
      mpl->in_fp = NULL;
      mpl->in_file = NULL;
      return;
}

/*----------------------------------------------------------------------
-- open_output - open output text file.
--
-- This routine opens the output text file for writing data produced by
-- display and printf statements. */

void open_output(MPL *mpl, char *file)
{     xassert(mpl->out_fp == NULL);
      if (file == NULL)
      {  file = "<stdout>";
         mpl->out_fp = (void *)stdout;
      }
      else
      {  mpl->out_fp = glp_open(file, "w");
         if (mpl->out_fp == NULL)
            error(mpl, "unable to create %s - %s", file, get_err_msg());
      }
      mpl->out_file = xmalloc(strlen(file)+1);
      strcpy(mpl->out_file, file);
      return;
}

/*----------------------------------------------------------------------
-- write_char - write next character to output text file.
--
-- This routine writes an ASCII character to the output text file. */

void write_char(MPL *mpl, int c)
{     xassert(mpl->out_fp != NULL);
      if (mpl->out_fp == (void *)stdout)
         xprintf("%c", c);
      else
         xfprintf(mpl->out_fp, "%c", c);
      return;
}

/*----------------------------------------------------------------------
-- write_text - format and write text to output text file.
--
-- This routine formats a text using the format control string and then
-- writes this text to the output text file. */

void write_text(MPL *mpl, char *fmt, ...)
{     va_list arg;
      char buf[OUTBUF_SIZE], *c;
      va_start(arg, fmt);
      vsprintf(buf, fmt, arg);
      xassert(strlen(buf) < sizeof(buf));
      va_end(arg);
      for (c = buf; *c != '\0'; c++) write_char(mpl, *c);
      return;
}

/*----------------------------------------------------------------------
-- flush_output - finalize writing data to output text file.
--
-- This routine finalizes writing data to the output text file. */

void flush_output(MPL *mpl)
{     xassert(mpl->out_fp != NULL);
      if (mpl->out_fp != (void *)stdout)
      {
#if 0 /* FIXME */
         xfflush(mpl->out_fp);
#endif
         if (glp_ioerr(mpl->out_fp))
            error(mpl, "write error on %s - %s", mpl->out_file,
               get_err_msg());
      }
      return;
}

/**********************************************************************/
/* * *                      SOLVER INTERFACE                      * * */
/**********************************************************************/

/*----------------------------------------------------------------------
-- error - print error message and terminate model processing.
--
-- This routine formats and prints an error message and then terminates
-- model processing. */

void error(MPL *mpl, char *fmt, ...)
{     va_list arg;
      char msg[4095+1];
      va_start(arg, fmt);
      vsprintf(msg, fmt, arg);
      xassert(strlen(msg) < sizeof(msg));
      va_end(arg);
      switch (mpl->phase)
      {  case 1:
         case 2:
            /* translation phase */
            xprintf("%s:%d: %s\n",
               mpl->in_file == NULL ? "(unknown)" : mpl->in_file,
               mpl->line, msg);
            print_context(mpl);
            break;
         case 3:
            /* generation/postsolve phase */
            xprintf("%s:%d: %s\n",
               mpl->mod_file == NULL ? "(unknown)" : mpl->mod_file,
               mpl->stmt == NULL ? 0 : mpl->stmt->line, msg);
            break;
         default:
            xassert(mpl != mpl);
      }
      mpl->phase = 4;
      longjmp(mpl->jump, 1);
      /* no return */
}

/*----------------------------------------------------------------------
-- warning - print warning message and continue model processing.
--
-- This routine formats and prints a warning message and returns to the
-- calling program. */

void warning(MPL *mpl, char *fmt, ...)
{     va_list arg;
      char msg[4095+1];
      va_start(arg, fmt);
      vsprintf(msg, fmt, arg);
      xassert(strlen(msg) < sizeof(msg));
      va_end(arg);
      switch (mpl->phase)
      {  case 1:
         case 2:
            /* translation phase */
            xprintf("%s:%d: warning: %s\n",
               mpl->in_file == NULL ? "(unknown)" : mpl->in_file,
               mpl->line, msg);
            break;
         case 3:
            /* generation/postsolve phase */
            xprintf("%s:%d: warning: %s\n",
               mpl->mod_file == NULL ? "(unknown)" : mpl->mod_file,
               mpl->stmt == NULL ? 0 : mpl->stmt->line, msg);
            break;
         default:
            xassert(mpl != mpl);
      }
      return;
}

/*----------------------------------------------------------------------
-- mpl_initialize - create and initialize translator database.
--
-- *Synopsis*
--
-- #include "glpmpl.h"
-- MPL *mpl_initialize(void);
--
-- *Description*
--
-- The routine mpl_initialize creates and initializes the database used
-- by the GNU MathProg translator.
--
-- *Returns*
--
-- The routine returns a pointer to the database created. */

MPL *mpl_initialize(void)
{     MPL *mpl;
      mpl = xmalloc(sizeof(MPL));
      /* scanning segment */
      mpl->line = 0;
      mpl->c = 0;
      mpl->token = 0;
      mpl->imlen = 0;
      mpl->image = xcalloc(MAX_LENGTH+1, sizeof(char));
      mpl->image[0] = '\0';
      mpl->value = 0.0;
      mpl->b_token = 0;
      mpl->b_imlen = 0;
      mpl->b_image = xcalloc(MAX_LENGTH+1, sizeof(char));
      mpl->b_image[0] = '\0';
      mpl->b_value = 0.0;
      mpl->f_dots = 0;
      mpl->f_scan = 0;
      mpl->f_token = 0;
      mpl->f_imlen = 0;
      mpl->f_image = xcalloc(MAX_LENGTH+1, sizeof(char));
      mpl->f_image[0] = '\0';
      mpl->f_value = 0.0;
      mpl->context = xcalloc(CONTEXT_SIZE, sizeof(char));
      memset(mpl->context, ' ', CONTEXT_SIZE);
      mpl->c_ptr = 0;
      mpl->flag_d = 0;
      /* translating segment */
      mpl->pool = dmp_create_poolx(0);
      mpl->tree = avl_create_tree(avl_strcmp, NULL);
      mpl->model = NULL;
      mpl->flag_x = 0;
      mpl->as_within = 0;
      mpl->as_in = 0;
      mpl->as_binary = 0;
      mpl->flag_s = 0;
      /* common segment */
      mpl->strings = dmp_create_poolx(sizeof(STRING));
      mpl->symbols = dmp_create_poolx(sizeof(SYMBOL));
      mpl->tuples = dmp_create_poolx(sizeof(TUPLE));
      mpl->arrays = dmp_create_poolx(sizeof(ARRAY));
      mpl->members = dmp_create_poolx(sizeof(MEMBER));
      mpl->elemvars = dmp_create_poolx(sizeof(ELEMVAR));
      mpl->formulae = dmp_create_poolx(sizeof(FORMULA));
      mpl->elemcons = dmp_create_poolx(sizeof(ELEMCON));
      mpl->a_list = NULL;
      mpl->sym_buf = xcalloc(255+1, sizeof(char));
      mpl->sym_buf[0] = '\0';
      mpl->tup_buf = xcalloc(255+1, sizeof(char));
      mpl->tup_buf[0] = '\0';
      /* generating/postsolving segment */
      mpl->rand = rng_create_rand();
      mpl->flag_p = 0;
      mpl->stmt = NULL;
#if 1 /* 11/II-2008 */
      mpl->dca = NULL;
#endif
      mpl->m = 0;
      mpl->n = 0;
      mpl->row = NULL;
      mpl->col = NULL;
      /* input/output segment */
      mpl->in_fp = NULL;
      mpl->in_file = NULL;
      mpl->out_fp = NULL;
      mpl->out_file = NULL;
      mpl->prt_fp = NULL;
      mpl->prt_file = NULL;
      /* solver interface segment */
      if (setjmp(mpl->jump)) xassert(mpl != mpl);
      mpl->phase = 0;
      mpl->mod_file = NULL;
      mpl->mpl_buf = xcalloc(255+1, sizeof(char));
      mpl->mpl_buf[0] = '\0';
      return mpl;
}

/*----------------------------------------------------------------------
-- mpl_read_model - read model section and optional data section.
--
-- *Synopsis*
--
-- #include "glpmpl.h"
-- int mpl_read_model(MPL *mpl, char *file, int skip_data);
--
-- *Description*
--
-- The routine mpl_read_model reads model section and optionally data
-- section, which may follow the model section, from the text file,
-- whose name is the character string file, performs translating model
-- statements and data blocks, and stores all the information in the
-- translator database.
--
-- The parameter skip_data is a flag. If the input file contains the
-- data section and this flag is set, the data section is not read as
-- if there were no data section and a warning message is issued. This
-- allows reading the data section from another input file.
--
-- This routine should be called once after the routine mpl_initialize
-- and before other API routines.
--
-- *Returns*
--
-- The routine mpl_read_model returns one the following codes:
--
-- 1 - translation successful. The input text file contains only model
--     section. In this case the calling program may call the routine
--     mpl_read_data to read data section from another file.
-- 2 - translation successful. The input text file contains both model
--     and data section.
-- 4 - processing failed due to some errors. In this case the calling
--     program should call the routine mpl_terminate to terminate model
--     processing. */

int mpl_read_model(MPL *mpl, char *file, int skip_data)
{     if (mpl->phase != 0)
         xfault("mpl_read_model: invalid call sequence\n");
      if (file == NULL)
         xfault("mpl_read_model: no input filename specified\n");
      /* set up error handler */
      if (setjmp(mpl->jump)) goto done;
      /* translate model section */
      mpl->phase = 1;
      xprintf("Reading model section from %s...\n", file);
      open_input(mpl, file);
      model_section(mpl);
      if (mpl->model == NULL)
         error(mpl, "empty model section not allowed");
      /* save name of the input text file containing model section for
         error diagnostics during the generation phase */
      mpl->mod_file = xcalloc(strlen(file)+1, sizeof(char));
      strcpy(mpl->mod_file, mpl->in_file);
      /* allocate content arrays for all model objects */
      alloc_content(mpl);
      /* optional data section may begin with the keyword 'data' */
      if (is_keyword(mpl, "data"))
      {  if (skip_data)
         {  warning(mpl, "data section ignored");
            goto skip;
         }
         mpl->flag_d = 1;
         get_token(mpl /* data */);
         if (mpl->token != T_SEMICOLON)
            error(mpl, "semicolon missing where expected");
         get_token(mpl /* ; */);
         /* translate data section */
         mpl->phase = 2;
         xprintf("Reading data section from %s...\n", file);
         data_section(mpl);
      }
      /* process end statement */
      end_statement(mpl);
skip: xprintf("%d line%s were read\n",
         mpl->line, mpl->line == 1 ? "" : "s");
      close_input(mpl);
done: /* return to the calling program */
      return mpl->phase;
}

/*----------------------------------------------------------------------
-- mpl_read_data - read data section.
--
-- *Synopsis*
--
-- #include "glpmpl.h"
-- int mpl_read_data(MPL *mpl, char *file);
--
-- *Description*
--
-- The routine mpl_read_data reads data section from the text file,
-- whose name is the character string file, performs translating data
-- blocks, and stores the data read in the translator database.
--
-- If this routine is used, it should be called once after the routine
-- mpl_read_model and if the latter returned the code 1.
--
-- *Returns*
--
-- The routine mpl_read_data returns one of the following codes:
--
-- 2 - data section has been successfully processed.
-- 4 - processing failed due to some errors. In this case the calling
--     program should call the routine mpl_terminate to terminate model
--     processing. */

int mpl_read_data(MPL *mpl, char *file)
#if 0 /* 02/X-2008 */
{     if (mpl->phase != 1)
#else
{     if (!(mpl->phase == 1 || mpl->phase == 2))
#endif
         xfault("mpl_read_data: invalid call sequence\n");
      if (file == NULL)
         xfault("mpl_read_data: no input filename specified\n");
      /* set up error handler */
      if (setjmp(mpl->jump)) goto done;
      /* process data section */
      mpl->phase = 2;
      xprintf("Reading data section from %s...\n", file);
      mpl->flag_d = 1;
      open_input(mpl, file);
      /* in this case the keyword 'data' is optional */
      if (is_literal(mpl, "data"))
      {  get_token(mpl /* data */);
         if (mpl->token != T_SEMICOLON)
            error(mpl, "semicolon missing where expected");
         get_token(mpl /* ; */);
      }
      data_section(mpl);
      /* process end statement */
      end_statement(mpl);
      xprintf("%d line%s were read\n",
         mpl->line, mpl->line == 1 ? "" : "s");
      close_input(mpl);
done: /* return to the calling program */
      return mpl->phase;
}

/*----------------------------------------------------------------------
-- mpl_generate - generate model.
--
-- *Synopsis*
--
-- #include "glpmpl.h"
-- int mpl_generate(MPL *mpl, char *file);
--
-- *Description*
--
-- The routine mpl_generate generates the model using its description
-- stored in the translator database. This phase means generating all
-- variables, constraints, and objectives, executing check and display
-- statements, which precede the solve statement (if it is presented),
-- and building the problem instance.
--
-- The character string file specifies the name of output text file, to
-- which output produced by display statements should be written. It is
-- allowed to specify NULL, in which case the output goes to stdout via
-- the routine print.
--
-- This routine should be called once after the routine mpl_read_model
-- or mpl_read_data and if one of the latters returned the code 2.
--
-- *Returns*
--
-- The routine mpl_generate returns one of the following codes:
--
-- 3 - model has been successfully generated. In this case the calling
--     program may call other api routines to obtain components of the
--     problem instance from the translator database.
-- 4 - processing failed due to some errors. In this case the calling
--     program should call the routine mpl_terminate to terminate model
--     processing. */

int mpl_generate(MPL *mpl, char *file)
{     if (!(mpl->phase == 1 || mpl->phase == 2))
         xfault("mpl_generate: invalid call sequence\n");
      /* set up error handler */
      if (setjmp(mpl->jump)) goto done;
      /* generate model */
      mpl->phase = 3;
      open_output(mpl, file);
      generate_model(mpl);
      flush_output(mpl);
      /* build problem instance */
      build_problem(mpl);
      /* generation phase has been finished */
      xprintf("Model has been successfully generated\n");
done: /* return to the calling program */
      return mpl->phase;
}

/*----------------------------------------------------------------------
-- mpl_get_prob_name - obtain problem (model) name.
--
-- *Synopsis*
--
-- #include "glpmpl.h"
-- char *mpl_get_prob_name(MPL *mpl);
--
-- *Returns*
--
-- The routine mpl_get_prob_name returns a pointer to internal buffer,
-- which contains symbolic name of the problem (model).
--
-- *Note*
--
-- Currently MathProg has no feature to assign a symbolic name to the
-- model. Therefore the routine mpl_get_prob_name tries to construct
-- such name using the name of input text file containing model section,
-- although this is not a good idea (due to portability problems). */

char *mpl_get_prob_name(MPL *mpl)
{     char *name = mpl->mpl_buf;
      char *file = mpl->mod_file;
      int k;
      if (mpl->phase != 3)
         xfault("mpl_get_prob_name: invalid call sequence\n");
      for (;;)
      {  if (strchr(file, '/') != NULL)
            file = strchr(file, '/') + 1;
         else if (strchr(file, '\\') != NULL)
            file = strchr(file, '\\') + 1;
         else if (strchr(file, ':') != NULL)
            file = strchr(file, ':') + 1;
         else
            break;
      }
      for (k = 0; ; k++)
      {  if (k == 255) break;
         if (!(isalnum((unsigned char)*file) || *file == '_')) break;
         name[k] = *file++;
      }
      if (k == 0)
         strcpy(name, "Unknown");
      else
         name[k] = '\0';
      xassert(strlen(name) <= 255);
      return name;
}

/*----------------------------------------------------------------------
-- mpl_get_num_rows - determine number of rows.
--
-- *Synopsis*
--
-- #include "glpmpl.h"
-- int mpl_get_num_rows(MPL *mpl);
--
-- *Returns*
--
-- The routine mpl_get_num_rows returns total number of rows in the
-- problem, where each row is an individual constraint or objective. */

int mpl_get_num_rows(MPL *mpl)
{     if (mpl->phase != 3)
         xfault("mpl_get_num_rows: invalid call sequence\n");
      return mpl->m;
}

/*----------------------------------------------------------------------
-- mpl_get_num_cols - determine number of columns.
--
-- *Synopsis*
--
-- #include "glpmpl.h"
-- int mpl_get_num_cols(MPL *mpl);
--
-- *Returns*
--
-- The routine mpl_get_num_cols returns total number of columns in the
-- problem, where each column is an individual variable. */

int mpl_get_num_cols(MPL *mpl)
{     if (mpl->phase != 3)
         xfault("mpl_get_num_cols: invalid call sequence\n");
      return mpl->n;
}

/*----------------------------------------------------------------------
-- mpl_get_row_name - obtain row name.
--
-- *Synopsis*
--
-- #include "glpmpl.h"
-- char *mpl_get_row_name(MPL *mpl, int i);
--
-- *Returns*
--
-- The routine mpl_get_row_name returns a pointer to internal buffer,
-- which contains symbolic name of i-th row of the problem. */

char *mpl_get_row_name(MPL *mpl, int i)
{     char *name = mpl->mpl_buf, *t;
      int len;
      if (mpl->phase != 3)
         xfault("mpl_get_row_name: invalid call sequence\n");
      if (!(1 <= i && i <= mpl->m))
         xfault("mpl_get_row_name: i = %d; row number out of range\n",
            i);
      strcpy(name, mpl->row[i]->con->name);
      len = strlen(name);
      xassert(len <= 255);
      t = format_tuple(mpl, '[', mpl->row[i]->memb->tuple);
      while (*t)
      {  if (len == 255) break;
         name[len++] = *t++;
      }
      name[len] = '\0';
      if (len == 255) strcpy(name+252, "...");
      xassert(strlen(name) <= 255);
      return name;
}

/*----------------------------------------------------------------------
-- mpl_get_row_kind - determine row kind.
--
-- *Synopsis*
--
-- #include "glpmpl.h"
-- int mpl_get_row_kind(MPL *mpl, int i);
--
-- *Returns*
--
-- The routine mpl_get_row_kind returns the kind of i-th row, which can
-- be one of the following:
--
-- MPL_ST  - non-free (constraint) row;
-- MPL_MIN - free (objective) row to be minimized;
-- MPL_MAX - free (objective) row to be maximized. */

int mpl_get_row_kind(MPL *mpl, int i)
{     int kind;
      if (mpl->phase != 3)
         xfault("mpl_get_row_kind: invalid call sequence\n");
      if (!(1 <= i && i <= mpl->m))
         xfault("mpl_get_row_kind: i = %d; row number out of range\n",
            i);
      switch (mpl->row[i]->con->type)
      {  case A_CONSTRAINT:
            kind = MPL_ST; break;
         case A_MINIMIZE:
            kind = MPL_MIN; break;
         case A_MAXIMIZE:
            kind = MPL_MAX; break;
         default:
            xassert(mpl != mpl);
      }
      return kind;
}

/*----------------------------------------------------------------------
-- mpl_get_row_bnds - obtain row bounds.
--
-- *Synopsis*
--
-- #include "glpmpl.h"
-- int mpl_get_row_bnds(MPL *mpl, int i, double *lb, double *ub);
--
-- *Description*
--
-- The routine mpl_get_row_bnds stores lower and upper bounds of i-th
-- row of the problem to the locations, which the parameters lb and ub
-- point to, respectively. Besides the routine returns the type of the
-- i-th row.
--
-- If some of the parameters lb and ub is NULL, the corresponding bound
-- value is not stored.
--
-- Types and bounds have the following meaning:
--
--     Type           Bounds          Note
--    -----------------------------------------------------------
--    MPL_FR   -inf <  f(x) <  +inf   Free linear form
--    MPL_LO     lb <= f(x) <  +inf   Inequality f(x) >= lb
--    MPL_UP   -inf <  f(x) <=  ub    Inequality f(x) <= ub
--    MPL_DB     lb <= f(x) <=  ub    Inequality lb <= f(x) <= ub
--    MPL_FX           f(x)  =  lb    Equality f(x) = lb
--
-- where f(x) is the corresponding linear form of the i-th row.
--
-- If the row has no lower bound, *lb is set to zero; if the row has
-- no upper bound, *ub is set to zero; and if the row is of fixed type,
-- both *lb and *ub are set to the same value.
--
-- *Returns*
--
-- The routine returns the type of the i-th row as it is stated in the
-- table above. */

int mpl_get_row_bnds(MPL *mpl, int i, double *_lb, double *_ub)
{     ELEMCON *con;
      int type;
      double lb, ub;
      if (mpl->phase != 3)
         xfault("mpl_get_row_bnds: invalid call sequence\n");
      if (!(1 <= i && i <= mpl->m))
         xfault("mpl_get_row_bnds: i = %d; row number out of range\n",
            i);
      con = mpl->row[i];
#if 0 /* 21/VII-2006 */
      if (con->con->lbnd == NULL && con->con->ubnd == NULL)
         type = MPL_FR, lb = ub = 0.0;
      else if (con->con->ubnd == NULL)
         type = MPL_LO, lb = con->lbnd, ub = 0.0;
      else if (con->con->lbnd == NULL)
         type = MPL_UP, lb = 0.0, ub = con->ubnd;
      else if (con->con->lbnd != con->con->ubnd)
         type = MPL_DB, lb = con->lbnd, ub = con->ubnd;
      else
         type = MPL_FX, lb = ub = con->lbnd;
#else
      lb = (con->con->lbnd == NULL ? -DBL_MAX : con->lbnd);
      ub = (con->con->ubnd == NULL ? +DBL_MAX : con->ubnd);
      if (lb == -DBL_MAX && ub == +DBL_MAX)
         type = MPL_FR, lb = ub = 0.0;
      else if (ub == +DBL_MAX)
         type = MPL_LO, ub = 0.0;
      else if (lb == -DBL_MAX)
         type = MPL_UP, lb = 0.0;
      else if (con->con->lbnd != con->con->ubnd)
         type = MPL_DB;
      else
         type = MPL_FX;
#endif
      if (_lb != NULL) *_lb = lb;
      if (_ub != NULL) *_ub = ub;
      return type;
}

/*----------------------------------------------------------------------
-- mpl_get_mat_row - obtain row of the constraint matrix.
--
-- *Synopsis*
--
-- #include "glpmpl.h"
-- int mpl_get_mat_row(MPL *mpl, int i, int ndx[], double val[]);
--
-- *Description*
--
-- The routine mpl_get_mat_row stores column indices and numeric values
-- of constraint coefficients for the i-th row to locations ndx[1], ...,
-- ndx[len] and val[1], ..., val[len], respectively, where 0 <= len <= n
-- is number of (structural) non-zero constraint coefficients, and n is
-- number of columns in the problem.
--
-- If the parameter ndx is NULL, column indices are not stored. If the
-- parameter val is NULL, numeric values are not stored.
--
-- Note that free rows may have constant terms, which are not part of
-- the constraint matrix and therefore not reported by this routine. The
-- constant term of a particular row can be obtained, if necessary, via
-- the routine mpl_get_row_c0.
--
-- *Returns*
--
-- The routine mpl_get_mat_row returns len, which is length of i-th row
-- of the constraint matrix (i.e. number of non-zero coefficients). */

int mpl_get_mat_row(MPL *mpl, int i, int ndx[], double val[])
{     FORMULA *term;
      int len = 0;
      if (mpl->phase != 3)
         xfault("mpl_get_mat_row: invalid call sequence\n");
      if (!(1 <= i && i <= mpl->m))
         xfault("mpl_get_mat_row: i = %d; row number out of range\n",
            i);
      for (term = mpl->row[i]->form; term != NULL; term = term->next)
      {  xassert(term->var != NULL);
         len++;
         xassert(len <= mpl->n);
         if (ndx != NULL) ndx[len] = term->var->j;
         if (val != NULL) val[len] = term->coef;
      }
      return len;
}

/*----------------------------------------------------------------------
-- mpl_get_row_c0 - obtain constant term of free row.
--
-- *Synopsis*
--
-- #include "glpmpl.h"
-- double mpl_get_row_c0(MPL *mpl, int i);
--
-- *Returns*
--
-- The routine mpl_get_row_c0 returns numeric value of constant term of
-- i-th row.
--
-- Note that only free rows may have non-zero constant terms. Therefore
-- if i-th row is not free, the routine returns zero. */

double mpl_get_row_c0(MPL *mpl, int i)
{     ELEMCON *con;
      double c0;
      if (mpl->phase != 3)
         xfault("mpl_get_row_c0: invalid call sequence\n");
      if (!(1 <= i && i <= mpl->m))
         xfault("mpl_get_row_c0: i = %d; row number out of range\n",
            i);
      con = mpl->row[i];
      if (con->con->lbnd == NULL && con->con->ubnd == NULL)
         c0 = - con->lbnd;
      else
         c0 = 0.0;
      return c0;
}

/*----------------------------------------------------------------------
-- mpl_get_col_name - obtain column name.
--
-- *Synopsis*
--
-- #include "glpmpl.h"
-- char *mpl_get_col_name(MPL *mpl, int j);
--
-- *Returns*
--
-- The routine mpl_get_col_name returns a pointer to internal buffer,
-- which contains symbolic name of j-th column of the problem. */

char *mpl_get_col_name(MPL *mpl, int j)
{     char *name = mpl->mpl_buf, *t;
      int len;
      if (mpl->phase != 3)
         xfault("mpl_get_col_name: invalid call sequence\n");
      if (!(1 <= j && j <= mpl->n))
         xfault("mpl_get_col_name: j = %d; column number out of range\n"
            , j);
      strcpy(name, mpl->col[j]->var->name);
      len = strlen(name);
      xassert(len <= 255);
      t = format_tuple(mpl, '[', mpl->col[j]->memb->tuple);
      while (*t)
      {  if (len == 255) break;
         name[len++] = *t++;
      }
      name[len] = '\0';
      if (len == 255) strcpy(name+252, "...");
      xassert(strlen(name) <= 255);
      return name;
}

/*----------------------------------------------------------------------
-- mpl_get_col_kind - determine column kind.
--
-- *Synopsis*
--
-- #include "glpmpl.h"
-- int mpl_get_col_kind(MPL *mpl, int j);
--
-- *Returns*
--
-- The routine mpl_get_col_kind returns the kind of j-th column, which
-- can be one of the following:
--
-- MPL_NUM - continuous variable;
-- MPL_INT - integer variable;
-- MPL_BIN - binary variable.
--
-- Note that column kinds are defined independently on type and bounds
-- (reported by the routine mpl_get_col_bnds) of corresponding columns.
-- This means, in particular, that bounds of an integer column may be
-- fractional, or a binary column may have lower and upper bounds that
-- are not 0 and 1 (or it may have no lower/upper bound at all). */

int mpl_get_col_kind(MPL *mpl, int j)
{     int kind;
      if (mpl->phase != 3)
         xfault("mpl_get_col_kind: invalid call sequence\n");
      if (!(1 <= j && j <= mpl->n))
         xfault("mpl_get_col_kind: j = %d; column number out of range\n"
            , j);
      switch (mpl->col[j]->var->type)
      {  case A_NUMERIC:
            kind = MPL_NUM; break;
         case A_INTEGER:
            kind = MPL_INT; break;
         case A_BINARY:
            kind = MPL_BIN; break;
         default:
            xassert(mpl != mpl);
      }
      return kind;
}

/*----------------------------------------------------------------------
-- mpl_get_col_bnds - obtain column bounds.
--
-- *Synopsis*
--
-- #include "glpmpl.h"
-- int mpl_get_col_bnds(MPL *mpl, int j, double *lb, double *ub);
--
-- *Description*
--
-- The routine mpl_get_col_bnds stores lower and upper bound of j-th
-- column of the problem to the locations, which the parameters lb and
-- ub point to, respectively. Besides the routine returns the type of
-- the j-th column.
--
-- If some of the parameters lb and ub is NULL, the corresponding bound
-- value is not stored.
--
-- Types and bounds have the following meaning:
--
--     Type         Bounds         Note
--    ------------------------------------------------------
--    MPL_FR   -inf <  x <  +inf   Free (unbounded) variable
--    MPL_LO     lb <= x <  +inf   Variable with lower bound
--    MPL_UP   -inf <  x <=  ub    Variable with upper bound
--    MPL_DB     lb <= x <=  ub    Double-bounded variable
--    MPL_FX           x  =  lb    Fixed variable
--
-- where x is individual variable corresponding to the j-th column.
--
-- If the column has no lower bound, *lb is set to zero; if the column
-- has no upper bound, *ub is set to zero; and if the column is of fixed
-- type, both *lb and *ub are set to the same value.
--
-- *Returns*
--
-- The routine returns the type of the j-th column as it is stated in
-- the table above. */

int mpl_get_col_bnds(MPL *mpl, int j, double *_lb, double *_ub)
{     ELEMVAR *var;
      int type;
      double lb, ub;
      if (mpl->phase != 3)
         xfault("mpl_get_col_bnds: invalid call sequence\n");
      if (!(1 <= j && j <= mpl->n))
         xfault("mpl_get_col_bnds: j = %d; column number out of range\n"
            , j);
      var = mpl->col[j];
#if 0 /* 21/VII-2006 */
      if (var->var->lbnd == NULL && var->var->ubnd == NULL)
         type = MPL_FR, lb = ub = 0.0;
      else if (var->var->ubnd == NULL)
         type = MPL_LO, lb = var->lbnd, ub = 0.0;
      else if (var->var->lbnd == NULL)
         type = MPL_UP, lb = 0.0, ub = var->ubnd;
      else if (var->var->lbnd != var->var->ubnd)
         type = MPL_DB, lb = var->lbnd, ub = var->ubnd;
      else
         type = MPL_FX, lb = ub = var->lbnd;
#else
      lb = (var->var->lbnd == NULL ? -DBL_MAX : var->lbnd);
      ub = (var->var->ubnd == NULL ? +DBL_MAX : var->ubnd);
      if (lb == -DBL_MAX && ub == +DBL_MAX)
         type = MPL_FR, lb = ub = 0.0;
      else if (ub == +DBL_MAX)
         type = MPL_LO, ub = 0.0;
      else if (lb == -DBL_MAX)
         type = MPL_UP, lb = 0.0;
      else if (var->var->lbnd != var->var->ubnd)
         type = MPL_DB;
      else
         type = MPL_FX;
#endif
      if (_lb != NULL) *_lb = lb;
      if (_ub != NULL) *_ub = ub;
      return type;
}

/*----------------------------------------------------------------------
-- mpl_has_solve_stmt - check if model has solve statement.
--
-- *Synopsis*
--
-- #include "glpmpl.h"
-- int mpl_has_solve_stmt(MPL *mpl);
--
-- *Returns*
--
-- If the model has the solve statement, the routine returns non-zero,
-- otherwise zero is returned. */

int mpl_has_solve_stmt(MPL *mpl)
{     if (mpl->phase != 3)
         xfault("mpl_has_solve_stmt: invalid call sequence\n");
      return mpl->flag_s;
}

#if 1 /* 15/V-2010 */
void mpl_put_row_soln(MPL *mpl, int i, int stat, double prim,
      double dual)
{     /* store row (constraint/objective) solution components */
      xassert(mpl->phase == 3);
      xassert(1 <= i && i <= mpl->m);
      mpl->row[i]->stat = stat;
      mpl->row[i]->prim = prim;
      mpl->row[i]->dual = dual;
      return;
}
#endif

#if 1 /* 15/V-2010 */
void mpl_put_col_soln(MPL *mpl, int j, int stat, double prim,
      double dual)
{     /* store column (variable) solution components */
      xassert(mpl->phase == 3);
      xassert(1 <= j && j <= mpl->n);
      mpl->col[j]->stat = stat;
      mpl->col[j]->prim = prim;
      mpl->col[j]->dual = dual;
      return;
}
#endif

#if 0 /* 15/V-2010 */
/*----------------------------------------------------------------------
-- mpl_put_col_value - store column value.
--
-- *Synopsis*
--
-- #include "glpmpl.h"
-- void mpl_put_col_value(MPL *mpl, int j, double val);
--
-- *Description*
--
-- The routine mpl_put_col_value stores numeric value of j-th column
-- into the translator database. It is assumed that the column value is
-- provided by the solver. */

void mpl_put_col_value(MPL *mpl, int j, double val)
{     if (mpl->phase != 3)
         xfault("mpl_put_col_value: invalid call sequence\n");
      if (!(1 <= j && j <= mpl->n))
         xfault(
         "mpl_put_col_value: j = %d; column number out of range\n", j);
      mpl->col[j]->prim = val;
      return;
}
#endif

/*----------------------------------------------------------------------
-- mpl_postsolve - postsolve model.
--
-- *Synopsis*
--
-- #include "glpmpl.h"
-- int mpl_postsolve(MPL *mpl);
--
-- *Description*
--
-- The routine mpl_postsolve performs postsolving of the model using
-- its description stored in the translator database. This phase means
-- executing statements, which follow the solve statement.
--
-- If this routine is used, it should be called once after the routine
-- mpl_generate and if the latter returned the code 3.
--
-- *Returns*
--
-- The routine mpl_postsolve returns one of the following codes:
--
-- 3 - model has been successfully postsolved.
-- 4 - processing failed due to some errors. In this case the calling
--     program should call the routine mpl_terminate to terminate model
--     processing. */

int mpl_postsolve(MPL *mpl)
{     if (!(mpl->phase == 3 && !mpl->flag_p))
         xfault("mpl_postsolve: invalid call sequence\n");
      /* set up error handler */
      if (setjmp(mpl->jump)) goto done;
      /* perform postsolving */
      postsolve_model(mpl);
      flush_output(mpl);
      /* postsolving phase has been finished */
      xprintf("Model has been successfully processed\n");
done: /* return to the calling program */
      return mpl->phase;
}

/*----------------------------------------------------------------------
-- mpl_terminate - free all resources used by translator.
--
-- *Synopsis*
--
-- #include "glpmpl.h"
-- void mpl_terminate(MPL *mpl);
--
-- *Description*
--
-- The routine mpl_terminate frees all the resources used by the GNU
-- MathProg translator. */

void mpl_terminate(MPL *mpl)
{     if (setjmp(mpl->jump)) xassert(mpl != mpl);
      switch (mpl->phase)
      {  case 0:
         case 1:
         case 2:
         case 3:
            /* there were no errors; clean the model content */
            clean_model(mpl);
            xassert(mpl->a_list == NULL);
#if 1 /* 11/II-2008 */
            xassert(mpl->dca == NULL);
#endif
            break;
         case 4:
            /* model processing has been finished due to error; delete
               search trees, which may be created for some arrays */
            {  ARRAY *a;
               for (a = mpl->a_list; a != NULL; a = a->next)
                  if (a->tree != NULL) avl_delete_tree(a->tree);
            }
#if 1 /* 11/II-2008 */
            free_dca(mpl);
#endif
            break;
         default:
            xassert(mpl != mpl);
      }
      /* delete the translator database */
      xfree(mpl->image);
      xfree(mpl->b_image);
      xfree(mpl->f_image);
      xfree(mpl->context);
      dmp_delete_pool(mpl->pool);
      avl_delete_tree(mpl->tree);
      dmp_delete_pool(mpl->strings);
      dmp_delete_pool(mpl->symbols);
      dmp_delete_pool(mpl->tuples);
      dmp_delete_pool(mpl->arrays);
      dmp_delete_pool(mpl->members);
      dmp_delete_pool(mpl->elemvars);
      dmp_delete_pool(mpl->formulae);
      dmp_delete_pool(mpl->elemcons);
      xfree(mpl->sym_buf);
      xfree(mpl->tup_buf);
      rng_delete_rand(mpl->rand);
      if (mpl->row != NULL) xfree(mpl->row);
      if (mpl->col != NULL) xfree(mpl->col);
      if (mpl->in_fp != NULL) glp_close(mpl->in_fp);
      if (mpl->out_fp != NULL && mpl->out_fp != (void *)stdout)
         glp_close(mpl->out_fp);
      if (mpl->out_file != NULL) xfree(mpl->out_file);
      if (mpl->prt_fp != NULL) glp_close(mpl->prt_fp);
      if (mpl->prt_file != NULL) xfree(mpl->prt_file);
      if (mpl->mod_file != NULL) xfree(mpl->mod_file);
      xfree(mpl->mpl_buf);
      xfree(mpl);
      return;
}

/* eof */
