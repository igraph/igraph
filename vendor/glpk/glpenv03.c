/* glpenv03.c (terminal output) */

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
#pragma clang diagnostic ignored "-Wunused-label"
#pragma clang diagnostic ignored "-Wunused-variable"
#endif

#include "glpapi.h"

/***********************************************************************
*  NAME
*
*  glp_printf - write formatted output to terminal
*
*  SYNOPSIS
*
*  void glp_printf(const char *fmt, ...);
*
*  DESCRIPTION
*
*  The routine glp_printf uses the format control string fmt to format
*  its parameters and writes the formatted output to the terminal. */

void glp_printf(const char *fmt, ...)
{     va_list arg;
      va_start(arg, fmt);
      xvprintf(fmt, arg);
      va_end(arg);
      return;
}

/***********************************************************************
*  NAME
*
*  glp_vprintf - write formatted output to terminal
*
*  SYNOPSIS
*
*  void glp_vprintf(const char *fmt, va_list arg);
*
*  DESCRIPTION
*
*  The routine glp_vprintf uses the format control string fmt to format
*  its parameters specified by the list arg and writes the formatted
*  output to the terminal. */

void glp_vprintf(const char *fmt, va_list arg)
{     ENV *env = get_env_ptr();
      /* if terminal output is disabled, do nothing */
      if (!env->term_out) goto skip;
      /* format the output */
      vsprintf(env->term_buf, fmt, arg);
      /* pass the output to the user-defined routine */
      if (env->term_hook != NULL)
      {  if (env->term_hook(env->term_info, env->term_buf) != 0)
            goto skip;
      }
#if 0
      /* send the output to the terminal */
      fputs(env->term_buf, stdout);
      fflush(stdout);
      /* copy the output to the text file */
      if (env->tee_file != NULL)
      {  fputs(env->term_buf, env->tee_file);
         fflush(env->tee_file);
      }
#endif
skip: return;
}

/***********************************************************************
*  NAME
*
*  glp_term_out - enable/disable terminal output
*
*  SYNOPSIS
*
*  int glp_term_out(int flag);
*
*  DESCRIPTION
*
*  Depending on the parameter flag the routine glp_term_out enables or
*  disables terminal output performed by glpk routines:
*
*  GLP_ON  - enable terminal output;
*  GLP_OFF - disable terminal output.
*
*  RETURNS
*
*  The routine glp_term_out returns the previous value of the terminal
*  output flag. */

int glp_term_out(int flag)
{     ENV *env = get_env_ptr();
      int old = env->term_out;
      if (!(flag == GLP_ON || flag == GLP_OFF))
         xerror("glp_term_out: flag = %d; invalid value\n", flag);
      env->term_out = flag;
      return old;
}

/***********************************************************************
*  NAME
*
*  glp_term_hook - install hook to intercept terminal output
*
*  SYNOPSIS
*
*  void glp_term_hook(int (*func)(void *info, const char *s),
*     void *info);
*
*  DESCRIPTION
*
*  The routine glp_term_hook installs a user-defined hook routine to
*  intercept all terminal output performed by glpk routines.
*
*  This feature can be used to redirect the terminal output to other
*  destination, for example to a file or a text window.
*
*  The parameter func specifies the user-defined hook routine. It is
*  called from an internal printing routine, which passes to it two
*  parameters: info and s. The parameter info is a transit pointer,
*  specified in the corresponding call to the routine glp_term_hook;
*  it may be used to pass some information to the hook routine. The
*  parameter s is a pointer to the null terminated character string,
*  which is intended to be written to the terminal. If the hook routine
*  returns zero, the printing routine writes the string s to the
*  terminal in a usual way; otherwise, if the hook routine returns
*  non-zero, no terminal output is performed.
*
*  To uninstall the hook routine the parameters func and info should be
*  specified as NULL. */

void glp_term_hook(int (*func)(void *info, const char *s), void *info)
{     ENV *env = get_env_ptr();
      if (func == NULL)
      {  env->term_hook = NULL;
         env->term_info = NULL;
      }
      else
      {  env->term_hook = func;
         env->term_info = info;
      }
      return;
}

/***********************************************************************
*  NAME
*
*  glp_open_tee - start copying terminal output to text file
*
*  SYNOPSIS
*
*  int glp_open_tee(const char *fname);
*
*  DESCRIPTION
*
*  The routine glp_open_tee starts copying all the terminal output to
*  an output text file, whose name is specified by the character string
*  fname.
*
*  RETURNS
*
*  0 - operation successful
*  1 - copying terminal output is already active
*  2 - unable to create output file */

int glp_open_tee(const char *fname)
{     ENV *env = get_env_ptr();
      if (env->tee_file != NULL)
      {  /* copying terminal output is already active */
         return 1;
      }
      env->tee_file = fopen(fname, "w");
      if (env->tee_file == NULL)
      {  /* unable to create output file */
         return 2;
      }
      return 0;
}

/***********************************************************************
*  NAME
*
*  glp_close_tee - stop copying terminal output to text file
*
*  SYNOPSIS
*
*  int glp_close_tee(void);
*
*  DESCRIPTION
*
*  The routine glp_close_tee stops copying the terminal output to the
*  output text file previously open by the routine glp_open_tee closing
*  that file.
*
*  RETURNS
*
*  0 - operation successful
*  1 - copying terminal output was not started */

int glp_close_tee(void)
{     ENV *env = get_env_ptr();
      if (env->tee_file == NULL)
      {  /* copying terminal output was not started */
         return 1;
      }
      fclose(env->tee_file);
      env->tee_file = NULL;
      return 0;
}

/* eof */
