/* glpenv04.c (error handling) */

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

#include "glpapi.h"

/***********************************************************************
*  NAME
*
*  glp_error - display error message and terminate execution
*
*  SYNOPSIS
*
*  void glp_error(const char *fmt, ...);
*
*  DESCRIPTION
*
*  The routine glp_error (implemented as a macro) formats its
*  parameters using the format control string fmt, writes the formatted
*  message to the terminal, and abnormally terminates the program. */

static void error(const char *fmt, ...)
{     ENV *env = get_env_ptr();
      va_list arg;
      env->term_out = GLP_ON;
      va_start(arg, fmt);
      xvprintf(fmt, arg);
      va_end(arg);
      xprintf("Error detected in file %s at line %d\n", env->err_file,
         env->err_line);
      if (env->err_hook != NULL)
         env->err_hook(env->err_info);
      abort();
      exit(EXIT_FAILURE);
      /* no return */
}

_glp_error glp_error_(const char *file, int line)
{     ENV *env = get_env_ptr();
      env->err_file = file;
      env->err_line = line;
      return error;
}

/***********************************************************************
*  NAME
*
*  glp_assert - check for logical condition
*
*  SYNOPSIS
*
*  #include "glplib.h"
*  void glp_assert(int expr);
*
*  DESCRIPTION
*
*  The routine glp_assert (implemented as a macro) checks for a logical
*  condition specified by the parameter expr. If the condition is false
*  (i.e. the value of expr is zero), the routine writes a message to
*  the terminal and abnormally terminates the program. */

void glp_assert_(const char *expr, const char *file, int line)
{     glp_error_(file, line)("Assertion failed: %s\n", expr);
      /* no return */
}

/***********************************************************************
*  NAME
*
*  glp_error_hook - install hook to intercept abnormal termination
*
*  SYNOPSIS
*
*  void glp_error_hook(void (*func)(void *info), void *info);
*
*  DESCRIPTION
*
*  The routine glp_error_hook installs a user-defined hook routine to
*  intercept abnormal termination.
*
*  The parameter func specifies the user-defined hook routine. It is
*  called from the routine glp_error before the latter calls the abort
*  function to abnormally terminate the application program because of
*  fatal error. The parameter info is a transit pointer, specified in
*  the corresponding call to the routine glp_error_hook; it may be used
*  to pass some information to the hook routine.
*
*  To uninstall the hook routine the parameters func and info should be
*  specified as NULL. */

void glp_error_hook(void (*func)(void *info), void *info)
{     ENV *env = get_env_ptr();
      if (func == NULL)
      {  env->err_hook = NULL;
         env->err_info = NULL;
      }
      else
      {  env->err_hook = func;
         env->err_info = info;
      }
      return;
}

/* eof */
