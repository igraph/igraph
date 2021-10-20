/* error.c (error handling) */

/***********************************************************************
*  This code is part of GLPK (GNU Linear Programming Kit).
*  Copyright (C) 2000-2015 Free Software Foundation, Inc.
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

/***********************************************************************
*  NAME
*
*  glp_error - display fatal error message and terminate execution
*
*  SYNOPSIS
*
*  void glp_error(const char *fmt, ...);
*
*  DESCRIPTION
*
*  The routine glp_error (implemented as a macro) formats its
*  parameters using the format control string fmt, writes the formatted
*  message on the terminal, and abnormally terminates the program. */

static void errfunc(const char *fmt, ...)
{     ENV *env = get_env_ptr();
      va_list arg;
#if 1 /* 07/XI-2015 */
      env->err_st = 1;
#endif
      env->term_out = GLP_ON;
      va_start(arg, fmt);
      xvprintf(fmt, arg);
      va_end(arg);
      xprintf("Error detected in file %s at line %d\n",
         env->err_file, env->err_line);
      if (env->err_hook != NULL)
         env->err_hook(env->err_info);
      abort();
      exit(EXIT_FAILURE);
      /* no return */
}

glp_errfunc glp_error_(const char *file, int line)
{     ENV *env = get_env_ptr();
      env->err_file = file;
      env->err_line = line;
      return errfunc;
}

#if 1 /* 07/XI-2015 */
/***********************************************************************
*  NAME
*
*  glp_at_error - check for error state
*
*  SYNOPSIS
*
*  int glp_at_error(void);
*
*  DESCRIPTION
*
*  The routine glp_at_error checks if the GLPK environment is at error
*  state, i.e. if the call to the routine is (indirectly) made from the
*  glp_error routine via an user-defined hook routine.
*
*  RETURNS
*
*  If the GLPK environment is at error state, the routine glp_at_error
*  returns non-zero, otherwise zero. */

int glp_at_error(void)
{     ENV *env = get_env_ptr();
      return env->err_st;
}
#endif

/***********************************************************************
*  NAME
*
*  glp_assert - check for logical condition
*
*  SYNOPSIS
*
*  void glp_assert(int expr);
*
*  DESCRIPTION
*
*  The routine glp_assert (implemented as a macro) checks for a logical
*  condition specified by the parameter expr. If the condition is false
*  (i.e. the value of expr is zero), the routine writes a message on
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
*  both specified as NULL. */

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

/***********************************************************************
*  NAME
*
*  put_err_msg - provide error message string
*
*  SYNOPSIS
*
*  #include "env.h"
*  void put_err_msg(const char *msg);
*
*  DESCRIPTION
*
*  The routine put_err_msg stores an error message string pointed to by
*  msg to the environment block. */

void put_err_msg(const char *msg)
{     ENV *env = get_env_ptr();
      int len;
      len = strlen(msg);
      if (len >= EBUF_SIZE)
         len = EBUF_SIZE - 1;
      memcpy(env->err_buf, msg, len);
      if (len > 0 && env->err_buf[len-1] == '\n')
         len--;
      env->err_buf[len] = '\0';
      return;
}

/***********************************************************************
*  NAME
*
*  get_err_msg - obtain error message string
*
*  SYNOPSIS
*
*  #include "env.h"
*  const char *get_err_msg(void);
*
*  RETURNS
*
*  The routine get_err_msg returns a pointer to an error message string
*  previously stored by the routine put_err_msg. */

const char *get_err_msg(void)
{     ENV *env = get_env_ptr();
      return env->err_buf;
}

/* eof */
