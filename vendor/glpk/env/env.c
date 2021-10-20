/* env.c (GLPK environment initialization/termination) */

/***********************************************************************
*  This code is part of GLPK (GNU Linear Programming Kit).
*  Copyright (C) 2000-2017 Free Software Foundation, Inc.
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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "glpk.h"
#include "env.h"

/***********************************************************************
*  NAME
*
*  glp_init_env - initialize GLPK environment
*
*  SYNOPSIS
*
*  int glp_init_env(void);
*
*  DESCRIPTION
*
*  The routine glp_init_env initializes the GLPK environment. Normally
*  the application program does not need to call this routine, because
*  it is called automatically on the first call to any API routine.
*
*  RETURNS
*
*  The routine glp_init_env returns one of the following codes:
*
*  0 - initialization successful;
*  1 - environment has been already initialized;
*  2 - initialization failed (insufficient memory);
*  3 - initialization failed (unsupported programming model). */

int glp_init_env(void)
{     ENV *env;
      int ok;
      /* check if the programming model is supported */
      ok = (CHAR_BIT == 8 && sizeof(char) == 1 &&
         sizeof(short) == 2 && sizeof(int) == 4 &&
         (sizeof(void *) == 4 || sizeof(void *) == 8));
      if (!ok)
         return 3;
      /* check if the environment is already initialized */
      if (tls_get_ptr() != NULL)
         return 1;
      /* allocate and initialize the environment block */
      env = malloc(sizeof(ENV));
      if (env == NULL)
         return 2;
      memset(env, 0, sizeof(ENV));
#if 0 /* 14/I-2017 */
      sprintf(env->version, "%d.%d",
         GLP_MAJOR_VERSION, GLP_MINOR_VERSION);
#endif
      env->self = env;
      env->term_buf = malloc(TBUF_SIZE);
      if (env->term_buf == NULL)
      {  free(env);
         return 2;
      }
      env->term_out = GLP_ON;
      env->term_hook = NULL;
      env->term_info = NULL;
      env->tee_file = NULL;
#if 1 /* 23/XI-2015 */
      env->err_st = 0;
#endif
      env->err_file = NULL;
      env->err_line = 0;
      env->err_hook = NULL;
      env->err_info = NULL;
      env->err_buf = malloc(EBUF_SIZE);
      if (env->err_buf == NULL)
      {  free(env->term_buf);
         free(env);
         return 2;
      }
      env->err_buf[0] = '\0';
      env->mem_limit = SIZE_T_MAX;
      env->mem_ptr = NULL;
      env->mem_count = env->mem_cpeak = 0;
      env->mem_total = env->mem_tpeak = 0;
#if 1 /* 23/XI-2015 */
      env->gmp_pool = NULL;
      env->gmp_size = 0;
      env->gmp_work = NULL;
#endif
      env->h_odbc = env->h_mysql = NULL;
      /* save pointer to the environment block */
      tls_set_ptr(env);
      /* initialization successful */
      return 0;
}

/***********************************************************************
*  NAME
*
*  get_env_ptr - retrieve pointer to environment block
*
*  SYNOPSIS
*
*  #include "env.h"
*  ENV *get_env_ptr(void);
*
*  DESCRIPTION
*
*  The routine get_env_ptr retrieves and returns a pointer to the GLPK
*  environment block.
*
*  If the GLPK environment has not been initialized yet, the routine
*  performs initialization. If initialization fails, the routine prints
*  an error message to stderr and terminates the program.
*
*  RETURNS
*
*  The routine returns a pointer to the environment block. */

ENV *get_env_ptr(void)
{     ENV *env = tls_get_ptr();
      /* check if the environment has been initialized */
      if (env == NULL)
      {  /* not initialized yet; perform initialization */
         if (glp_init_env() != 0)
         {  /* initialization failed; display an error message */
            fprintf(stderr, "GLPK initialization failed\n");
            fflush(stderr);
            /* and abnormally terminate the program */
            abort();
         }
         /* initialization successful; retrieve the pointer */
         env = tls_get_ptr();
      }
      /* check if the environment block is valid */
      if (env->self != env)
      {  fprintf(stderr, "Invalid GLPK environment\n");
         fflush(stderr);
         abort();
      }
      return env;
}

/***********************************************************************
*  NAME
*
*  glp_version - determine library version
*
*  SYNOPSIS
*
*  const char *glp_version(void);
*
*  RETURNS
*
*  The routine glp_version returns a pointer to a null-terminated
*  character string, which specifies the version of the GLPK library in
*  the form "X.Y", where X is the major version number, and Y is the
*  minor version number, for example, "4.16". */

#define str(s) # s
#define xstr(s) str(s)

const char *glp_version(void)
#if 0 /* 14/I-2017 */
{     ENV *env = get_env_ptr();
      return env->version;
}
#else /* suggested by Heinrich */
{     return
         xstr(GLP_MAJOR_VERSION) "." xstr(GLP_MINOR_VERSION);
}
#endif

/***********************************************************************
*  NAME
*
*  glp_config - determine library configuration
*
*  SYNOPSIS
*
*  const char *glp_config(const char *option);
*
*  DESCRIPTION
*
*  The routine glp_config determines some options which were specified
*  on configuring the GLPK library.
*
*  RETURNS
*
*  The routine glp_config returns a pointer to a null-terminating
*  string depending on the option inquired.
*
*  For option = "TLS" the routine returns the thread local storage
*  class specifier used (e.g. "_Thread_local") if the GLPK library was
*  configured to run in multi-threaded environment, or NULL otherwise.
*
*  For option = "ODBC_DLNAME" the routine returns the name of ODBC
*  shared library if this option was enabled, or NULL otherwise.
*
*  For option = "MYSQL_DLNAME" the routine returns the name of MySQL
*  shared library if this option was enabled, or NULL otherwise. */

const char *glp_config(const char *option)
{     const char *s;
      if (strcmp(option, "TLS") == 0)
#ifndef TLS
         s = NULL;
#else
         s = xstr(TLS);
#endif
      else if (strcmp(option, "ODBC_DLNAME") == 0)
#ifndef ODBC_DLNAME
         s = NULL;
#else
         s = ODBC_DLNAME;
#endif
      else if (strcmp(option, "MYSQL_DLNAME") == 0)
#ifndef MYSQL_DLNAME
         s = NULL;
#else
         s = MYSQL_DLNAME;
#endif
      else
      {  /* invalid option is always disabled */
         s = NULL;
      }
      return s;
}

/***********************************************************************
*  NAME
*
*  glp_free_env - free GLPK environment
*
*  SYNOPSIS
*
*  int glp_free_env(void);
*
*  DESCRIPTION
*
*  The routine glp_free_env frees all resources used by GLPK routines
*  (memory blocks, etc.) which are currently still in use.
*
*  Normally the application program does not need to call this routine,
*  because GLPK routines always free all unused resources. However, if
*  the application program even has deleted all problem objects, there
*  will be several memory blocks still allocated for the library needs.
*  For some reasons the application program may want GLPK to free this
*  memory, in which case it should call glp_free_env.
*
*  Note that a call to glp_free_env invalidates all problem objects as
*  if no GLPK routine were called.
*
*  RETURNS
*
*  0 - termination successful;
*  1 - environment is inactive (was not initialized). */

int glp_free_env(void)
{     ENV *env = tls_get_ptr();
      MBD *desc;
      /* check if the environment is active */
      if (env == NULL)
         return 1;
      /* check if the environment block is valid */
      if (env->self != env)
      {  fprintf(stderr, "Invalid GLPK environment\n");
         fflush(stderr);
         abort();
      }
      /* close handles to shared libraries */
      if (env->h_odbc != NULL)
         xdlclose(env->h_odbc);
      if (env->h_mysql != NULL)
         xdlclose(env->h_mysql);
      /* free memory blocks which are still allocated */
      while (env->mem_ptr != NULL)
      {  desc = env->mem_ptr;
         env->mem_ptr = desc->next;
         free(desc);
      }
      /* close text file used for copying terminal output */
      if (env->tee_file != NULL)
         fclose(env->tee_file);
      /* invalidate the environment block */
      env->self = NULL;
      /* free memory allocated to the environment block */
      free(env->term_buf);
      free(env->err_buf);
      free(env);
      /* reset a pointer to the environment block */
      tls_set_ptr(NULL);
      /* termination successful */
      return 0;
}

/* eof */
