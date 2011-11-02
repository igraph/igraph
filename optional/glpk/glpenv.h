/* glpenv.h (GLPK environment) */

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

#ifndef GLPENV_H
#define GLPENV_H

#include "glpstd.h"
#include "glplib.h"

typedef struct ENV ENV;
typedef struct MEM MEM;
typedef struct XFILE XFILE;

#define ENV_MAGIC 0x454E5631
/* environment block magic value */

#define TERM_BUF_SIZE 4096
/* terminal output buffer size, in bytes */

#define IOERR_MSG_SIZE 1024
/* i/o error message buffer size, in bytes */

#define MEM_MAGIC 0x4D454D31
/* memory block descriptor magic value */

struct ENV
{     /* environment block */
      int magic;
      /* magic value used for debugging */
      char version[7+1];
      /* version string returned by the routine glp_version */
      /*--------------------------------------------------------------*/
      /* terminal output */
      char *term_buf; /* char term_buf[TERM_BUF_SIZE]; */
      /* terminal output buffer */
      int term_out;
      /* flag to enable/disable terminal output */
      int (*term_hook)(void *info, const char *s);
      /* user-defined routine to intercept terminal output */
      void *term_info;
      /* transit pointer (cookie) passed to the routine term_hook */
      FILE *tee_file;
      /* output stream used to copy terminal output */
      /*--------------------------------------------------------------*/
      /* error handling */
      const char *err_file;
      /* value of the __FILE__ macro passed to glp_error */
      int err_line;
      /* value of the __LINE__ macro passed to glp_error */
      void (*err_hook)(void *info);
      /* user-defined routine to intercept abnormal termination */
      void *err_info;
      /* transit pointer (cookie) passed to the routine err_hook */
      /*--------------------------------------------------------------*/
      /* memory allocation */
      glp_long mem_limit;
      /* maximal amount of memory (in bytes) available for dynamic
         allocation */
      MEM *mem_ptr;
      /* pointer to the linked list of allocated memory blocks */
      int mem_count;
      /* total number of currently allocated memory blocks */
      int mem_cpeak;
      /* peak value of mem_count */
      glp_long mem_total;
      /* total amount of currently allocated memory (in bytes; is the
         sum of the size field over all memory block descriptors) */
      glp_long mem_tpeak;
      /* peak value of mem_total */
      /*--------------------------------------------------------------*/
      /* stream input/output */
      XFILE *file_ptr;
      /* pointer to the linked list of active stream descriptors */
      char *ioerr_msg; /* char ioerr_msg[IOERR_MSG_SIZE]; */
      /* input/output error message buffer */
      /*--------------------------------------------------------------*/
      /* shared libraries support */
      void *h_odbc;
      /* handle to ODBC shared library */
      void *h_mysql;
      /* handle to MySQL shared library */
};

struct MEM
{     /* memory block descriptor */
      int flag;
      /* descriptor flag */
      int size;
      /* size of block (in bytes, including descriptor) */
      MEM *prev;
      /* pointer to previous memory block descriptor */
      MEM *next;
      /* pointer to next memory block descriptor */
};

struct XFILE
{     /* input/output stream descriptor */
      int type;
      /* stream handle type: */
#define FH_FILE   0x11  /* FILE   */
#define FH_ZLIB   0x22  /* gzFile */
      void *fh;
      /* pointer to stream handle */
      XFILE *prev;
      /* pointer to previous stream descriptor */
      XFILE *next;
      /* pointer to next stream descriptor */
};

#define XEOF (-1)

#define get_env_ptr _glp_get_env_ptr
ENV *get_env_ptr(void);
/* retrieve pointer to environment block */

#define tls_set_ptr _glp_tls_set_ptr
void tls_set_ptr(void *ptr);
/* store global pointer in TLS */

#define tls_get_ptr _glp_tls_get_ptr
void *tls_get_ptr(void);
/* retrieve global pointer from TLS */

#define xprintf glp_printf
void glp_printf(const char *fmt, ...);
/* write formatted output to the terminal */

#define xvprintf glp_vprintf
void glp_vprintf(const char *fmt, va_list arg);
/* write formatted output to the terminal */

#ifndef GLP_ERROR_DEFINED
#define GLP_ERROR_DEFINED
typedef void (*_glp_error)(const char *fmt, ...);
#endif

#define xerror glp_error_(__FILE__, __LINE__)
_glp_error glp_error_(const char *file, int line);
/* display error message and terminate execution */

#define xassert(expr) \
      ((void)((expr) || (glp_assert_(#expr, __FILE__, __LINE__), 1)))
void glp_assert_(const char *expr, const char *file, int line);
/* check for logical condition */

#define xmalloc glp_malloc
void *glp_malloc(int size);
/* allocate memory block */

#define xcalloc glp_calloc
void *glp_calloc(int n, int size);
/* allocate memory block */

#define xfree glp_free
void glp_free(void *ptr);
/* free memory block */

#define xtime glp_time
glp_long glp_time(void);
/* determine current universal time */

#define xdifftime glp_difftime
double glp_difftime(glp_long t1, glp_long t0);
/* compute difference between two time values, in seconds */

#define lib_err_msg _glp_lib_err_msg
void lib_err_msg(const char *msg);

#define xerrmsg _glp_lib_xerrmsg
const char *xerrmsg(void);

#define xfopen _glp_lib_xfopen
XFILE *xfopen(const char *fname, const char *mode);

#define xferror _glp_lib_xferror
int xferror(XFILE *file);

#define xfeof _glp_lib_xfeof
int xfeof(XFILE *file);

#define xfgetc _glp_lib_xfgetc
int xfgetc(XFILE *file);

#define xfputc _glp_lib_xfputc
int xfputc(int c, XFILE *file);

#define xfflush _glp_lib_xfflush
int xfflush(XFILE *fp);

#define xfclose _glp_lib_xfclose
int xfclose(XFILE *file);

#define xfprintf _glp_lib_xfprintf
int xfprintf(XFILE *file, const char *fmt, ...);

#define xdlopen _glp_xdlopen
void *xdlopen(const char *module);

#define xdlsym _glp_xdlsym
void *xdlsym(void *h, const char *symbol);

#define xdlclose _glp_xdlclose
void xdlclose(void *h);

#endif

/* eof */
