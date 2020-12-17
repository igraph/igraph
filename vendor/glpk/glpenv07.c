/* glpenv07.c (stream input/output) */

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
#pragma clang diagnostic ignored "-Wshorten-64-to-32"
#pragma clang diagnostic ignored "-Wsometimes-uninitialized"
#endif

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "glpenv.h"

/***********************************************************************
*  NAME
*
*  lib_err_msg - save error message string
*
*  SYNOPSIS
*
*  #include "glpenv.h"
*  void lib_err_msg(const char *msg);
*
*  DESCRIPTION
*
*  The routine lib_err_msg saves an error message string specified by
*  the parameter msg. The message is obtained by some library routines
*  with a call to strerror(errno). */

void lib_err_msg(const char *msg)
{     ENV *env = get_env_ptr();
      int len = strlen(msg);
      if (len >= IOERR_MSG_SIZE)
         len = IOERR_MSG_SIZE - 1;
      memcpy(env->ioerr_msg, msg, len);
      if (len > 0 && env->ioerr_msg[len-1] == '\n') len--;
      env->ioerr_msg[len] = '\0';
      return;
}

/***********************************************************************
*  NAME
*
*  xerrmsg - retrieve error message string
*
*  SYNOPSIS
*
*  #include "glpenv.h"
*  const char *xerrmsg(void);
*
*  RETURNS
*
*  The routine xerrmsg returns a pointer to an error message string
*  previously set by some library routine to indicate an error. */

const char *xerrmsg(void)
{     ENV *env = get_env_ptr();
      return env->ioerr_msg;
}

/***********************************************************************
*  NAME
*
*  xfopen - open a stream
*
*  SYNOPSIS
*
*  #include "glpenv.h"
*  XFILE *xfopen(const char *fname, const char *mode);
*
*  DESCRIPTION
*
*  The routine xfopen opens the file whose name is a string pointed to
*  by fname and associates a stream with it.
*
*  The parameter mode points to a string, which indicates the open mode
*  and should be one of the following:
*
*  "r"   open text file for reading;
*  "w"   truncate to zero length or create text file for writing;
*  "rb"  open binary file for reading;
*  "wb"  truncate to zero length or create binary file for writing.
*
*  RETURNS
*
*  The routine xfopen returns a pointer to the object controlling the
*  stream. If the open operation fails, xfopen returns NULL. */

static void *c_fopen(const char *fname, const char *mode);
static void *z_fopen(const char *fname, const char *mode);

static int is_gz_file(const char *fname)
{     char *ext = strrchr(fname, '.');
      return ext != NULL && strcmp(ext, ".gz") == 0;
}

XFILE *xfopen(const char *fname, const char *mode)
{     ENV *env = get_env_ptr();
      XFILE *fp;
      int type;
      void *fh;
      if (!is_gz_file(fname))
      {  type = FH_FILE;
         fh = c_fopen(fname, mode);
      }
      else
      {  type = FH_ZLIB;
         fh = z_fopen(fname, mode);
      }
      if (fh == NULL)
      {  fp = NULL;
         goto done;
      }
      fp = xmalloc(sizeof(XFILE));
      fp->type = type;
      fp->fh = fh;
      fp->prev = NULL;
      fp->next = env->file_ptr;
      if (fp->next != NULL) fp->next->prev = fp;
      env->file_ptr = fp;
done: return fp;
}

/***********************************************************************
*  NAME
*
*  xfgetc - read character from the stream
*
*  SYNOPSIS
*
*  #include "glpenv.h"
*  int xfgetc(XFILE *fp);
*
*  DESCRIPTION
*
*  If the end-of-file indicator for the input stream pointed to by fp
*  is not set and a next character is present, the routine xfgetc
*  obtains that character as an unsigned char converted to an int and
*  advances the associated file position indicator for the stream (if
*  defined).
*
*  RETURNS
*
*  If the end-of-file indicator for the stream is set, or if the
*  stream is at end-of-file, the end-of-file indicator for the stream
*  is set and the routine xfgetc returns XEOF. Otherwise, the routine
*  xfgetc returns the next character from the input stream pointed to
*  by fp. If a read error occurs, the error indicator for the stream is
*  set and the xfgetc routine returns XEOF.
*
*  Note: An end-of-file and a read error can be distinguished by use of
*  the routines xfeof and xferror. */

static int c_fgetc(void *fh);
static int z_fgetc(void *fh);

int xfgetc(XFILE *fp)
{     int c;
      switch (fp->type)
      {  case FH_FILE:
            c = c_fgetc(fp->fh);
            break;
         case FH_ZLIB:
            c = z_fgetc(fp->fh);
            break;
         default:
            xassert(fp != fp);
      }
      return c;
}

/***********************************************************************
*  NAME
*
*  xfputc - write character to the stream
*
*  SYNOPSIS
*
*  #include "glpenv.h"
*  int xfputc(int c, XFILE *fp);
*
*  DESCRIPTION
*
*  The routine xfputc writes the character specified by c (converted
*  to an unsigned char) to the output stream pointed to by fp, at the
*  position indicated by the associated file position indicator (if
*  defined), and advances the indicator appropriately.
*
*  RETURNS
*
*  The routine xfputc returns the character written. If a write error
*  occurs, the error indicator for the stream is set and xfputc returns
*  XEOF. */

static int c_fputc(int c, void *fh);
static int z_fputc(int c, void *fh);

int xfputc(int c, XFILE *fp)
{     switch (fp->type)
      {  case FH_FILE:
            c = c_fputc(c, fp->fh);
            break;
         case FH_ZLIB:
            c = z_fputc(c, fp->fh);
            break;
         default:
            xassert(fp != fp);
      }
      return c;
}

/***********************************************************************
*  NAME
*
*  xferror - test error indicator for the stream
*
*  SYNOPSIS
*
*  #include "glpenv.h"
*  int xferror(XFILE *fp);
*
*  DESCRIPTION
*
*  The routine xferror tests the error indicator for the stream
*  pointed to by fp.
*
*  RETURNS
*
*  The routine xferror returns non-zero if and only if the error
*  indicator is set for the stream. */

static int c_ferror(void *fh);
static int z_ferror(void *fh);

int xferror(XFILE *fp)
{     int ret;
      switch (fp->type)
      {  case FH_FILE:
            ret = c_ferror(fp->fh);
            break;
         case FH_ZLIB:
            ret = z_ferror(fp->fh);
            break;
         default:
            xassert(fp != fp);
      }
      return ret;
}

/***********************************************************************
*  NAME
*
*  xfeof - test end-of-file indicator for the stream
*
*  SYNOPSIS
*
*  #include "glpenv.h"
*  int xfeof(XFILE *fp);
*
*  DESCRIPTION
*
*  The routine xfeof tests the end-of-file indicator for the stream
*  pointed to by fp.
*
*  RETURNS
*
*  The routine xfeof returns non-zero if and only if the end-of-file
*  indicator is set for the stream. */

static int c_feof(void *fh);
static int z_feof(void *fh);

int xfeof(XFILE *fp)
{     int ret;
      switch (fp->type)
      {  case FH_FILE:
            ret = c_feof(fp->fh);
            break;
         case FH_ZLIB:
            ret = z_feof(fp->fh);
            break;
         default:
            xassert(fp != fp);
      }
      return ret;
}

int xfprintf(XFILE *file, const char *fmt, ...)
{     ENV *env = get_env_ptr();
      int cnt, j;
      va_list arg;
      va_start(arg, fmt);
      cnt = vsprintf(env->term_buf, fmt, arg);
      va_end(arg);
      for (j = 0; j < cnt; j++)
      {  if (xfputc(env->term_buf[j], file) < 0)
         {  cnt = -1;
            break;
         }
      }
      return cnt;
}

/***********************************************************************
*  NAME
*
*  xfflush - flush the stream
*
*  SYNOPSIS
*
*  #include "glpenv.h"
*  int xfflush(XFILE *fp);
*
*  DESCRIPTION
*
*  The routine xfflush causes any unwritten data for the output stream
*  pointed to by fp to be written to the associated file.
*
*  RETURNS
*
*  The routine xfflush returns zero if the stream was successfully
*  flushed. Otherwise, xfflush sets the error indicator for the stream
*  and returns XEOF. */

static int c_fflush(void *fh);
static int z_fflush(void *fh);

int xfflush(XFILE *fp)
{     int ret;
      switch (fp->type)
      {  case FH_FILE:
            ret = c_fflush(fp->fh);
            break;
         case FH_ZLIB:
            ret = z_fflush(fp->fh);
            break;
         default:
            xassert(fp != fp);
      }
      return ret;
}

/***********************************************************************
*  NAME
*
*  xfclose - close the stream
*
*  SYNOPSIS
*
*  #include "glpenv.h"
*  int xfclose(XFILE *fp);
*
*  DESCRIPTION
*
*  A successful call to the routine xfclose causes the stream pointed
*  to by fp to be flushed and the associated file to be closed. Whether
*  or not the call succeeds, the stream is disassociated from the file.
*
*  RETURNS
*
*  The routine xfclose returns zero if the stream was successfully
*  closed, or XEOF if any errors were detected. */

static int c_fclose(void *fh);
static int z_fclose(void *fh);

int xfclose(XFILE *fp)
{     ENV *env = get_env_ptr();
      int ret;
      switch (fp->type)
      {  case FH_FILE:
            ret = c_fclose(fp->fh);
            break;
         case FH_ZLIB:
            ret = z_fclose(fp->fh);
            break;
         default:
            xassert(fp != fp);
      }
      fp->type = 0xF00BAD;
      if (fp->prev == NULL)
         env->file_ptr = fp->next;
      else
         fp->prev->next = fp->next;
      if (fp->next == NULL)
         ;
      else
         fp->next->prev = fp->prev;
      xfree(fp);
      return ret;
}

/***********************************************************************
*  The following routines implement stream input/output based on the
*  standard C streams. */

static void *c_fopen(const char *fname, const char *mode)
{     FILE *fh;
      /* if (strcmp(fname, "/dev/stdin") == 0) */
      /*    fh = stdin; */
      /* else if (strcmp(fname, "/dev/stdout") == 0) */
      /*    fh = stdout; */
      /* else if (strcmp(fname, "/dev/stderr") == 0) */
      /*    fh = stderr; */
      /* else */
         fh = fopen(fname, mode);
      if (fh == NULL)
         lib_err_msg(strerror(errno));
      return fh;
}

static int c_fgetc(void *_fh)
{     FILE *fh = _fh;
      int c;
      if (ferror(fh) || feof(fh))
      {  c = XEOF;
         goto done;
      }
      c = fgetc(fh);
      if (ferror(fh))
      {  lib_err_msg(strerror(errno));
         c = XEOF;
      }
      else if (feof(fh))
         c = XEOF;
      else
         xassert(0x00 <= c && c <= 0xFF);
done: return c;
}

static int c_fputc(int c, void *_fh)
{     FILE *fh = _fh;
      if (ferror(fh))
      {  c = XEOF;
         goto done;
      }
      c = (unsigned char)c;
      fputc(c, fh);
      if (ferror(fh))
      {  lib_err_msg(strerror(errno));
         c = XEOF;
      }
done: return c;
}

static int c_ferror(void *_fh)
{     FILE *fh = _fh;
      return ferror(fh);
}

static int c_feof(void *_fh)
{     FILE *fh = _fh;
      return feof(fh);
}

static int c_fflush(void *_fh)
{     FILE *fh = _fh;
      int ret;
      ret = fflush(fh);
      if (ret != 0)
      {  lib_err_msg(strerror(errno));
         ret = XEOF;
      }
      return ret;
}

static int c_fclose(void *_fh)
{     FILE *fh = _fh;
      int ret;
      /* if (fh == stdin) */
      /*    ret = 0; */
      /* else if (fh == stdout || fh == stderr) */
      /*    fflush(fh), ret = 0; */
      /* else */
         ret = fclose(fh);
      if (ret != 0)
      {  lib_err_msg(strerror(errno));
         ret = XEOF;
      }
      return ret;
}

/***********************************************************************
*  The following routines implement stream input/output based on the
*  zlib library, which provides processing .gz files "on the fly". */

#ifndef HAVE_ZLIB

static void *z_fopen(const char *fname, const char *mode)
{     xassert(fname == fname);
      xassert(mode == mode);
      lib_err_msg("Compressed files not supported");
      return NULL;
}

static int z_fgetc(void *fh)
{     xassert(fh != fh);
      return 0;
}

static int z_fputc(int c, void *fh)
{     xassert(c != c);
      xassert(fh != fh);
      return 0;
}

static int z_ferror(void *fh)
{     xassert(fh != fh);
      return 0;
}

static int z_feof(void *fh)
{     xassert(fh != fh);
      return 0;
}

static int z_fflush(void *fh)
{     xassert(fh != fh);
      return 0;
}

static int z_fclose(void *fh)
{     xassert(fh != fh);
      return 0;
}

#else

#include <zlib.h>

struct z_file
{     /* .gz file handle */
      gzFile file;
      /* pointer to .gz stream */
      int err;
      /* i/o error indicator */
      int eof;
      /* end-of-file indicator */
};

static void *z_fopen(const char *fname, const char *mode)
{     struct z_file *fh;
      gzFile file;
      if (strcmp(mode, "r") == 0 || strcmp(mode, "rb") == 0)
         mode = "rb";
      else if (strcmp(mode, "w") == 0 || strcmp(mode, "wb") == 0)
         mode = "wb";
      else
      {  lib_err_msg("Invalid open mode");
         fh = NULL;
         goto done;
      }
      file = gzopen(fname, mode);
      if (file == NULL)
      {  lib_err_msg(strerror(errno));
         fh = NULL;
         goto done;
      }
      fh = xmalloc(sizeof(struct z_file));
      fh->file = file;
      fh->err = fh->eof = 0;
done: return fh;
}

static int z_fgetc(void *_fh)
{     struct z_file *fh = _fh;
      int c;
      if (fh->err || fh->eof)
      {  c = XEOF;
         goto done;
      }
      c = gzgetc(fh->file);
      if (c < 0)
      {  int errnum;
         const char *msg;
         msg = gzerror(fh->file, &errnum);
         if (errnum == Z_STREAM_END)
            fh->eof = 1;
         else if (errnum == Z_ERRNO)
         {  fh->err = 1;
            lib_err_msg(strerror(errno));
         }
         else
         {  fh->err = 1;
            lib_err_msg(msg);
         }
         c = XEOF;
      }
      else
         xassert(0x00 <= c && c <= 0xFF);
done: return c;
}

static int z_fputc(int c, void *_fh)
{     struct z_file *fh = _fh;
      if (fh->err)
      {  c = XEOF;
         goto done;
      }
      c = (unsigned char)c;
      if (gzputc(fh->file, c) < 0)
      {  int errnum;
         const char *msg;
         fh->err = 1;
         msg = gzerror(fh->file, &errnum);
         if (errnum == Z_ERRNO)
            lib_err_msg(strerror(errno));
         else
            lib_err_msg(msg);
         c = XEOF;
      }
done: return c;
}

static int z_ferror(void *_fh)
{     struct z_file *fh = _fh;
      return fh->err;
}

static int z_feof(void *_fh)
{     struct z_file *fh = _fh;
      return fh->eof;
}

static int z_fflush(void *_fh)
{     struct z_file *fh = _fh;
      int ret;
      ret = gzflush(fh->file, Z_FINISH);
      if (ret == Z_OK)
         ret = 0;
      else
      {  int errnum;
         const char *msg;
         fh->err = 1;
         msg = gzerror(fh->file, &errnum);
         if (errnum == Z_ERRNO)
            lib_err_msg(strerror(errno));
         else
            lib_err_msg(msg);
         ret = XEOF;
      }
      return ret;
}

static int z_fclose(void *_fh)
{     struct z_file *fh = _fh;
      gzclose(fh->file);
      xfree(fh);
      return 0;
}

#endif

/* eof */
