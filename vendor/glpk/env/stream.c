/* stream.c (stream input/output) */

/***********************************************************************
*  This code is part of GLPK (GNU Linear Programming Kit).
*  Copyright (C) 2008-2017 Free Software Foundation, Inc.
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
#include "zlib.h"

struct glp_file
{     /* sequential stream descriptor */
      char *base;
      /* pointer to buffer */
      int size;
      /* size of buffer, in bytes */
      char *ptr;
      /* pointer to next byte in buffer */
      int cnt;
      /* count of bytes in buffer */
      int flag;
      /* stream flags: */
#define IONULL 0x01 /* null file */
#define IOSTD  0x02 /* standard stream */
#define IOGZIP 0x04 /* gzipped file */
#define IOWRT  0x08 /* output stream */
#define IOEOF  0x10 /* end of file */
#define IOERR  0x20 /* input/output error */
      void *file;
      /* pointer to underlying control object */
};

/***********************************************************************
*  NAME
*
*  glp_open - open stream
*
*  SYNOPSIS
*
*  glp_file *glp_open(const char *name, const char *mode);
*
*  DESCRIPTION
*
*  The routine glp_open opens a file whose name is a string pointed to
*  by name and associates a stream with it.
*
*  The following special filenames are recognized by the routine (this
*  feature is platform independent):
*
*  "/dev/null"    empty (null) file;
*  "/dev/stdin"   standard input stream;
*  "/dev/stdout"  standard output stream;
*  "/dev/stderr"  standard error stream.
*
*  If the specified filename is ended with ".gz", it is assumed that
*  the file is in gzipped format. In this case the file is compressed
*  or decompressed by the I/O routines "on the fly".
*
*  The parameter mode points to a string, which indicates the open mode
*  and should be one of the following:
*
*  "r"   open text file for reading;
*  "w"   truncate to zero length or create text file for writing;
*  "a"   append, open or create text file for writing at end-of-file;
*  "rb"  open binary file for reading;
*  "wb"  truncate to zero length or create binary file for writing;
*  "ab"  append, open or create binary file for writing at end-of-file.
*
*  RETURNS
*
*  The routine glp_open returns a pointer to the object controlling the
*  stream. If the operation fails, the routine returns NULL. */

glp_file *glp_open(const char *name, const char *mode)
{     glp_file *f;
      int flag;
      void *file;
      if (strcmp(mode, "r") == 0 || strcmp(mode, "rb") == 0)
         flag = 0;
      else if (strcmp(mode, "w") == 0 || strcmp(mode, "wb") == 0)
         flag = IOWRT;
#if 1 /* 08/V-2014 */
      else if (strcmp(mode, "a") == 0 || strcmp(mode, "ab") == 0)
         flag = IOWRT;
#endif
      else
         xerror("glp_open: invalid mode string\n");
      if (strcmp(name, "/dev/null") == 0)
      {  flag |= IONULL;
         file = NULL;
      }
      else if (strcmp(name, "/dev/stdin") == 0)
      {  flag |= IOSTD;
         file = stdin;
      }
      else if (strcmp(name, "/dev/stdout") == 0)
      {  flag |= IOSTD;
         file = stdout;
      }
      else if (strcmp(name, "/dev/stderr") == 0)
      {  flag |= IOSTD;
         file = stderr;
      }
      else
      {  char *ext = strrchr(name, '.');
         if (ext == NULL || strcmp(ext, ".gz") != 0)
         {  file = fopen(name, mode);
            if (file == NULL)
#if 0 /* 29/I-2017 */
            {  put_err_msg(strerror(errno));
#else
            {  put_err_msg(xstrerr(errno));
#endif
               return NULL;
            }
         }
         else
         {  flag |= IOGZIP;
            if (strcmp(mode, "r") == 0)
               mode = "rb";
            else if (strcmp(mode, "w") == 0)
               mode = "wb";
#if 1 /* 08/V-2014; this mode seems not to work */
            else if (strcmp(mode, "a") == 0)
               mode = "ab";
#endif
            file = gzopen(name, mode);
            if (file == NULL)
#if 0 /* 29/I-2017 */
            {  put_err_msg(strerror(errno));
#else
            {  put_err_msg(xstrerr(errno));
#endif
               return NULL;
            }
         }
      }
      f = talloc(1, glp_file);
      f->base = talloc(BUFSIZ, char);
      f->size = BUFSIZ;
      f->ptr = f->base;
      f->cnt = 0;
      f->flag = flag;
      f->file = file;
      return f;
}

/***********************************************************************
*  NAME
*
*  glp_eof - test end-of-file indicator
*
*  SYNOPSIS
*
*  int glp_eof(glp_file *f);
*
*  DESCRIPTION
*
*  The routine glp_eof tests the end-of-file indicator for the stream
*  pointed to by f.
*
*  RETURNS
*
*  The routine glp_eof returns non-zero if and only if the end-of-file
*  indicator is set for the specified stream. */

int glp_eof(glp_file *f)
{     return
         f->flag & IOEOF;
}

/***********************************************************************
*  NAME
*
*  glp_ioerr - test I/O error indicator
*
*  SYNOPSIS
*
*  int glp_ioerr(glp_file *f);
*
*  DESCRIPTION
*
*  The routine glp_ioerr tests the I/O error indicator for the stream
*  pointed to by f.
*
*  RETURNS
*
*  The routine glp_ioerr returns non-zero if and only if the I/O error
*  indicator is set for the specified stream. */

int glp_ioerr(glp_file *f)
{     return
         f->flag & IOERR;
}

/***********************************************************************
*  NAME
*
*  glp_read - read data from stream
*
*  SYNOPSIS
*
*  int glp_read(glp_file *f, void *buf, int nnn);
*
*  DESCRIPTION
*
*  The routine glp_read reads, into the buffer pointed to by buf, up to
*  nnn bytes, from the stream pointed to by f.
*
*  RETURNS
*
*  The routine glp_read returns the number of bytes successfully read
*  (which may be less than nnn). If an end-of-file is encountered, the
*  end-of-file indicator for the stream is set and glp_read returns
*  zero. If a read error occurs, the error indicator for the stream is
*  set and glp_read returns a negative value. */

int glp_read(glp_file *f, void *buf, int nnn)
{     int nrd, cnt;
      if (f->flag & IOWRT)
         xerror("glp_read: attempt to read from output stream\n");
      if (nnn < 1)
         xerror("glp_read: nnn = %d; invalid parameter\n", nnn);
      for (nrd = 0; nrd < nnn; nrd += cnt)
      {  if (f->cnt == 0)
         {  /* buffer is empty; fill it */
            if (f->flag & IONULL)
               cnt = 0;
            else if (!(f->flag & IOGZIP))
            {  cnt = fread(f->base, 1, f->size, (FILE *)(f->file));
               if (ferror((FILE *)(f->file)))
               {  f->flag |= IOERR;
#if 0 /* 29/I-2017 */
                  put_err_msg(strerror(errno));
#else
                  put_err_msg(xstrerr(errno));
#endif
                  return EOF;
               }
            }
            else
            {  int errnum;
               const char *msg;
               cnt = gzread((gzFile)(f->file), f->base, f->size);
               if (cnt < 0)
               {  f->flag |= IOERR;
                  msg = gzerror((gzFile)(f->file), &errnum);
                  if (errnum == Z_ERRNO)
#if 0 /* 29/I-2017 */
                     put_err_msg(strerror(errno));
#else
                     put_err_msg(xstrerr(errno));
#endif
                  else
                     put_err_msg(msg);
                  return EOF;
               }
            }
            if (cnt == 0)
            {  if (nrd == 0)
                  f->flag |= IOEOF;
               break;
            }
            f->ptr = f->base;
            f->cnt = cnt;
         }
         cnt = nnn - nrd;
         if (cnt > f->cnt)
            cnt = f->cnt;
         memcpy((char *)buf + nrd, f->ptr, cnt);
         f->ptr += cnt;
         f->cnt -= cnt;
      }
      return nrd;
}

/***********************************************************************
*  NAME
*
*  glp_getc - read character from stream
*
*  SYNOPSIS
*
*  int glp_getc(glp_file *f);
*
*  DESCRIPTION
*
*  The routine glp_getc obtains a next character as an unsigned char
*  converted to an int from the input stream pointed to by f.
*
*  RETURNS
*
*  The routine glp_getc returns the next character obtained. However,
*  if an end-of-file is encountered or a read error occurs, the routine
*  returns EOF. (An end-of-file and a read error can be distinguished
*  by use of the routines glp_eof and glp_ioerr.) */

int glp_getc(glp_file *f)
{     unsigned char buf[1];
      if (f->flag & IOWRT)
         xerror("glp_getc: attempt to read from output stream\n");
      if (glp_read(f, buf, 1) != 1)
         return EOF;
      return buf[0];
}

/***********************************************************************
*  do_flush - flush output stream
*
*  This routine causes buffered data for the specified output stream to
*  be written to the associated file.
*
*  If the operation was successful, the routine returns zero, otherwise
*  non-zero. */

static int do_flush(glp_file *f)
{     xassert(f->flag & IOWRT);
      if (f->cnt > 0)
      {  if (f->flag & IONULL)
            ;
         else if (!(f->flag & IOGZIP))
         {  if ((int)fwrite(f->base, 1, f->cnt, (FILE *)(f->file))
               != f->cnt)
            {  f->flag |= IOERR;
#if 0 /* 29/I-2017 */
               put_err_msg(strerror(errno));
#else
               put_err_msg(xstrerr(errno));
#endif
               return EOF;
            }
         }
         else
         {  int errnum;
            const char *msg;
            if (gzwrite((gzFile)(f->file), f->base, f->cnt) != f->cnt)
            {  f->flag |= IOERR;
               msg = gzerror((gzFile)(f->file), &errnum);
               if (errnum == Z_ERRNO)
#if 0 /* 29/I-2017 */
                  put_err_msg(strerror(errno));
#else
                  put_err_msg(xstrerr(errno));
#endif
               else
                  put_err_msg(msg);
               return EOF;
            }
         }
      }
      f->ptr = f->base;
      f->cnt = 0;
      return 0;
}

/***********************************************************************
*  NAME
*
*  glp_write - write data to stream
*
*  SYNOPSIS
*
*  int glp_write(glp_file *f, const void *buf, int nnn);
*
*  DESCRIPTION
*
*  The routine glp_write writes, from the buffer pointed to by buf, up
*  to nnn bytes, to the stream pointed to by f.
*
*  RETURNS
*
*  The routine glp_write returns the number of bytes successfully
*  written (which is equal to nnn). If a write error occurs, the error
*  indicator for the stream is set and glp_write returns a negative
*  value. */

int glp_write(glp_file *f, const void *buf, int nnn)
{     int nwr, cnt;
      if (!(f->flag & IOWRT))
         xerror("glp_write: attempt to write to input stream\n");
      if (nnn < 1)
         xerror("glp_write: nnn = %d; invalid parameter\n", nnn);
      for (nwr = 0; nwr < nnn; nwr += cnt)
      {  cnt = nnn - nwr;
         if (cnt > f->size - f->cnt)
            cnt = f->size - f->cnt;
         memcpy(f->ptr, (const char *)buf + nwr, cnt);
         f->ptr += cnt;
         f->cnt += cnt;
         if (f->cnt == f->size)
         {  /* buffer is full; flush it */
            if (do_flush(f) != 0)
               return EOF;
         }
      }
      return nwr;
}

/***********************************************************************
*  NAME
*
*  glp_format - write formatted data to stream
*
*  SYNOPSIS
*
*  int glp_format(glp_file *f, const char *fmt, ...);
*
*  DESCRIPTION
*
*  The routine glp_format writes formatted data to the stream pointed
*  to by f. The format control string pointed to by fmt specifies how
*  subsequent arguments are converted for output.
*
*  RETURNS
*
*  The routine glp_format returns the number of characters written, or
*  a negative value if an output error occurs. */

int glp_format(glp_file *f, const char *fmt, ...)
{     ENV *env = get_env_ptr();
      va_list arg;
      int nnn;
      if (!(f->flag & IOWRT))
         xerror("glp_format: attempt to write to input stream\n");
      va_start(arg, fmt);
      nnn = vsprintf(env->term_buf, fmt, arg);
      xassert(0 <= nnn && nnn < TBUF_SIZE);
      va_end(arg);
      return nnn == 0 ? 0 : glp_write(f, env->term_buf, nnn);
}

/***********************************************************************
*  NAME
*
*  glp_close - close stream
*
*  SYNOPSIS
*
*  int glp_close(glp_file *f);
*
*  DESCRIPTION
*
*  The routine glp_close closes the stream pointed to by f.
*
*  RETURNS
*
*  If the operation was successful, the routine returns zero, otherwise
*  non-zero. */

int glp_close(glp_file *f)
{     int ret = 0;
      if (f->flag & IOWRT)
      {  if (do_flush(f) != 0)
            ret = EOF;
      }
      if (f->flag & (IONULL | IOSTD))
         ;
      else if (!(f->flag & IOGZIP))
      {  if (fclose((FILE *)(f->file)) != 0)
         {  if (ret == 0)
#if 0 /* 29/I-2017 */
            {  put_err_msg(strerror(errno));
#else
            {  put_err_msg(xstrerr(errno));
#endif
               ret = EOF;
            }
         }
      }
      else
      {  int errnum;
         errnum = gzclose((gzFile)(f->file));
         if (errnum == Z_OK)
            ;
         else if (errnum == Z_ERRNO)
         {  if (ret == 0)
#if 0 /* 29/I-2017 */
            {  put_err_msg(strerror(errno));
#else
            {  put_err_msg(xstrerr(errno));
#endif
               ret = EOF;
            }
         }
#if 1 /* FIXME */
         else
         {  if (ret == 0)
            {  ENV *env = get_env_ptr();
               sprintf(env->term_buf, "gzclose returned %d", errnum);
               put_err_msg(env->term_buf);
               ret = EOF;
            }
         }
#endif
      }
      tfree(f->base);
      tfree(f);
      return ret;
}

/* eof */
