/* glpsdf.c (plain data file reading routines) */

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
#endif

#define GLPSDF_H

#define GLP_DATA_DEFINED
typedef struct glp_data glp_data;

#include "glpapi.h"

struct glp_data
{     /* plain data file */
      char *fname;
      /* name of data file */
      XFILE *fp;
      /* stream assigned to data file */
      void *jump; /* jmp_buf jump; */
      /* label for go to in case of error */
      int count;
      /* line count */
      int c;
      /* current character of XEOF */
      char item[255+1];
      /* current data item */
};

static void next_char(glp_data *data);

glp_data *glp_sdf_open_file(const char *fname)
{     /* open plain data file */
      glp_data *data = NULL;
      XFILE *fp;
      jmp_buf jump;
      fp = xfopen(fname, "r");
      if (fp == NULL)
      {  xprintf("Unable to open `%s' - %s\n", fname, xerrmsg());
         goto done;
      }
      data = xmalloc(sizeof(glp_data));
      data->fname = xmalloc(strlen(fname)+1);
      strcpy(data->fname, fname);
      data->fp = fp;
      data->jump = NULL;
      data->count = 0;
      data->c = '\n';
      data->item[0] = '\0';
      /* read the very first character */
      if (setjmp(jump))
      {  glp_sdf_close_file(data);
         data = NULL;
         goto done;
      }
      data->jump = jump;
      next_char(data);
      data->jump = NULL;
done: return data;
}

void glp_sdf_set_jump(glp_data *data, void *jump)
{     /* set up error handling */
      data->jump = jump;
      return;
}

void glp_sdf_error(glp_data *data, const char *fmt, ...)
{     /* print error message */
      va_list arg;
      xprintf("%s:%d: ", data->fname, data->count);
      va_start(arg, fmt);
      xvprintf(fmt, arg);
      va_end(arg);
      if (data->jump == NULL)
         xerror("");
      else
         longjmp(data->jump, 1);
      /* no return */
}

void glp_sdf_warning(glp_data *data, const char *fmt, ...)
{     /* print warning message */
      va_list arg;
      xprintf("%s:%d: warning: ", data->fname, data->count);
      va_start(arg, fmt);
      xvprintf(fmt, arg);
      va_end(arg);
      return;
}

static void next_char(glp_data *data)
{     /* read next character */
      int c;
      if (data->c == XEOF)
         glp_sdf_error(data, "unexpected end of file\n");
      else if (data->c == '\n')
         data->count++;
      c = xfgetc(data->fp);
      if (c < 0)
      {  if (xferror(data->fp))
            glp_sdf_error(data, "read error - %s\n", xerrmsg());
         else if (data->c == '\n')
            c = XEOF;
         else
         {  glp_sdf_warning(data, "missing final end of line\n");
            c = '\n';
         }
      }
      else if (c == '\n')
         ;
      else if (isspace(c))
         c = ' ';
      else if (iscntrl(c))
         glp_sdf_error(data, "invalid control character 0x%02X\n", c);
      data->c = c;
      return;
}

static void skip_pad(glp_data *data)
{     /* skip uninteresting characters and comments */
loop: while (data->c == ' ' || data->c == '\n')
         next_char(data);
      if (data->c == '/')
      {  next_char(data);
         if (data->c != '*')
            glp_sdf_error(data, "invalid use of slash\n");
         next_char(data);
         for (;;)
         {  if (data->c == '*')
            {  next_char(data);
               if (data->c == '/')
               {  next_char(data);
                  break;
               }
            }
            next_char(data);
         }
         goto loop;
      }
      return;
}

static void next_item(glp_data *data)
{     /* read next item */
      int len;
      skip_pad(data);
      len = 0;
      while (!(data->c == ' ' || data->c == '\n'))
      {  data->item[len++] = (char)data->c;
         if (len == sizeof(data->item))
            glp_sdf_error(data, "data item `%.31s...' too long\n",
               data->item);
         next_char(data);
      }
      data->item[len] = '\0';
      return;
}

int glp_sdf_read_int(glp_data *data)
{     /* read integer number */
      int x;
      next_item(data);
      switch (str2int(data->item, &x))
      {  case 0:
            break;
         case 1:
            glp_sdf_error(data, "integer `%s' out of range\n",
               data->item);
         case 2:
            glp_sdf_error(data, "cannot convert `%s' to integer\n",
               data->item);
         default:
            xassert(data != data);
      }
      return x;
}

double glp_sdf_read_num(glp_data *data)
{     /* read floating-point number */
      double x;
      next_item(data);
      switch (str2num(data->item, &x))
      {  case 0:
            break;
         case 1:
            glp_sdf_error(data, "number `%s' out of range\n",
               data->item);
         case 2:
            glp_sdf_error(data, "cannot convert `%s' to number\n",
               data->item);
         default:
            xassert(data != data);
      }
      return x;
}

const char *glp_sdf_read_item(glp_data *data)
{     /* read data item */
      next_item(data);
      return data->item;
}

const char *glp_sdf_read_text(glp_data *data)
{     /* read text until end of line */
      int c, len = 0;
      for (;;)
      {  c = data->c;
         next_char(data);
         if (c == ' ')
         {  /* ignore initial spaces */
            if (len == 0) continue;
            /* and multiple ones */
            if (data->item[len-1] == ' ') continue;
         }
         else if (c == '\n')
         {  /* remove trailing space */
            if (len > 0 && data->item[len-1] == ' ') len--;
            /* and stop reading */
            break;
         }
         /* add current character to the buffer */
         data->item[len++] = (char)c;
         if (len == sizeof(data->item))
            glp_sdf_error(data, "line too long\n", data->item);
      }
      data->item[len] = '\0';
      return data->item;
}

int glp_sdf_line(glp_data *data)
{     /* determine current line number */
      return data->count;
}

void glp_sdf_close_file(glp_data *data)
{     /* close plain data file */
      xfclose(data->fp);
      xfree(data->fname);
      xfree(data);
      return;
}

/* eof */
