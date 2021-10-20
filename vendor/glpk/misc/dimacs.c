/* dimacs.c (reading data in DIMACS format) */

/***********************************************************************
*  This code is part of GLPK (GNU Linear Programming Kit).
*  Copyright (C) 2009-2015 Free Software Foundation, Inc.
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

#include "dimacs.h"

void dmx_error(DMX *csa, const char *fmt, ...)
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

void dmx_warning(DMX *csa, const char *fmt, ...)
{     /* print warning message and continue processing */
      va_list arg;
      xprintf("%s:%d: warning: ", csa->fname, csa->count);
      va_start(arg, fmt);
      xvprintf(fmt, arg);
      va_end(arg);
      xprintf("\n");
      return;
}

void dmx_read_char(DMX *csa)
{     /* read character from input text file */
      int c;
      if (csa->c == '\n') csa->count++;
      c = glp_getc(csa->fp);
      if (c < 0)
      {  if (glp_ioerr(csa->fp))
            dmx_error(csa, "read error - %s", get_err_msg());
         else if (csa->c == '\n')
            dmx_error(csa, "unexpected end of file");
         else
         {  dmx_warning(csa, "missing final end of line");
            c = '\n';
         }
      }
      else if (c == '\n')
         ;
      else if (isspace(c))
         c = ' ';
      else if (iscntrl(c))
         dmx_error(csa, "invalid control character 0x%02X", c);
      csa->c = c;
      return;
}

void dmx_read_designator(DMX *csa)
{     /* read one-character line designator */
      xassert(csa->c == '\n');
      dmx_read_char(csa);
      for (;;)
      {  /* skip preceding white-space characters */
         while (csa->c == ' ')
            dmx_read_char(csa);
         if (csa->c == '\n')
         {  /* ignore empty line */
            if (!csa->empty)
            {  dmx_warning(csa, "empty line ignored");
               csa->empty = 1;
            }
            dmx_read_char(csa);
         }
         else if (csa->c == 'c')
         {  /* skip comment line */
            while (csa->c != '\n')
               dmx_read_char(csa);
            dmx_read_char(csa);
         }
         else
         {  /* hmm... looks like a line designator */
            csa->field[0] = (char)csa->c, csa->field[1] = '\0';
            /* check that it is followed by a white-space character */
            dmx_read_char(csa);
            if (!(csa->c == ' ' || csa->c == '\n'))
               dmx_error(csa, "line designator missing or invalid");
            break;
         }
      }
      return;
}

void dmx_read_field(DMX *csa)
{     /* read data field */
      int len = 0;
      /* skip preceding white-space characters */
      while (csa->c == ' ')
         dmx_read_char(csa);
      /* scan data field */
      if (csa->c == '\n')
         dmx_error(csa, "unexpected end of line");
      while (!(csa->c == ' ' || csa->c == '\n'))
      {  if (len == sizeof(csa->field)-1)
            dmx_error(csa, "data field '%.15s...' too long",
               csa->field);
         csa->field[len++] = (char)csa->c;
         dmx_read_char(csa);
      }
      csa->field[len] = '\0';
      return;
}

void dmx_end_of_line(DMX *csa)
{     /* skip white-space characters until end of line */
      while (csa->c == ' ')
         dmx_read_char(csa);
      if (csa->c != '\n')
         dmx_error(csa, "too many data fields specified");
      return;
}

void dmx_check_int(DMX *csa, double num)
{     /* print a warning if non-integer data are detected */
      if (!csa->nonint && num != floor(num))
      {  dmx_warning(csa, "non-integer data detected");
         csa->nonint = 1;
      }
      return;
}

/* eof */
