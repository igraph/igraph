/* dimacs.h (reading data in DIMACS format) */

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

#ifndef DIMACS_H
#define DIMACS_H

#include "env.h"

typedef struct DMX DMX;

struct DMX
{     /* DIMACS data reader */
      jmp_buf jump;
      /* label for go to in case of error */
      const char *fname;
      /* name of input text file */
      glp_file *fp;
      /* stream assigned to input text file */
      int count;
      /* line count */
      int c;
      /* current character */
      char field[255+1];
      /* data field */
      int empty;
      /* warning 'empty line ignored' was printed */
      int nonint;
      /* warning 'non-integer data detected' was printed */
};

#define dmx_error _glp_dmx_error
void dmx_error(DMX *csa, const char *fmt, ...);
/* print error message and terminate processing */

#define dmx_warning _glp_dmx_warning
void dmx_warning(DMX *csa, const char *fmt, ...);
/* print warning message and continue processing */

#define dmx_read_char _glp_dmx_read_char
void dmx_read_char(DMX *csa);
/* read character from input text file */

#define dmx_read_designator _glp_dmx_read_designator
void dmx_read_designator(DMX *csa);
/* read one-character line designator */

#define dmx_read_field _glp_dmx_read_field
void dmx_read_field(DMX *csa);
/* read data field */

#define dmx_end_of_line _glp_dmx_end_of_line
void dmx_end_of_line(DMX *csa);
/* skip white-space characters until end of line */

#define dmx_check_int _glp_dmx_check_int
void dmx_check_int(DMX *csa, double num);
/* print a warning if non-integer data are detected */

#endif

/* eof */
