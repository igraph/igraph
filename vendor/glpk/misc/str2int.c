/* str2int.c (convert string to int) */

/***********************************************************************
*  This code is part of GLPK (GNU Linear Programming Kit).
*  Copyright (C) 2000 Free Software Foundation, Inc.
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

#include "misc.h"
#include "stdc.h"

/***********************************************************************
*  NAME
*
*  str2int - convert character string to value of int type
*
*  SYNOPSIS
*
*  #include "misc.h"
*  int str2int(const char *str, int *val);
*
*  DESCRIPTION
*
*  The routine str2int converts the character string str to a value of
*  integer type and stores the value into location, which the parameter
*  val points to (in the case of error content of this location is not
*  changed).
*
*  RETURNS
*
*  The routine returns one of the following error codes:
*
*  0 - no error;
*  1 - value out of range;
*  2 - character string is syntactically incorrect. */

int str2int(const char *str, int *val_)
{     int d, k, s, val = 0;
      /* scan optional sign */
      if (str[0] == '+')
         s = +1, k = 1;
      else if (str[0] == '-')
         s = -1, k = 1;
      else
         s = +1, k = 0;
      /* check for the first digit */
      if (!isdigit((unsigned char)str[k]))
         return 2;
      /* scan digits */
      while (isdigit((unsigned char)str[k]))
      {  d = str[k++] - '0';
         if (s > 0)
         {  if (val > INT_MAX / 10)
               return 1;
            val *= 10;
            if (val > INT_MAX - d)
               return 1;
            val += d;
         }
         else /* s < 0 */
         {  if (val < INT_MIN / 10)
               return 1;
            val *= 10;
            if (val < INT_MIN + d)
               return 1;
            val -= d;
         }
      }
      /* check for terminator */
      if (str[k] != '\0')
         return 2;
      /* conversion has been done */
      *val_ = val;
      return 0;
}

/* eof */
