/* str2num.c (convert string to double) */

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
*  str2num - convert character string to value of double type
*
*  SYNOPSIS
*
*  #include "misc.h"
*  int str2num(const char *str, double *val);
*
*  DESCRIPTION
*
*  The routine str2num converts the character string str to a value of
*  double type and stores the value into location, which the parameter
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

int str2num(const char *str, double *val_)
{     int k;
      double val;
      /* scan optional sign */
      k = (str[0] == '+' || str[0] == '-' ? 1 : 0);
      /* check for decimal point */
      if (str[k] == '.')
      {  k++;
         /* a digit should follow it */
         if (!isdigit((unsigned char)str[k]))
            return 2;
         k++;
         goto frac;
      }
      /* integer part should start with a digit */
      if (!isdigit((unsigned char)str[k]))
         return 2;
      /* scan integer part */
      while (isdigit((unsigned char)str[k]))
         k++;
      /* check for decimal point */
      if (str[k] == '.') k++;
frac: /* scan optional fraction part */
      while (isdigit((unsigned char)str[k]))
         k++;
      /* check for decimal exponent */
      if (str[k] == 'E' || str[k] == 'e')
      {  k++;
         /* scan optional sign */
         if (str[k] == '+' || str[k] == '-')
            k++;
         /* a digit should follow E, E+ or E- */
         if (!isdigit((unsigned char)str[k]))
            return 2;
      }
      /* scan optional exponent part */
      while (isdigit((unsigned char)str[k]))
         k++;
      /* check for terminator */
      if (str[k] != '\0')
         return 2;
      /* perform conversion */
      {  char *endptr;
         val = strtod(str, &endptr);
         if (*endptr != '\0')
            return 2;
      }
      /* check for overflow */
      if (!(-DBL_MAX <= val && val <= +DBL_MAX))
         return 1;
      /* check for underflow */
      if (-DBL_MIN < val && val < +DBL_MIN)
         val = 0.0;
      /* conversion has been done */
      *val_ = val;
      return 0;
}

/* eof */
