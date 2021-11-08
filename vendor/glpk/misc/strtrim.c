/* strtrim.c (remove trailing spaces from string) */

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
*  strtrim - remove trailing spaces from character string
*
*  SYNOPSIS
*
*  #include "misc.h"
*  char *strtrim(char *str);
*
*  DESCRIPTION
*
*  The routine strtrim removes trailing spaces from the character
*  string str.
*
*  RETURNS
*
*  The routine returns a pointer to the character string.
*
*  EXAMPLES
*
*  strtrim("Errare humanum est   ") => "Errare humanum est"
*
*  strtrim("      ")                => ""                             */

char *strtrim(char *str)
{     char *t;
      for (t = strrchr(str, '\0') - 1; t >= str; t--)
      {  if (*t != ' ')
            break;
         *t = '\0';
      }
      return str;
}

/* eof */
