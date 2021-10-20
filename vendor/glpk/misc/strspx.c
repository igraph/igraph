/* strspx.c (remove all spaces from string) */

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

/***********************************************************************
*  NAME
*
*  strspx - remove all spaces from character string
*
*  SYNOPSIS
*
*  #include "misc.h"
*  char *strspx(char *str);
*
*  DESCRIPTION
*
*  The routine strspx removes all spaces from the character string str.
*
*  RETURNS
*
*  The routine returns a pointer to the character string.
*
*  EXAMPLES
*
*  strspx("   Errare   humanum   est   ") => "Errarehumanumest"
*
*  strspx("      ")                       => ""                       */

char *strspx(char *str)
{     char *s, *t;
      for (s = t = str; *s; s++)
      {  if (*s != ' ')
            *t++ = *s;
      }
      *t = '\0';
      return str;
}

/* eof */
