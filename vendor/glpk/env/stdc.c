/* stdc.c (replacements for standard non-thread-safe functions) */

/***********************************************************************
*  This code is part of GLPK (GNU Linear Programming Kit).
*  Copyright (C) 2017 Free Software Foundation, Inc.
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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/* portable ANSI C version ********************************************/

#if !defined(TLS)

#define ENABLE_NON_SAFE
#include "stdc.h"

struct tm *xgmtime(const time_t *timer)
{     return
         gmtime(timer);
}

char *xstrerr(int errnum)
{     return
         strerror(errnum);
}

char *xstrtok(char *s1, const char *s2)
{     return
         strtok(s1, s2);
}

/* MS Windows version *************************************************/

#elif defined(__WOE__)

#include "stdc.h"

struct tm *xgmtime(const time_t *timer)
{     static TLS struct tm result;
      gmtime_s(&result, timer);
      return &result;
}

char *xstrerr(int errnum)
{     static TLS char s[1023+1];
      strerror_s(s, sizeof(s), errnum);
      return s;
}

char *xstrtok(char *s1, const char *s2)
{     static TLS char *ptr;
      return strtok_s(s1, s2, &ptr);
}

/* GNU/Linux version **************************************************/

#else

#include "stdc.h"

struct tm *xgmtime(const time_t *timer)
{     static TLS struct tm result;
      gmtime_r(timer, &result);
      return &result;
}

char *xstrerr(int errnum)
{     static TLS char s[1023+1];
      strerror_r(errnum, s, sizeof(s));
      return s;
}

char *xstrtok(char *s1, const char *s2)
{     static TLS char *ptr;
      return strtok_r(s1, s2, &ptr);
}

#endif

/* eof */
