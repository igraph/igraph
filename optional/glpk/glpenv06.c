/* glpenv06.c (standard time) */

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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "glpapi.h"

/***********************************************************************
*  NAME
*
*  glp_time - determine current universal time
*
*  SYNOPSIS
*
*  glp_long glp_time(void);
*
*  RETURNS
*
*  The routine glp_time returns the current universal time (UTC), in
*  milliseconds, elapsed since 00:00:00 GMT January 1, 1970. */

static const int epoch = 2440588; /* = jday(1, 1, 1970) */

/* POSIX version ******************************************************/

#if defined(HAVE_SYS_TIME_H) && defined(HAVE_GETTIMEOFDAY)

#include <sys/time.h>
#include <time.h>

glp_long glp_time(void)
{     struct timeval tv;
      struct tm *tm;
      glp_long t;
      int j;
      gettimeofday(&tv, NULL);
      tm = gmtime(&tv.tv_sec);
      j = jday(tm->tm_mday, tm->tm_mon + 1, 1900 + tm->tm_year);
      xassert(j >= 0);
      t = xlset(j - epoch);
      t = xlmul(t, xlset(24));
      t = xladd(t, xlset(tm->tm_hour));
      t = xlmul(t, xlset(60));
      t = xladd(t, xlset(tm->tm_min));
      t = xlmul(t, xlset(60));
      t = xladd(t, xlset(tm->tm_sec));
      t = xlmul(t, xlset(1000));
      t = xladd(t, xlset(tv.tv_usec / 1000));
      return t;
}

/* Windows version ****************************************************/

#elif defined(__WOE__)

#include <windows.h>

glp_long glp_time(void)
{     SYSTEMTIME st;
      glp_long t;
      int j;
      GetSystemTime(&st);
      j = jday(st.wDay, st.wMonth, st.wYear);
      xassert(j >= 0);
      t = xlset(j - epoch);
      t = xlmul(t, xlset(24));
      t = xladd(t, xlset(st.wHour));
      t = xlmul(t, xlset(60));
      t = xladd(t, xlset(st.wMinute));
      t = xlmul(t, xlset(60));
      t = xladd(t, xlset(st.wSecond));
      t = xlmul(t, xlset(1000));
      t = xladd(t, xlset(st.wMilliseconds));
      return t;
}

/* portable ISO C version *********************************************/

#else

#include <time.h>

glp_long glp_time(void)
{     time_t timer;
      struct tm *tm;
      glp_long t;
      int j;
      timer = time(NULL);
      tm = gmtime(&timer);
      j = jday(tm->tm_mday, tm->tm_mon + 1, 1900 + tm->tm_year);
      xassert(j >= 0);
      t = xlset(j - epoch);
      t = xlmul(t, xlset(24));
      t = xladd(t, xlset(tm->tm_hour));
      t = xlmul(t, xlset(60));
      t = xladd(t, xlset(tm->tm_min));
      t = xlmul(t, xlset(60));
      t = xladd(t, xlset(tm->tm_sec));
      t = xlmul(t, xlset(1000));
      return t;
}

#endif

/***********************************************************************
*  NAME
*
*  glp_difftime - compute difference between two time values
*
*  SYNOPSIS
*
*  double glp_difftime(glp_long t1, glp_long t0);
*
*  RETURNS
*
*  The routine glp_difftime returns the difference between two time
*  values t1 and t0, expressed in seconds. */

double glp_difftime(glp_long t1, glp_long t0)
{     return
         xltod(xlsub(t1, t0)) / 1000.0;
}

/**********************************************************************/

#if 0
int main(void)
{     glp_long t;
      glp_ldiv d;
      int ttt, ss, mm, hh, day, month, year;
      char s[50];
      t = glp_time();
      xprintf("t = %s\n", xltoa(t, s));
      d = xldiv(t, xlset(1000));
      ttt = d.rem.lo, t = d.quot;
      d = xldiv(t, xlset(60));
      ss = d.rem.lo, t = d.quot;
      d = xldiv(t, xlset(60));
      mm = d.rem.lo, t = d.quot;
      d = xldiv(t, xlset(24));
      hh = d.rem.lo, t = d.quot;
      xassert(jdate(t.lo + epoch, &day, &month, &year) == 0);
      xprintf("%04d-%02d-%02d %02d:%02d:%02d.%03d\n", year, month, day,
         hh, mm, ss, ttt);
      return 0;
}
#endif

/* eof */
