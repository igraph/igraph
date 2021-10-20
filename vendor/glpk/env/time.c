/* time.c (standard time) */

/***********************************************************************
*  This code is part of GLPK (GNU Linear Programming Kit).
*  Copyright (C) 2000-2017 Free Software Foundation, Inc.
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

#include "env.h"
#include "jd.h"

/***********************************************************************
*  NAME
*
*  glp_time - determine current universal time
*
*  SYNOPSIS
*
*  double glp_time(void);
*
*  RETURNS
*
*  The routine glp_time returns the current universal time (UTC), in
*  milliseconds, elapsed since 00:00:00 GMT January 1, 1970. */

#define EPOCH 2440588 /* = jday(1, 1, 1970) */

/* POSIX version ******************************************************/

#if defined(HAVE_SYS_TIME_H) && defined(HAVE_GETTIMEOFDAY)

#if 0 /* 29/VI-2017 */
#include <sys/time.h>
#include <time.h>

double glp_time(void)
{     struct timeval tv;
      struct tm *tm;
      int j;
      double t;
      gettimeofday(&tv, NULL);
#if 0 /* 29/I-2017 */
      tm = gmtime(&tv.tv_sec);
#else
      tm = xgmtime(&tv.tv_sec);
#endif
      j = jday(tm->tm_mday, tm->tm_mon + 1, 1900 + tm->tm_year);
      xassert(j >= 0);
      t = ((((double)(j - EPOCH) * 24.0 + (double)tm->tm_hour) * 60.0 +
         (double)tm->tm_min) * 60.0 + (double)tm->tm_sec) * 1000.0 +
         (double)(tv.tv_usec / 1000);
      return t;
}
#else
#include <sys/time.h>

double glp_time(void)
{     struct timeval tv;
      double t;
      gettimeofday(&tv, NULL);
      t = (double)tv.tv_sec + (double)(tv.tv_usec) / 1e6;
      xassert(0.0 <= t && t < 4294967296.0);
      return 1000.0 * t;
}
#endif

/* MS Windows version *************************************************/

#elif defined(__WOE__)

#include <windows.h>

double glp_time(void)
{     SYSTEMTIME st;
      int j;
      double t;
      GetSystemTime(&st);
      j = jday(st.wDay, st.wMonth, st.wYear);
      xassert(j >= 0);
      t = ((((double)(j - EPOCH) * 24.0 + (double)st.wHour) * 60.0 +
         (double)st.wMinute) * 60.0 + (double)st.wSecond) * 1000.0 +
         (double)st.wMilliseconds;
      return t;
}

/* portable ANSI C version ********************************************/

#else

#include <time.h>

double glp_time(void)
{     time_t timer;
      struct tm *tm;
      int j;
      double t;
      timer = time(NULL);
#if 0 /* 29/I-2017 */
      tm = gmtime(&timer);
#else
      tm = xgmtime(&timer);
#endif
      j = jday(tm->tm_mday, tm->tm_mon + 1, 1900 + tm->tm_year);
      xassert(j >= 0);
      t = ((((double)(j - EPOCH) * 24.0 + (double)tm->tm_hour) * 60.0 +
         (double)tm->tm_min) * 60.0 + (double)tm->tm_sec) * 1000.0;
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
*  double glp_difftime(double t1, double t0);
*
*  RETURNS
*
*  The routine glp_difftime returns the difference between two time
*  values t1 and t0, expressed in seconds. */

double glp_difftime(double t1, double t0)
{     return
         (t1 - t0) / 1000.0;
}

/* eof */
