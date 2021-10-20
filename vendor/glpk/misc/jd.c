/* jd.c (conversions between calendar date and Julian day number) */

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

#include <stddef.h>
#include "jd.h"

/***********************************************************************
*  NAME
*
*  jday - convert calendar date to Julian day number
*
*  SYNOPSIS
*
*  #include "jd.h"
*  int jday(int d, int m, int y);
*
*  DESCRIPTION
*
*  The routine jday converts a calendar date, Gregorian calendar, to
*  corresponding Julian day number j.
*
*  From the given day d, month m, and year y, the Julian day number j
*  is computed without using tables.
*
*  The routine is valid for 1 <= y <= 4000.
*
*  RETURNS
*
*  The routine jday returns the Julian day number, or negative value if
*  the specified date is incorrect.
*
*  REFERENCES
*
*  R. G. Tantzen, Algorithm 199: conversions between calendar date and
*  Julian day number, Communications of the ACM, vol. 6, no. 8, p. 444,
*  Aug. 1963. */

int jday(int d, int m, int y)
{     int c, ya, j, dd;
      if (!(1 <= d && d <= 31 &&
            1 <= m && m <= 12 &&
            1 <= y && y <= 4000))
         return -1;
      if (m >= 3)
         m -= 3;
      else
         m += 9, y--;
      c = y / 100;
      ya = y - 100 * c;
      j = (146097 * c) / 4 + (1461 * ya) / 4 + (153 * m + 2) / 5 + d +
         1721119;
      jdate(j, &dd, NULL, NULL);
      if (d != dd)
         return -1;
      return j;
}

/***********************************************************************
*  NAME
*
*  jdate - convert Julian day number to calendar date
*
*  SYNOPSIS
*
*  #include "jd.h"
*  int jdate(int j, int *d, int *m, int *y);
*
*  DESCRIPTION
*
*  The routine jdate converts a Julian day number j to corresponding
*  calendar date, Gregorian calendar.
*
*  The day d, month m, and year y are computed without using tables and
*  stored in corresponding locations.
*
*  The routine is valid for 1721426 <= j <= 3182395.
*
*  RETURNS
*
*  If the conversion is successful, the routine returns zero, otherwise
*  non-zero.
*
*  REFERENCES
*
*  R. G. Tantzen, Algorithm 199: conversions between calendar date and
*  Julian day number, Communications of the ACM, vol. 6, no. 8, p. 444,
*  Aug. 1963. */

int jdate(int j, int *d_, int *m_, int *y_)
{     int d, m, y;
      if (!(1721426 <= j && j <= 3182395))
         return 1;
      j -= 1721119;
      y = (4 * j - 1) / 146097;
      j = (4 * j - 1) % 146097;
      d = j / 4;
      j = (4 * d + 3) / 1461;
      d = (4 * d + 3) % 1461;
      d = (d + 4) / 4;
      m = (5 * d - 3) / 153;
      d = (5 * d - 3) % 153;
      d = (d + 5) / 5;
      y = 100 * y + j;
      if (m <= 9)
         m += 3;
      else m -= 9,
         y++;
      if (d_ != NULL) *d_ = d;
      if (m_ != NULL) *m_ = m;
      if (y_ != NULL) *y_ = y;
      return 0;
}

#ifdef GLP_TEST
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

int main(void)
{     int jbeg, jend, j, d, m, y;
      jbeg = jday(1, 1, 1);
      jend = jday(31, 12, 4000);
      for (j = jbeg; j <= jend; j++)
      {  assert(jdate(j, &d, &m, &y) == 0);
         assert(jday(d, m, y) == j);
      }
      printf("Routines jday and jdate work correctly.\n");
      return 0;
}
#endif

/* eof */
