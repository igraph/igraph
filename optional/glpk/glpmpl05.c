/* glpmpl05.c */

/***********************************************************************
*  This code is part of GLPK (GNU Linear Programming Kit).
*
*  Authors: Andrew Makhorin <mao@gnu.org>
*           Heinrich Schuchardt <xypron.glpk@gmx.de>
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

#ifdef __clang__
#pragma clang diagnostic ignored "-Wsign-conversion"
#endif

#define _GLPSTD_STDIO
#define _GLPSTD_TIME
#include "glpmpl.h"

double fn_gmtime(MPL *mpl)
{     /* obtain the current calendar time (UTC) */
      time_t timer;
      struct tm *tm;
      int j;
      time(&timer);
      if (timer == (time_t)(-1))
err:     error(mpl, "gmtime(); unable to obtain current calendar time");
      tm = gmtime(&timer);
      if (tm == NULL) goto err;
      j = jday(tm->tm_mday, tm->tm_mon + 1, 1900 + tm->tm_year);
      if (j < 0) goto err;
      return (((double)(j - jday(1, 1, 1970)) * 24.0 +
         (double)tm->tm_hour) * 60.0 + (double)tm->tm_min) * 60.0 +
         (double)tm->tm_sec;
}

static char *week[] = { "Monday", "Tuesday", "Wednesday", "Thursday",
      "Friday", "Saturday", "Sunday" };

static char *moon[] = { "January", "February", "March", "April", "May",
      "June", "July", "August", "September", "October", "November",
      "December" };

static void error1(MPL *mpl, const char *str, const char *s,
      const char *fmt, const char *f, const char *msg)
{     xprintf("Input string passed to str2time:\n");
      xprintf("%s\n", str);
      xprintf("%*s\n", (s - str) + 1, "^");
      xprintf("Format string passed to str2time:\n");
      xprintf("%s\n", fmt);
      xprintf("%*s\n", (f - fmt) + 1, "^");
      error(mpl, "%s", msg);
      /* no return */
}

double fn_str2time(MPL *mpl, const char *str, const char *fmt)
{     /* convert character string to the calendar time */
      int j, year, month, day, hh, mm, ss, zone;
      const char *s, *f;
      year = month = day = hh = mm = ss = -1, zone = INT_MAX;
      s = str;
      for (f = fmt; *f != '\0'; f++)
      {  if (*f == '%')
         {  f++;
            if (*f == 'b' || *f == 'h')
            {  /* the abbreviated month name */
               int k;
               char *name;
               if (month >= 0)
                  error1(mpl, str, s, fmt, f, "month multiply specified"
                     );
               while (*s == ' ') s++;
               for (month = 1; month <= 12; month++)
               {  name = moon[month-1];
                  for (k = 0; k <= 2; k++)
                  {  if (toupper((unsigned char)s[k]) !=
                         toupper((unsigned char)name[k])) goto next;
                  }
                  s += 3;
                  for (k = 3; name[k] != '\0'; k++)
                  {  if (toupper((unsigned char)*s) !=
                         toupper((unsigned char)name[k])) break;
                     s++;
                  }
                  break;
next:             ;
               }
               if (month > 12)
                  error1(mpl, str, s, fmt, f, "abbreviated month name m"
                     "issing or invalid");
            }
            else if (*f == 'd')
            {  /* the day of the month as a decimal number (01..31) */
               if (day >= 0)
                  error1(mpl, str, s, fmt, f, "day multiply specified");
               while (*s == ' ') s++;
               if (!('0' <= *s && *s <= '9'))
                  error1(mpl, str, s, fmt, f, "day missing or invalid");
               day = (*s++) - '0';
               if ('0' <= *s && *s <= '9')
                  day = 10 * day + ((*s++) - '0');
               if (!(1 <= day && day <= 31))
                  error1(mpl, str, s, fmt, f, "day out of range");
            }
            else if (*f == 'H')
            {  /* the hour as a decimal number, using a 24-hour clock
                  (00..23) */
               if (hh >= 0)
                  error1(mpl, str, s, fmt, f, "hour multiply specified")
                     ;
               while (*s == ' ') s++;
               if (!('0' <= *s && *s <= '9'))
                  error1(mpl, str, s, fmt, f, "hour missing or invalid")
                     ;
               hh = (*s++) - '0';
               if ('0' <= *s && *s <= '9')
                  hh = 10 * hh + ((*s++) - '0');
               if (!(0 <= hh && hh <= 23))
                  error1(mpl, str, s, fmt, f, "hour out of range");
            }
            else if (*f == 'm')
            {  /* the month as a decimal number (01..12) */
               if (month >= 0)
                  error1(mpl, str, s, fmt, f, "month multiply specified"
                     );
               while (*s == ' ') s++;
               if (!('0' <= *s && *s <= '9'))
                  error1(mpl, str, s, fmt, f, "month missing or invalid"
                     );
               month = (*s++) - '0';
               if ('0' <= *s && *s <= '9')
                  month = 10 * month + ((*s++) - '0');
               if (!(1 <= month && month <= 12))
                  error1(mpl, str, s, fmt, f, "month out of range");
            }
            else if (*f == 'M')
            {  /* the minute as a decimal number (00..59) */
               if (mm >= 0)
                  error1(mpl, str, s, fmt, f, "minute multiply specifie"
                     "d");
               while (*s == ' ') s++;
               if (!('0' <= *s && *s <= '9'))
                  error1(mpl, str, s, fmt, f, "minute missing or invali"
                     "d");
               mm = (*s++) - '0';
               if ('0' <= *s && *s <= '9')
                  mm = 10 * mm + ((*s++) - '0');
               if (!(0 <= mm && mm <= 59))
                  error1(mpl, str, s, fmt, f, "minute out of range");
            }
            else if (*f == 'S')
            {  /* the second as a decimal number (00..60) */
               if (ss >= 0)
                  error1(mpl, str, s, fmt, f, "second multiply specifie"
                     "d");
               while (*s == ' ') s++;
               if (!('0' <= *s && *s <= '9'))
                  error1(mpl, str, s, fmt, f, "second missing or invali"
                     "d");
               ss = (*s++) - '0';
               if ('0' <= *s && *s <= '9')
                  ss = 10 * ss + ((*s++) - '0');
               if (!(0 <= ss && ss <= 60))
                  error1(mpl, str, s, fmt, f, "second out of range");
            }
            else if (*f == 'y')
            {  /* the year without a century as a decimal number
                  (00..99); the values 00 to 68 mean the years 2000 to
                  2068 while the values 69 to 99 mean the years 1969 to
                  1999 */
               if (year >= 0)
                  error1(mpl, str, s, fmt, f, "year multiply specified")
                     ;
               while (*s == ' ') s++;
               if (!('0' <= *s && *s <= '9'))
                  error1(mpl, str, s, fmt, f, "year missing or invalid")
                     ;
               year = (*s++) - '0';
               if ('0' <= *s && *s <= '9')
                  year = 10 * year + ((*s++) - '0');
               year += (year >= 69 ? 1900 : 2000);
            }
            else if (*f == 'Y')
            {  /* the year as a decimal number, using the Gregorian
                  calendar */
               if (year >= 0)
                  error1(mpl, str, s, fmt, f, "year multiply specified")
                     ;
               while (*s == ' ') s++;
               if (!('0' <= *s && *s <= '9'))
                  error1(mpl, str, s, fmt, f, "year missing or invalid")
                     ;
               year = 0;
               for (j = 1; j <= 4; j++)
               {  if (!('0' <= *s && *s <= '9')) break;
                  year = 10 * year + ((*s++) - '0');
               }
               if (!(1 <= year && year <= 4000))
                  error1(mpl, str, s, fmt, f, "year out of range");
            }
            else if (*f == 'z')
            {  /* time zone offset in the form zhhmm */
               int z, hh, mm;
               if (zone != INT_MAX)
                  error1(mpl, str, s, fmt, f, "time zone offset multipl"
                     "y specified");
               while (*s == ' ') s++;
               if (*s == 'Z')
               {  z = hh = mm = 0, s++;
                  goto skip;
               }
               if (*s == '+')
                  z = +1, s++;
               else if (*s == '-')
                  z = -1, s++;
               else
                  error1(mpl, str, s, fmt, f, "time zone offset sign mi"
                     "ssing");
               hh = 0;
               for (j = 1; j <= 2; j++)
               {  if (!('0' <= *s && *s <= '9'))
err1:                error1(mpl, str, s, fmt, f, "time zone offset valu"
                        "e incomplete or invalid");
                  hh = 10 * hh + ((*s++) - '0');
               }
               if (hh > 23)
err2:             error1(mpl, str, s, fmt, f, "time zone offset value o"
                     "ut of range");
               if (*s == ':')
               {  s++;
                  if (!('0' <= *s && *s <= '9')) goto err1;
               }
               mm = 0;
               if (!('0' <= *s && *s <= '9')) goto skip;
               for (j = 1; j <= 2; j++)
               {  if (!('0' <= *s && *s <= '9')) goto err1;
                  mm = 10 * mm + ((*s++) - '0');
               }
               if (mm > 59) goto err2;
skip:          zone = z * (60 * hh + mm);
            }
            else if (*f == '%')
            {  /* literal % character */
               goto test;
            }
            else
               error1(mpl, str, s, fmt, f, "invalid conversion specifie"
                  "r");
         }
         else if (*f == ' ')
            ;
         else
test:    {  /* check a matching character in the input string */
            if (*s != *f)
               error1(mpl, str, s, fmt, f, "character mismatch");
            s++;
         }
      }
      if (year < 0) year = 1970;
      if (month < 0) month = 1;
      if (day < 0) day = 1;
      if (hh < 0) hh = 0;
      if (mm < 0) mm = 0;
      if (ss < 0) ss = 0;
      if (zone == INT_MAX) zone = 0;
      j = jday(day, month, year);
      xassert(j >= 0);
      return (((double)(j - jday(1, 1, 1970)) * 24.0 + (double)hh) *
         60.0 + (double)mm) * 60.0 + (double)ss - 60.0 * (double)zone;
}

static void error2(MPL *mpl, const char *fmt, const char *f,
      const char *msg)
{     xprintf("Format string passed to time2str:\n");
      xprintf("%s\n", fmt);
      xprintf("%*s\n", (f - fmt) + 1, "^");
      error(mpl, "%s", msg);
      /* no return */
}

static int weekday(int j)
{     /* determine weekday number (1 = Mon, ..., 7 = Sun) */
      return (j + jday(1, 1, 1970)) % 7 + 1;
}

static int firstday(int year)
{     /* determine the first day of the first week for a specified year
         according to ISO 8601 */
      int j;
      /* if 1 January is Monday, Tuesday, Wednesday or Thursday, it is
         in week 01; if 1 January is Friday, Saturday or Sunday, it is
         in week 52 or 53 of the previous year */
      j = jday(1, 1, year) - jday(1, 1, 1970);
      switch (weekday(j))
      {  case 1: /* 1 Jan is Mon */ j += 0; break;
         case 2: /* 1 Jan is Tue */ j -= 1; break;
         case 3: /* 1 Jan is Wed */ j -= 2; break;
         case 4: /* 1 Jan is Thu */ j -= 3; break;
         case 5: /* 1 Jan is Fri */ j += 3; break;
         case 6: /* 1 Jan is Sat */ j += 2; break;
         case 7: /* 1 Jan is Sun */ j += 1; break;
         default: xassert(j != j);
      }
      /* the first day of the week must be Monday */
      xassert(weekday(j) == 1);
      return j;
}

void fn_time2str(MPL *mpl, char *str, double t, const char *fmt)
{     /* convert the calendar time to character string */
      int j, year, month, day, hh, mm, ss, len;
      double temp;
      const char *f;
      char buf[MAX_LENGTH+1];
      if (!(-62135596800.0 <= t && t <= 64092211199.0))
         error(mpl, "time2str(%.*g,...); argument out of range",
            DBL_DIG, t);
      t = floor(t + 0.5);
      temp = fabs(t) / 86400.0;
      j = (int)floor(temp);
      if (t < 0.0)
      {  if (temp == floor(temp))
            j = - j;
         else
            j = - (j + 1);
      }
      xassert(jdate(j + jday(1, 1, 1970), &day, &month, &year) == 0);
      ss = (int)(t - 86400.0 * (double)j);
      xassert(0 <= ss && ss < 86400);
      mm = ss / 60, ss %= 60;
      hh = mm / 60, mm %= 60;
      len = 0;
      for (f = fmt; *f != '\0'; f++)
      {  if (*f == '%')
         {  f++;
            if (*f == 'a')
            {  /* the abbreviated weekday name */
               memcpy(buf, week[weekday(j)-1], 3), buf[3] = '\0';
            }
            else if (*f == 'A')
            {  /* the full weekday name */
               strcpy(buf, week[weekday(j)-1]);
            }
            else if (*f == 'b' || *f == 'h')
            {  /* the abbreviated month name */
               memcpy(buf, moon[month-1], 3), buf[3] = '\0';
            }
            else if (*f == 'B')
            {  /* the full month name */
               strcpy(buf, moon[month-1]);
            }
            else if (*f == 'C')
            {  /* the century of the year */
               sprintf(buf, "%02d", year / 100);
            }
            else if (*f == 'd')
            {  /* the day of the month as a decimal number (01..31) */
               sprintf(buf, "%02d", day);
            }
            else if (*f == 'D')
            {  /* the date using the format %m/%d/%y */
               sprintf(buf, "%02d/%02d/%02d", month, day, year % 100);
            }
            else if (*f == 'e')
            {  /* the day of the month like with %d, but padded with
                  blank (1..31) */
               sprintf(buf, "%2d", day);
            }
            else if (*f == 'F')
            {  /* the date using the format %Y-%m-%d */
               sprintf(buf, "%04d-%02d-%02d", year, month, day);
            }
            else if (*f == 'g')
            {  /* the year corresponding to the ISO week number, but
                  without the century (range 00 through 99); this has
                  the same format and value as %y, except that if the
                  ISO week number (see %V) belongs to the previous or
                  next year, that year is used instead */
               int iso;
               if (j < firstday(year))
                  iso = year - 1;
               else if (j < firstday(year + 1))
                  iso = year;
               else
                  iso = year + 1;
               sprintf(buf, "%02d", iso % 100);
            }
            else if (*f == 'G')
            {  /* the year corresponding to the ISO week number; this
                  has the same format and value as %Y, excepth that if
                  the ISO week number (see %V) belongs to the previous
                  or next year, that year is used instead */
               int iso;
               if (j < firstday(year))
                  iso = year - 1;
               else if (j < firstday(year + 1))
                  iso = year;
               else
                  iso = year + 1;
               sprintf(buf, "%04d", iso);
            }
            else if (*f == 'H')
            {  /* the hour as a decimal number, using a 24-hour clock
                  (00..23) */
               sprintf(buf, "%02d", hh);
            }
            else if (*f == 'I')
            {  /* the hour as a decimal number, using a 12-hour clock
                  (01..12) */
               sprintf(buf, "%02d",
                  hh == 0 ? 12 : hh <= 12 ? hh : hh - 12);
            }
            else if (*f == 'j')
            {  /* the day of the year as a decimal number (001..366) */
               sprintf(buf, "%03d",
                  jday(day, month, year) - jday(1, 1, year) + 1);
            }
            else if (*f == 'k')
            {  /* the hour as a decimal number, using a 24-hour clock
                  like %H, but padded with blank (0..23) */
               sprintf(buf, "%2d", hh);
            }
            else if (*f == 'l')
            {  /* the hour as a decimal number, using a 12-hour clock
                  like %I, but padded with blank (1..12) */
               sprintf(buf, "%2d",
                  hh == 0 ? 12 : hh <= 12 ? hh : hh - 12);
            }
            else if (*f == 'm')
            {  /* the month as a decimal number (01..12) */
               sprintf(buf, "%02d", month);
            }
            else if (*f == 'M')
            {  /* the minute as a decimal number (00..59) */
               sprintf(buf, "%02d", mm);
            }
            else if (*f == 'p')
            {  /* either AM or PM, according to the given time value;
                  noon is treated as PM and midnight as AM */
               strcpy(buf, hh <= 11 ? "AM" : "PM");
            }
            else if (*f == 'P')
            {  /* either am or pm, according to the given time value;
                  noon is treated as pm and midnight as am */
               strcpy(buf, hh <= 11 ? "am" : "pm");
            }
            else if (*f == 'r')
            {  /* the calendar time using the format %I:%M:%S %p */
               sprintf(buf, "%02d:%02d:%02d %s",
                  hh == 0 ? 12 : hh <= 12 ? hh : hh - 12,
                  mm, ss, hh <= 11 ? "AM" : "PM");
            }
            else if (*f == 'R')
            {  /* the hour and minute using the format %H:%M */
               sprintf(buf, "%02d:%02d", hh, mm);
            }
            else if (*f == 'S')
            {  /* the second as a decimal number (00..59) */
               sprintf(buf, "%02d", ss);
            }
            else if (*f == 'T')
            {  /* the time of day using the format %H:%M:%S */
               sprintf(buf, "%02d:%02d:%02d", hh, mm, ss);
            }
            else if (*f == 'u')
            {  /* the day of the week as a decimal number (1..7),
                  Monday being 1 */
               sprintf(buf, "%d", weekday(j));
            }
            else if (*f == 'U')
            {  /* the week number of the current year as a decimal
                  number (range 00 through 53), starting with the first
                  Sunday as the first day of the first week; days
                  preceding the first Sunday in the year are considered
                  to be in week 00 */
#if 1 /* 09/I-2009 */
#undef sun
/* causes compilation error in SunOS */
#endif
               int sun;
               /* sun = the first Sunday of the year */
               sun = jday(1, 1, year) - jday(1, 1, 1970);
               sun += (7 - weekday(sun));
               sprintf(buf, "%02d", (j + 7 - sun) / 7);
            }
            else if (*f == 'V')
            {  /* the ISO week number as a decimal number (range 01
                  through 53); ISO weeks start with Monday and end with
                  Sunday; week 01 of a year is the first week which has
                  the majority of its days in that year; week 01 of
                  a year can contain days from the previous year; the
                  week before week 01 of a year is the last week (52 or
                  53) of the previous year even if it contains days
                  from the new year */
               int iso;
               if (j < firstday(year))
                  iso = j - firstday(year - 1);
               else if (j < firstday(year + 1))
                  iso = j - firstday(year);
               else
                  iso = j - firstday(year + 1);
               sprintf(buf, "%02d", iso / 7 + 1);
            }
            else if (*f == 'w')
            {  /* the day of the week as a decimal number (0..6),
                  Sunday being 0 */
               sprintf(buf, "%d", weekday(j) % 7);
            }
            else if (*f == 'W')
            {  /* the week number of the current year as a decimal
                  number (range 00 through 53), starting with the first
                  Monday as the first day of the first week; days
                  preceding the first Monday in the year are considered
                  to be in week 00 */
               int mon;
               /* mon = the first Monday of the year */
               mon = jday(1, 1, year) - jday(1, 1, 1970);
               mon += (8 - weekday(mon)) % 7;
               sprintf(buf, "%02d", (j + 7 - mon) / 7);
            }
            else if (*f == 'y')
            {  /* the year without a century as a decimal number
                  (00..99) */
               sprintf(buf, "%02d", year % 100);
            }
            else if (*f == 'Y')
            {  /* the year as a decimal number, using the Gregorian
                  calendar */
               sprintf(buf, "%04d", year);
            }
            else if (*f == '%')
            {  /* a literal % character */
               buf[0] = '%', buf[1] = '\0';
            }
            else
               error2(mpl, fmt, f, "invalid conversion specifier");
         }
         else
            buf[0] = *f, buf[1] = '\0';
         if (len + strlen(buf) > MAX_LENGTH)
            error(mpl, "time2str; output string length exceeds %d chara"
               "cters", MAX_LENGTH);
         memcpy(str+len, buf, strlen(buf));
         len += strlen(buf);
      }
      str[len] = '\0';
      return;
}

/* eof */
