/*
Copyright (C) 2003-2006 Tommi Junttila

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License version 2
 as published by the Free Software Foundation.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA.
*/

/* FSF address fixed in the above notice on 1 Oct 2009 by Tamas Nepusz */

#include <unistd.h>
#ifdef HAVE_TIMES_H
#include <sys/times.h>
#else
#include <ctime>
#endif
#include "bliss_timer.hh"

using namespace std;

namespace igraph {

#ifdef HAVE_TIMES_H
static const double numTicksPerSec = (double)(sysconf(_SC_CLK_TCK));
#else
static const double numTicksPerSec = CLOCKS_PER_SEC;
#endif

Timer::Timer()
{
  start_time = 0.0;
  end_time = 0.0;
}

void Timer::start()
{

#ifdef HAVE_TIMES_H
  struct tms clkticks;

  times(&clkticks);
  start_time =
    ((double) clkticks.tms_utime + (double) clkticks.tms_stime) /
    numTicksPerSec;

#else

  clock_t clkticks;
  
  clkticks=clock();
  start_time = (double)clkticks / numTicksPerSec;

#endif
}

void Timer::stop()
{
#ifdef HAVE_TIMES_H
  struct tms clkticks;

  times(&clkticks);
  end_time =
    ((double) clkticks.tms_utime + (double) clkticks.tms_stime) /
    numTicksPerSec;
#else

  clock_t clkticks;
  
  clkticks=clock();
  end_time = (double)clkticks / numTicksPerSec;

#endif
}

double Timer::get_intermediate()
{
#ifdef HAVE_TIMES_H
  struct tms clkticks;

  times(&clkticks);
  double intermediate = 
    ((double) clkticks.tms_utime + (double) clkticks.tms_stime) /
    numTicksPerSec;

#else

  clock_t clkticks;
  
  clkticks=clock();
  double intermediate = (double)clkticks / numTicksPerSec;

#endif

  return intermediate - start_time;
}

double Timer::get_duration()
{
  return(end_time - start_time);
}

}
