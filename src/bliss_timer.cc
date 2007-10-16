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
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#include <unistd.h>
#include <sys/times.h>
#include "bliss_timer.hh"

namespace igraph {

static const double numTicksPerSec = (double)(sysconf(_SC_CLK_TCK));

Timer::Timer()
{
  start_time = 0.0;
  end_time = 0.0;
}

void Timer::start()
{
  struct tms clkticks;

  times(&clkticks);
  start_time =
    ((double) clkticks.tms_utime + (double) clkticks.tms_stime) /
    numTicksPerSec;
}

void Timer::stop()
{
  struct tms clkticks;

  times(&clkticks);
  end_time =
    ((double) clkticks.tms_utime + (double) clkticks.tms_stime) /
    numTicksPerSec;
}

double Timer::get_intermediate()
{
  struct tms clkticks;

  times(&clkticks);
  double intermediate = 
    ((double) clkticks.tms_utime + (double) clkticks.tms_stime) /
    numTicksPerSec;
  return intermediate - start_time;
}

double Timer::get_duration()
{
  return(end_time - start_time);
}

}
