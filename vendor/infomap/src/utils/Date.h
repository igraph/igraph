/*******************************************************************************
 Infomap software package for multi-level network clustering
 Copyright (c) 2013, 2014 Daniel Edler, Anton Holmgren, Martin Rosvall

 This file is part of the Infomap software package.
 See file LICENSE_GPLv3.txt for full license details.
 For more information, see <http://www.mapequation.org>
 ******************************************************************************/

#ifndef DATE_H_
#define DATE_H_

#include <ctime>
#include <cmath>
#include <ostream>

namespace infomap {

class ElapsedTime {
public:
  ElapsedTime(double elapsedTime = 0.0) : m_elapsedTime(elapsedTime) { }

  double getSeconds() const { return m_elapsedTime; }

  friend std::ostream& operator<<(std::ostream& out, const ElapsedTime& elapsedTime)
  {
    auto temp = static_cast<unsigned int>(std::floor(elapsedTime.getSeconds()));
    if (temp > 60) {
      if (temp > 3600) {
        if (temp > 86400) {
          out << temp / 86400 << "d ";
          temp %= 86400;
        }
        out << temp / 3600 << "h ";
        temp %= 3600;
      }
      out << temp / 60 << "m ";
      temp %= 60;
      out << temp << "s";
    } else {
      out << temp << "s";
    }
    return out;
  }

private:
  double m_elapsedTime;
};

class Date {
public:
  friend std::ostream& operator<<(std::ostream& out, const Date& date)
  {
    struct std::tm t = *localtime(&date.m_timeOfCreation);
    return out << "" << (t.tm_year + 1900) << (t.tm_mon < 9 ? "-0" : "-") << (t.tm_mon + 1) << (t.tm_mday < 10 ? "-0" : "-") << t.tm_mday << (t.tm_hour < 10 ? " 0" : " ") << t.tm_hour << (t.tm_min < 10 ? ":0" : ":") << t.tm_min << (t.tm_sec < 10 ? ":0" : ":") << t.tm_sec << "";
  }

  ElapsedTime operator-(const Date& date) const
  {
    return { difftime(m_timeOfCreation, date.m_timeOfCreation) };
  }

private:
  std::time_t m_timeOfCreation = time(nullptr);
};

} // namespace infomap

#endif // DATE_H_
