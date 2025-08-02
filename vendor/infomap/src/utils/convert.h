/*******************************************************************************
 Infomap software package for multi-level network clustering
 Copyright (c) 2013, 2014 Daniel Edler, Anton Holmgren, Martin Rosvall

 This file is part of the Infomap software package.
 See file LICENSE_GPLv3.txt for full license details.
 For more information, see <http://www.mapequation.org>
 ******************************************************************************/

#ifndef CONVERT_H_
#define CONVERT_H_

#include <stdexcept>
#include <iomanip>
#include <sstream>
#include <string>
#include <locale> // std::locale, std::tolower
#include <iostream>
#include <vector>

namespace infomap {

template <typename T>
struct TypeInfo {
  static bool isNumeric() { return false; }
};
template <>
struct TypeInfo<bool> {
  static bool isNumeric() { return false; }
};
template <>
struct TypeInfo<int> {
  static bool isNumeric() { return true; }
};
template <>
struct TypeInfo<unsigned int> {
  static bool isNumeric() { return true; }
};
template <>
struct TypeInfo<double> {
  static bool isNumeric() { return true; }
};

namespace io {

  inline std::string tolower(std::string str)
  {
    std::locale loc;
    for (char& c : str)
      c = std::tolower(c, loc);
    return str;
  }

  template <typename T>
  inline std::string stringify(T& x)
  {
    std::ostringstream o;
    if (!(o << x))
      throw std::runtime_error((o << "stringify(" << x << ")", o.str()));
    return o.str();
  }

  template <>
  inline std::string stringify(bool& x)
  {
    return x ? "true" : "false";
  }

  template <typename Container>
  inline std::string stringify(const Container& cont, const std::string& delimiter)
  {
    std::ostringstream o;
    if (cont.empty())
      return "";
    unsigned int maxIndex = cont.size() - 1;
    for (unsigned int i = 0; i < maxIndex; ++i) {
      if (!(o << cont[i]))
        throw std::runtime_error((o << "stringify(container[" << i << "])", o.str()));
      o << delimiter;
    }
    if (!(o << cont[maxIndex]))
      throw std::runtime_error((o << "stringify(container[" << maxIndex << "])", o.str()));
    return o.str();
  }

  struct InsensitiveCompare {
    bool operator()(const std::string& a, const std::string& b) const
    {
      auto lhs = a.begin();
      auto rhs = b.begin();

      std::locale loc;
      for (; lhs != a.end() && rhs != b.end(); ++lhs, ++rhs) {
        auto lhs_val = std::tolower(*lhs, loc);
        auto rhs_val = std::tolower(*rhs, loc);

        if (lhs_val != rhs_val)
          return lhs_val < rhs_val;
      }

      return (rhs != b.end());
    }
  };

  inline std::vector<std::string>& split(const std::string& s, char delim, std::vector<std::string>& items)
  {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
      if (item.length() > 0) {
        items.push_back(item);
      }
    }
    return items;
  }

  inline std::vector<std::string> split(const std::string& s, char delim)
  {
    std::vector<std::string> items;
    split(s, delim, items);
    return items;
  }

  class Str {
  public:
    Str() = default;
    template <class T>
    Str& operator<<(const T& t)
    {
      m_oss << stringify(t);
      return *this;
    }
    Str& operator<<(std::ostream& (*f)(std::ostream&))
    {
      m_oss << f;
      return *this;
    }
    operator std::string() const
    {
      return m_oss.str();
    }

  private:
    std::ostringstream m_oss;
  };

  template <typename T>
  inline bool stringToValue(std::string const& str, T& value)
  {
    std::istringstream istream(str);
    return !!(istream >> value);
  }

  template <>
  inline bool stringToValue(std::string const& str, unsigned int& value)
  {
    std::istringstream istream(str);
    int target = 0;
    istream >> target;
    if (target < 0) return false;
    value = target;
    return true;
  }

  template <>
  inline bool stringToValue(std::string const& str, unsigned long& value)
  {
    std::istringstream istream(str);
    int target = 0;
    istream >> target;
    if (target < 0) return false;
    value = target;
    return true;
  }

  inline std::string firstWord(const std::string& line)
  {
    std::istringstream ss;
    std::string buf;
    ss.str(line);
    ss >> buf;
    return buf;
  }

  template <typename T>
  inline std::string padValue(T value, const std::string::size_type size, bool rightAligned = true, const char paddingChar = ' ')
  {
    std::string valStr = stringify(value);
    if (size == valStr.size())
      return valStr;
    if (size < valStr.size())
      return valStr.substr(0, size);

    if (!rightAligned)
      return valStr.append(size - valStr.size(), paddingChar);

    return std::string(size - valStr.size(), paddingChar).append(valStr);
  }

  inline std::string toPrecision(double value, unsigned int precision = 10, bool fixed = false)
  {
    std::ostringstream o;
    if (fixed)
      o << std::fixed << std::setprecision(static_cast<int>(precision));
    else
      o << std::setprecision(static_cast<int>(precision));
    if (!(o << value))
      throw std::runtime_error((o << "stringify(" << value << ")", o.str()));
    return o.str();
  }

} // namespace io

} // namespace infomap

#endif // CONVERT_H_
