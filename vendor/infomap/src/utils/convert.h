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
#include <sstream>
#include <string>
#include <locale> // std::locale, std::tolower
#include <limits>
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
  inline std::string stringify(const T& x)
  {
    std::ostringstream o;
    if (!(o << x))
      throw std::runtime_error((o << "stringify(" << x << ")", o.str()));
    return o.str();
  }

  template <>
  inline std::string stringify(const bool& x)
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
      if (!item.empty()) {
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

  template <typename T>
  inline bool stringToValue(std::string const& str, T& value)
  {
    std::istringstream istream(str);
    return !!(istream >> value);
  }

  template <typename T>
  inline bool stringToUnsignedValue(std::string const& str, T& value)
  {
    std::istringstream istream(str);
    istream >> std::ws;
    if (!istream.good() || istream.peek() == '-' || istream.peek() == std::char_traits<char>::eof())
      return false;

    unsigned long long target = 0;
    if (!(istream >> target))
      return false;

    istream >> std::ws;
    if (!istream.eof() || target > std::numeric_limits<T>::max())
      return false;

    value = static_cast<T>(target);
    return true;
  }

  template <>
  inline bool stringToValue(std::string const& str, unsigned int& value)
  {
    return stringToUnsignedValue(str, value);
  }

  template <>
  inline bool stringToValue(std::string const& str, unsigned long& value)
  {
    return stringToUnsignedValue(str, value);
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

  // Whether `token` is exactly a non-negative decimal integer that fits an
  // unsigned int, storing it in `value` when it is.
  //
  // `std::istream >> unsigned` accepts a leading '-' and stores the value modulo
  // 2^32, so reading an id that way turned "-2" into 4294967294 instead of
  // failing: a phantom node in its own module, and the real node silently losing
  // its name. The link and *States rows never had that problem because they parse
  // ids by hand; this is the same rule, in a place every id reader can reach.
  //
  // An optional leading '+' is accepted because detail::parseUnsigned, which the link
  // rows use, accepts it -- the point of this function is to give every id reader the
  // same rule, so it has to match the reader that was already right. A '-' is refused,
  // and so is anything after the digits: the link rows refuse "1x" and "1.5" too, so
  // the *Vertices and .clu readers refusing them is the same convention rather than a
  // new one.
  //
  // The bound is checked after every digit, not once at the end, so `result` is at
  // most 4294967295 entering an iteration and the next step reaches at most
  // 42949672959 -- nowhere near the 64-bit ceiling. That is what keeps an arbitrarily
  // long digit string safe: it is rejected by the eleventh digit at the latest. Move
  // the check out of the loop and that stops being true.
  inline bool parseNonNegativeInteger(const std::string& token, unsigned int& value)
  {
    std::size_t index = 0;
    if (index < token.size() && token[index] == '+')
      ++index;
    if (index >= token.size())
      return false;

    unsigned long long result = 0;
    for (; index < token.size(); ++index) {
      const char c = token[index];
      if (c < '0' || c > '9')
        return false;
      result = result * 10 + static_cast<unsigned long long>(c - '0');
      if (result > std::numeric_limits<unsigned int>::max())
        return false;
    }

    value = static_cast<unsigned int>(result);
    return true;
  }

  // Defined in convert.cpp with fmt, to keep the heavy fmt header out of this
  // widely-included file. {:.{}g} reproduces std::setprecision(precision) with
  // the default floatfield, {:.{}f} reproduces std::fixed << setprecision(...);
  // both verified byte-for-byte against the previous iostream implementation.
  std::string toPrecision(double value, unsigned int precision = 10, bool fixed = false);

  //! Render a value as an RFC 4180 quoted CSV field: wrapped in quotes, with any quote in
  //! it written twice. Without the doubling a name containing a quote closes the field
  //! early and the row gains a column, which every CSV reader then mis-parses (#908).
  inline std::string csvQuoted(const std::string& value)
  {
    std::string quoted;
    quoted.reserve(value.size() + 2);
    quoted += '"';
    for (const char c : value) {
      if (c == '"')
        quoted += '"';
      quoted += c;
    }
    quoted += '"';
    return quoted;
  }

} // namespace io

} // namespace infomap

#endif // CONVERT_H_
