/*******************************************************************************
 Infomap software package for multi-level network clustering
 Copyright (c) 2013, 2014 Daniel Edler, Anton Holmgren, Martin Rosvall

 This file is part of the Infomap software package.
 See file LICENSE_GPLv3.txt for full license details.
 For more information, see <http://www.mapequation.org>
 ******************************************************************************/

#ifndef LOG_H_
#define LOG_H_

#include <iostream>
#include <limits>
#include <iomanip>
#include <type_traits>

namespace infomap {

struct hideIf;

class Log {
  using ostreamFuncPtr = std::add_pointer_t<std::ostream&(std::ostream&)>;

public:
  /**
   * Log when level is below or equal Log::verboseLevel()
   * and maxLevel is above or equal Log::verboseLevel()
   */
  explicit Log(unsigned int level = 0, unsigned int maxLevel = std::numeric_limits<int>::max())
      : m_level(level), m_maxLevel(maxLevel), m_visible(isVisible(m_level, m_maxLevel)) { }

  bool isVisible() const { return isVisible(m_level, m_maxLevel); }

  void hide(bool value) { m_visible = !value && isVisible(); }

  Log& operator<<(const hideIf&) { return *this; }

  template <typename T>
  Log& operator<<(const T& data)
  {
    if (m_visible)
      m_ostream << data;
    return *this;
  }

  Log& operator<<(ostreamFuncPtr f)
  {
    if (m_visible)
      m_ostream << f;
    return *this;
  }

  static void init(unsigned int verboseLevel, bool silent, unsigned int numberPrecision)
  {
    setVerboseLevel(verboseLevel);
    setSilent(silent);
    Log() << std::setprecision(static_cast<int>(numberPrecision));
  }

  static bool isVisible(unsigned int level, unsigned int maxLevel)
  {
    return !s_silent && s_verboseLevel >= level && s_verboseLevel <= maxLevel;
  }

  static void setVerboseLevel(unsigned int level) { s_verboseLevel = level; }

  static void setSilent(bool silent) { s_silent = silent; }

  static bool isSilent() { return s_silent; }

  static std::streamsize precision() { return std::cout.precision(); }

  static std::streamsize precision(std::streamsize precision)
  {
    return std::cout.precision(precision);
  }

private:
  unsigned int m_level;
  unsigned int m_maxLevel;
  bool m_visible;
  std::ostream& m_ostream = std::cout;

  static unsigned int s_verboseLevel;
  static bool s_silent;
};

struct hideIf {
  explicit hideIf(bool value) : hide(value) { }

  friend Log& operator<<(Log& out, const hideIf& manip)
  {
    out.hide(manip.hide);
    return out;
  }

  bool hide;
};

} // namespace infomap

#endif // LOG_H_
