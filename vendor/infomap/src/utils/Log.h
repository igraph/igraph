/*******************************************************************************
 Infomap software package for multi-level network clustering
 Copyright (c) 2013, 2014 Daniel Edler, Anton Holmgren, Martin Rosvall

 This file is part of the Infomap software package.
 See file LICENSE_GPLv3.txt for full license details.
 For more information, see <http://www.mapequation.org>
 ******************************************************************************/

#ifndef LOG_H_
#define LOG_H_

#include "format_core.h"

#include <ostream>
#include <limits>
#include <iomanip>
#include <string>
#include <type_traits>

namespace infomap {

struct hideIf;

/// Receiver for complete log lines, installed with Log::setSink. Each line is
/// delivered newline-stripped, tagged with the verbosity level of the Log
/// object that wrote it (0 = default output, >= 1 = -v/-vv detail).
/// writeLine may be called from any thread that logs (at elevated verbosity,
/// OpenMP task threads emit detail lines); implementations synchronize
/// themselves. Visibility filtering (silent, verbosity, thread-mute) happens
/// upstream, so the sink only ever sees lines that would have been printed.
class LogSink {
public:
  virtual ~LogSink() = default;
  virtual void writeLine(unsigned int level, const std::string& line) = 0;
};

class Log {
  using ostreamFuncPtr = std::add_pointer_t<std::ostream&(std::ostream&)>;

public:
  class ScopedMute {
  public:
    explicit ScopedMute(bool enabled = true) : m_enabled(enabled)
    {
      if (m_enabled)
        ++s_threadMuteDepth;
    }

    ~ScopedMute()
    {
      if (m_enabled)
        --s_threadMuteDepth;
    }

    ScopedMute(const ScopedMute&) = delete;
    ScopedMute& operator=(const ScopedMute&) = delete;

  private:
    bool m_enabled;
  };

  /**
   * Log when level is below or equal Log::verboseLevel()
   * and maxLevel is above or equal Log::verboseLevel()
   */
  explicit Log(unsigned int level = 0, unsigned int maxLevel = std::numeric_limits<int>::max())
      : m_level(level), m_maxLevel(maxLevel), m_visible(isVisible(m_level, m_maxLevel)) {}

  bool isVisible() const { return isVisible(m_level, m_maxLevel); }

  void hide(bool value) { m_visible = !value && isVisible(); }

  // Defined out-of-line below, once hideIf is complete; a single non-ref-qualified
  // member covers both `log << hideIf(c)` and `Log() << hideIf(c)` (temporaries).
  Log& operator<<(const hideIf& manip);

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

  /// fmt-native logging: Log(level).print("{} of {}", a, b).
  /// The format string is checked at compile time under C++20 (consteval);
  /// pre-C++20 it is validated at runtime. Rendering (the heavy fmt/format.h)
  /// lives in vprint() in Log.cpp so this header only needs fmt's core API
  /// (fmt/base.h, via utils/format_core.h).
  template <typename... Args>
  Log& print(fmt::format_string<Args...> format, Args&&... args)
  {
    if (m_visible)
      vprint(format, fmt::make_format_args(args...));
    return *this;
  }

  static void init(unsigned int verboseLevel, bool silent, unsigned int numberPrecision)
  {
    setVerboseLevel(verboseLevel);
    setSilent(silent);
    precision(static_cast<std::streamsize>(numberPrecision));
  }

  static bool isVisible(unsigned int level, unsigned int maxLevel)
  {
    return s_threadMuteDepth == 0 && !s_silent && s_verboseLevel >= level && s_verboseLevel <= maxLevel;
  }

  static void setVerboseLevel(unsigned int level) { s_verboseLevel = level; }

  static void setSilent(bool silent) { s_silent = silent; }

  static bool isSilent() { return s_silent; }

  static bool isThreadMuted() { return s_threadMuteDepth > 0; }

  /// Set a custom output stream for all Log output.
  /// The caller must ensure @p os outlives all subsequent logging,
  /// or remains valid until another stream is configured.
  static void setOutputStream(std::ostream& os)
  {
    s_ostream = &os;
  }
  static void setOutputStream(std::ostream&&) = delete;

  /// Guarantee zero output: redirect to a null sink and set silent.
  static void setNoOutput();

  /// Install a line sink that replaces stream output entirely: every visible
  /// log line is delivered to @p sink (tagged with its message level) instead
  /// of the configured output stream. Pass nullptr to uninstall and restore
  /// stream output. The caller owns the sink and must keep it alive until it
  /// is uninstalled; install/uninstall must not race with a running engine
  /// (Log state is process-global, like init()/setOutputStream()).
  static void setSink(LogSink* sink);

  static LogSink* sink() { return s_sink; }

  /// Flush the calling thread's partially assembled sink lines (text written
  /// without a trailing newline) through to the sink. Call after an engine
  /// run returns so a trailing partial line is not silently dropped.
  static void flushSinkLines();

  static std::ostream& getOutputStream() { return ostream(); }

  /// Whether the active Log sink is the process stdout, i.e. not redirected via
  /// setOutputStream. Defined in Log.cpp so the stdout-identity comparison stays
  /// out of shared headers (C++ stream policy). Used by the console layer to
  /// avoid emitting ANSI into a redirected sink.
  static bool isWritingToStdout();

  static std::streamsize precision()
  {
    return s_sink ? sinkPrecision() : ostream().precision();
  }

  static std::streamsize precision(std::streamsize precision)
  {
    // When muted (a worker thread under --parallel-trials), deliberately do NOT touch the
    // shared global stream — concurrent precision() calls across workers would race on it.
    // The returned value is therefore a no-op echo, not the stream's previous precision.
    // This is safe because muted callers' save/restore pairs are balanced and never reach
    // the stream; do not "fix" this by querying ostream() here (that reintroduces the race).
    if (s_threadMuteDepth > 0)
      return precision;
    if (s_sink)
      return sinkPrecision(precision);
    return ostream().precision(precision);
  }

private:
  // Type-erased renderer for print(); defined in Log.cpp so the heavy
  // fmt/format.h (vformat) stays out of this widely-included header.
  void vprint(fmt::string_view format, fmt::format_args args);

  unsigned int m_level;
  unsigned int m_maxLevel;
  bool m_visible;
  // With a sink installed, bind to the calling thread's per-level line
  // assembler so every line is tagged with the level it was written at and
  // concurrent threads cannot tear each other's lines. m_level is declared
  // before m_ostream, so it is initialized first.
  std::ostream& m_ostream = streamFor(m_level);

  static std::ostream& ostream()
  {
#ifdef INFOMAP_R
    return *s_ostream;
#else
    return s_ostream ? *s_ostream : defaultStream();
#endif
  }
  static std::ostream& streamFor(unsigned int level)
  {
    return s_sink ? sinkStream(level) : ostream();
  }
  // Thread-local per-level line assembler; defined in Log.cpp.
  static std::ostream& sinkStream(unsigned int level);
  // Set/get the line-assembler precision for the calling thread (and the
  // default applied to assemblers created later); defined in Log.cpp.
  static std::streamsize sinkPrecision(std::streamsize precision);
  static std::streamsize sinkPrecision();
#ifndef INFOMAP_R
  static std::ostream& defaultStream();
#endif
  static std::ostream* s_ostream;
  static LogSink* s_sink;
  static unsigned int s_verboseLevel;
  static bool s_silent;
  static thread_local unsigned int s_threadMuteDepth;
};

struct hideIf {
  explicit hideIf(bool value) : hide(value) {}

  bool hide;
};

inline Log& Log::operator<<(const hideIf& manip)
{
  hide(manip.hide);
  return *this;
}

} // namespace infomap

#endif // LOG_H_
