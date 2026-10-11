/*******************************************************************************
 Infomap software package for multi-level network clustering
 Copyright (c) 2013, 2014 Daniel Edler, Anton Holmgren, Martin Rosvall

 This file is part of the Infomap software package.
 See file LICENSE_GPLv3.txt for full license details.
 For more information, see <http://www.mapequation.org>
 ******************************************************************************/

#ifndef INFOMAP_CONSOLE_H_
#define INFOMAP_CONSOLE_H_

#include "Log.h" // verbosity gating + the Log sink; also pulls in fmt's core API

#include <limits>
#include <string>

namespace infomap {

// The console presentation layer: renders structured program state to the
// terminal as section headers, key/value metrics, status bullets, live
// progress, and verbose notes/warnings, with ANSI styling and TTY detection.
// This is the human-facing console output, distinct from the file writers in
// io/ (Output, OutputPlan, ...). Everything is routed through infomap::Log.
class Console {
public:
  Console();

  const char* reset() const { return m_ansi ? "\033[0m" : ""; }
  const char* bold() const { return m_ansi ? "\033[1m" : ""; }
  const char* dim() const { return m_ansi ? "\033[2m" : ""; }
  const char* cyan() const { return m_ansi ? "\033[36m" : ""; }
  const char* green() const { return m_ansi ? "\033[32m" : ""; }
  const char* yellow() const { return m_ansi ? "\033[33m" : ""; }

  const char* bullet() const { return m_ansi ? "•" : "-"; }
  const char* branch() const { return m_ansi ? "╰─" : "  "; }

  void section(const std::string& title) const;
  void metric(const std::string& label, const std::string& value) const;
  void status(const std::string& label, const std::string& value) const;

  class Progress;

  // A live, single-line progress indicator for a long phase. On a TTY at
  // verbosity 0 it redraws one "• Label  <live value>" line in place as the
  // phase advances, then freezes into the permanent status() line via finish().
  // When output is piped/redirected, or in verbose/silent mode, the live
  // redraws are skipped and only the final frozen line is emitted (see
  // Console::Progress). The label is conventionally the one passed to finish().
  Progress progress(const std::string& label) const;

  // Verbose (-v/-vv) diagnostic emitters. Each gates on Log visibility for the
  // given level BEFORE formatting, so an invisible verbose line (the common case
  // at verbosity 0) costs nothing — no fmt::format, no allocation. Format string
  // is checked at compile time (fmt::format_string), mirroring Log::print; the
  // actual rendering lives in the v* helpers in Console.cpp so this header needs
  // only fmt's core API, not the full renderer.
  //
  //   Console::detail(1, "generate sub network with {} nodes", numNodes);
  //
  // Styling matches the former string-returning helpers exactly (see styleNote
  // /styleWarn/styleDetail). In plain mode (NO_COLOR / non-tty) the escapes
  // collapse to "" so the output stays greppable.

  // Notes: greppable "  -> " marker, dimmed arrow on a terminal.
  template <typename... Args>
  static void note(unsigned int level, fmt::format_string<Args...> format, Args&&... args)
  {
    if (Log::isVisible(level, std::numeric_limits<int>::max()))
      vNote(level, format, fmt::make_format_args(args...));
  }

  // Warnings: yellow "Warning:" marker, set apart from informational notes.
  template <typename... Args>
  static void warn(unsigned int level, fmt::format_string<Args...> format, Args&&... args)
  {
    if (Log::isVisible(level, std::numeric_limits<int>::max()))
      vWarn(level, format, fmt::make_format_args(args...));
  }

  // Tier-3 verbose trace: an indented, dimmed line subordinate to the section
  // and status lines above it, so verbose output reads as a bright skimmable
  // spine (sections + • bullets) with quiet detail nested beneath.
  template <typename... Args>
  static void detail(unsigned int level, fmt::format_string<Args...> format, Args&&... args)
  {
    if (Log::isVisible(level, std::numeric_limits<int>::max()))
      vDetail(level, format, fmt::make_format_args(args...));
  }

  // detail with an explicit upper verbosity bound (e.g. Log(1, 3)): visible only
  // when level <= verbosity <= maxLevel.
  template <typename... Args>
  static void detail(unsigned int level, unsigned int maxLevel, fmt::format_string<Args...> format, Args&&... args)
  {
    if (Log::isVisible(level, maxLevel))
      vDetail(level, format, fmt::make_format_args(args...));
  }

  // The "→" status connector, gated like bullet()/branch(): a real arrow on a
  // color-capable TTY, plain "->" otherwise so piped/dumb-terminal output stays
  // ASCII. Static (no per-call-site object) and used inside the format strings of
  // the optimization status/detail lines.
  static std::string arrow();

  static std::string percent(double value);
  static std::string fixed(double value, unsigned int precision = 6);

private:
  // Renders the "  • Label  value" body shared by status() and Progress
  // (no leading carriage return, no trailing newline).
  std::string statusLine(const std::string& label, const std::string& value) const;

  // Styling primitives for the verbose emitters: take an already-rendered
  // message and return the complete styled line (trailing newline). Plain mode
  // collapses the escapes to "". Internal — call note()/warn()/detail() instead.
  static std::string styleNote(const std::string& message);
  static std::string styleWarn(const std::string& message);
  static std::string styleDetail(const std::string& message);

  // Type-erased renderers behind note()/warn()/detail(): vformat the args, apply
  // the matching style, and emit at the given level. Defined in Console.cpp so
  // the heavy fmt renderer (vformat) stays out of this header — only reached
  // after the public templates have confirmed the line is visible.
  static void vNote(unsigned int level, fmt::string_view format, fmt::format_args args);
  static void vWarn(unsigned int level, fmt::string_view format, fmt::format_args args);
  static void vDetail(unsigned int level, fmt::string_view format, fmt::format_args args);

  bool m_ansi;

  friend class Progress;
};

// Handle for an in-place progress line; see Console::progress(). Construct
// one per phase, call update() in the work loop, and finish() with the final
// summary. Owns the on-screen line for its scope.
class Console::Progress {
public:
  explicit Progress(std::string label);
  // Movable so Console::progress() can return one by value; the moved-from
  // instance relinquishes the on-screen line (m_live = false) so only the live
  // owner commits a newline on destruction. Copying would duplicate ownership.
  Progress(Progress&& other) noexcept;
  Progress(const Progress&) = delete;
  Progress& operator=(const Progress&) = delete;
  Progress& operator=(Progress&&) = delete;
  ~Progress();

  // Redraw the live line with the current metric (TTY + verbosity 0 only).
  void update(const std::string& value);
  // Freeze into the permanent status line (always emitted when visible).
  void finish(const std::string& value);

private:
  Console m_console;
  std::string m_label;
  bool m_live = false; // an uncommitted in-place line is currently on screen
};

} // namespace infomap

#endif // INFOMAP_CONSOLE_H_
