/*******************************************************************************
 Infomap software package for multi-level network clustering
 Copyright (c) 2013, 2014 Daniel Edler, Anton Holmgren, Martin Rosvall

 This file is part of the Infomap software package.
 See file LICENSE_GPLv3.txt for full license details.
 For more information, see <http://www.mapequation.org>
 ******************************************************************************/

#include "Console.h"
#include "Log.h"
#include "convert.h"
#include "format.h"

#include <cmath>
#include <cstdlib>

#ifndef _WIN32
#include <unistd.h>
#endif

namespace infomap {

namespace {

  bool consoleSupportsAnsi()
  {
    // isatty(STDOUT) and TERM are fixed for the process lifetime — probe once.
    static const bool terminalIsAnsiCapable = []() {
#if defined(_WIN32) || defined(__EMSCRIPTEN__)
      return false;
#else
      const char* term = std::getenv("TERM");
      return ::isatty(STDOUT_FILENO) && term != nullptr && std::string(term) != "dumb";
#endif
    }();

    if (!terminalIsAnsiCapable)
      return false;

    // NO_COLOR can be set/unset after the first probe (tests, embedders), so honor
    // it per call — a cheap getenv, and it keeps behavior order-independent.
    if (std::getenv("NO_COLOR") != nullptr)
      return false;

    // Console text is routed through infomap::Log, whose sink can be redirected
    // away from stdout (Log::setOutputStream — R's Rprintf bridge, the binding
    // streams, test capture buffers). Only emit ANSI when Log actually writes to
    // the process stdout; otherwise escapes would leak into the redirected sink
    // even though STDOUT itself is a tty. The stdout-identity comparison lives in
    // Log.cpp so this file stays free of direct std-stream usage (stream policy).
    return Log::isWritingToStdout();
  }

} // namespace

Console::Console() : m_ansi(consoleSupportsAnsi()) {}

void Console::section(const std::string& title) const
{
  Log() << fmt::format(FMT_STRING("\n{}{}{}{}\n"), cyan(), bold(), title, reset());
}

void Console::metric(const std::string& label, const std::string& value) const
{
  // {:<30} == std::left << std::setw(30) (min width, left-justified, space fill).
  // The label (and its padding) is dimmed; reset() lands before the value.
  Log() << fmt::format(FMT_STRING("  {}{:<30}{}{}\n"), dim(), label, reset(), value);
}

std::string Console::statusLine(const std::string& label, const std::string& value) const
{
  return fmt::format(FMT_STRING("  {} {:<14}{}"), bullet(), label, value);
}

void Console::status(const std::string& label, const std::string& value) const
{
  Log() << statusLine(label, value) << "\n";
}

Console::Progress Console::progress(const std::string& label) const
{
  return Progress(label);
}

Console::Progress::Progress(std::string label) : m_label(std::move(label)) {}

Console::Progress::Progress(Progress&& other) noexcept
    : m_console(other.m_console), m_label(std::move(other.m_label)), m_live(other.m_live)
{
  other.m_live = false;
}

Console::Progress::~Progress()
{
  // If updates were drawn but finish() was never reached (e.g. an early return
  // or exception), commit a newline so following output doesn't glue onto the
  // dangling in-place line. Log(0, 0) matches update()'s gating, so this is a
  // no-op unless a live line was actually shown.
  if (m_live)
    Log(0, 0) << "\n";
}

void Console::Progress::update(const std::string& value)
{
  // In-place redraw: carriage return, the line body, then clear-to-end-of-line
  // (\033[K) to wipe any leftover from a previously longer value. Gated on a
  // TTY (escapes would corrupt a file) and on Log(0, 0), which is visible only
  // at verbosity 0 and is muted in silent / parallel-trial workers — so verbose
  // mode keeps its detailed per-step lines instead of a clobbered live line.
  if (!m_console.m_ansi)
    return;
  // Only claim the on-screen line if this verbosity actually draws it. Log(0, 0)
  // is muted in verbose / silent / parallel-trial modes, where no live line ever
  // appears — setting m_live there would lie about screen state and make finish()
  // emit a stray carriage-return/clear over output it never wrote.
  if (!Log(0, 0).isVisible())
    return;
  Log(0, 0) << "\r" << m_console.statusLine(m_label, value) << "\033[K" << std::flush;
  m_live = true;
}

void Console::Progress::finish(const std::string& value)
{
  // Freeze into the permanent status line, visible at verbosity 0 and above.
  // Overwrite the live line in place (\r ... \033[K) only if update() actually
  // drew one (m_live); otherwise — verbose mode, or finish() called without a
  // prior update() — emit the plain line, identical to status(), with no stray
  // carriage-return / clear sequence.
  if (m_console.m_ansi && m_live)
    Log() << "\r" << m_console.statusLine(m_label, value) << "\033[K\n";
  else
    Log() << m_console.statusLine(m_label, value) << "\n";
  m_live = false;
}

std::string Console::styleNote(const std::string& message)
{
  // Keeps the recognizable "  -> " marker (greppable; byte-identical to the
  // previous verbose notes in plain mode) but dims the arrow in a terminal so
  // verbose diagnostics share the level-0 pretty vocabulary.
  Console c;
  return fmt::format(FMT_STRING("  {}->{} {}\n"), c.dim(), c.reset(), message);
}

std::string Console::styleWarn(const std::string& message)
{
  // Warnings get their own marker (yellow "Warning:") instead of the note
  // arrow, so they stand apart from informational notes.
  Console c;
  return fmt::format(FMT_STRING("  {}Warning:{} {}\n"), c.yellow(), c.reset(), message);
}

std::string Console::styleDetail(const std::string& message)
{
  // Indent 4 (one step below the indent-2 status/metric tier) and dim, so the
  // bright section headers and • bullets stay the visual spine and trace reads
  // as subordinate. Plain mode keeps the indent, drops the dim.
  Console c;
  return fmt::format(FMT_STRING("    {}{}{}\n"), c.dim(), message, c.reset());
}

void Console::vNote(unsigned int level, fmt::string_view format, fmt::format_args args)
{
  Log(level) << styleNote(fmt::vformat(format, args));
}

void Console::vWarn(unsigned int level, fmt::string_view format, fmt::format_args args)
{
  Log(level) << styleWarn(fmt::vformat(format, args));
}

void Console::vDetail(unsigned int level, fmt::string_view format, fmt::format_args args)
{
  Log(level) << styleDetail(fmt::vformat(format, args));
}

std::string Console::arrow()
{
  // Same TTY/ANSI gate as bullet()/branch() (probed once, cached), so the
  // optimization status lines degrade to "->" when piped or on a dumb terminal.
  return consoleSupportsAnsi() ? "→" : "->";
}

std::string Console::percent(double value)
{
  if (std::abs(value) < 1e-9)
    return "0%";
  return fmt::format(FMT_STRING("{}%"), io::toPrecision(value, 2, true));
}

std::string Console::fixed(double value, unsigned int precision)
{
  // {:.{}f} reproduces std::fixed << std::setprecision(precision) byte-for-byte
  // (verified across edge cases); see scripts/check_cpp_stream_policy.py for why
  // formatting goes through fmt rather than iostreams.
  return fmt::format(FMT_STRING("{:.{}f}"), value, precision);
}

} // namespace infomap
