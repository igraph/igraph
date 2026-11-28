/*******************************************************************************
 Infomap software package for multi-level network clustering
 Copyright (c) 2013, 2014 Daniel Edler, Anton Holmgren, Martin Rosvall

 This file is part of the Infomap software package.
 See file LICENSE_GPLv3.txt for full license details.
 For more information, see <http://www.mapequation.org>
 ******************************************************************************/

#include "Log.h"
#include "format.h"
#ifndef INFOMAP_R
#include <iostream>
#endif
#include <array>
#include <atomic>
#include <memory>
#include <streambuf>

namespace infomap {

namespace {

  class NullStreambuf : public std::streambuf {
  protected:
    int_type overflow(int_type c) override { return traits_type::not_eof(c); }
  };

  // Buffers characters until '\n' and delivers each complete line to the
  // installed sink, tagged with this buffer's level. One instance per
  // (thread, level) — see sinkStreamsForThisThread() — so concurrent
  // threads assemble their lines independently and never tear each other's
  // output. The sink pointer is re-read at delivery time, so an uninstall
  // between writes simply drops the remainder instead of dereferencing a
  // dead sink.
  class SinkLineBuf : public std::streambuf {
  public:
    explicit SinkLineBuf(unsigned int level) : m_level(level) {}

    void flushPending()
    {
      if (!m_line.empty())
        deliver();
    }

  protected:
    int_type overflow(int_type c) override
    {
      if (traits_type::eq_int_type(c, traits_type::eof()))
        return traits_type::not_eof(c);
      const char ch = traits_type::to_char_type(c);
      if (ch == '\n')
        deliver();
      else
        m_line.push_back(ch);
      return c;
    }

    std::streamsize xsputn(const char* s, std::streamsize n) override
    {
      for (std::streamsize i = 0; i < n; ++i)
        overflow(traits_type::to_int_type(s[i]));
      return n;
    }

  private:
    void deliver()
    {
      if (LogSink* sink = Log::sink())
        sink->writeLine(m_level, m_line);
      m_line.clear();
    }

    unsigned int m_level;
    std::string m_line;
  };

  struct SinkStream {
    explicit SinkStream(unsigned int level, std::streamsize precision)
        : buf(level), stream(&buf)
    {
      stream.precision(precision);
    }

    SinkLineBuf buf;
    std::ostream stream;
  };

  // Detail lines above -vvv share the deepest assembler; in practice levels
  // are 0..3.
  constexpr unsigned int kNumSinkLevels = 4;

  // The default precision applied to assemblers created after Log::init set
  // it. Atomic: worker threads read it when lazily creating their
  // thread-local assemblers (first -vv detail write) while the calling
  // thread may be inside a Log::precision save/restore pair.
  std::atomic<std::streamsize> g_sinkDefaultPrecision { 6 };

  using SinkStreams = std::array<std::unique_ptr<SinkStream>, kNumSinkLevels>;

  SinkStreams& sinkStreamsForThisThread()
  {
    thread_local SinkStreams streams;
    return streams;
  }

  SinkStream& sinkStreamAt(unsigned int level)
  {
    const auto index = level < kNumSinkLevels ? level : kNumSinkLevels - 1;
    auto& slot = sinkStreamsForThisThread()[index];
    if (!slot)
      slot = std::make_unique<SinkStream>(index, g_sinkDefaultPrecision);
    return *slot;
  }

} // namespace

#ifndef INFOMAP_R
std::ostream& Log::defaultStream()
{
  static std::ostream& s = std::cout;
  return s;
}
#endif

std::ostream* Log::s_ostream = nullptr;
LogSink* Log::s_sink = nullptr;
unsigned int Log::s_verboseLevel = 0;
bool Log::s_silent = false;
thread_local unsigned int Log::s_threadMuteDepth = 0;

void Log::setSink(LogSink* sink)
{
  s_sink = sink;
}

void Log::flushSinkLines()
{
  for (auto& slot : sinkStreamsForThisThread()) {
    if (slot)
      slot->buf.flushPending();
  }
}

std::ostream& Log::sinkStream(unsigned int level)
{
  return sinkStreamAt(level).stream;
}

std::streamsize Log::sinkPrecision()
{
  return sinkStreamAt(0).stream.precision();
}

std::streamsize Log::sinkPrecision(std::streamsize precision)
{
  // Apply to every assembler of the calling thread (save/restore pairs like
  // `auto old = Log::precision(9)` must round-trip regardless of the level
  // written next), and remember it for assemblers other threads create later.
  g_sinkDefaultPrecision = precision;
  const auto previous = sinkStreamAt(0).stream.precision();
  for (unsigned int level = 0; level < kNumSinkLevels; ++level)
    sinkStreamAt(level).stream.precision(precision);
  return previous;
}

void Log::setNoOutput()
{
  static NullStreambuf nullBuf;
  static std::ostream nullStream(&nullBuf);
  s_ostream = &nullStream;
  s_silent = true;
}

bool Log::isWritingToStdout()
{
  // With a sink installed, lines go to the sink, never the process stdout —
  // this keeps the Console layer from emitting ANSI or live-progress
  // redraws into structured log records.
  if (s_sink != nullptr)
    return false;
#ifdef INFOMAP_R
  // The R build always routes Log through the Rprintf bridge (s_ostream), never
  // the process stdout, and does not even include <iostream>.
  return false;
#else
  return &ostream() == &std::cout;
#endif
}

void Log::vprint(fmt::string_view format, fmt::format_args args)
{
  // Called only from print() under m_visible; render here so fmt/format.h stays
  // out of Log.h. Routes through the same m_ostream sink as operator<<, so the
  // R/Python redirect and verbosity/silent gating all still apply.
  m_ostream << fmt::vformat(format, args);
}

} // namespace infomap
