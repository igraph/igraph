/*******************************************************************************
 Infomap software package for multi-level network clustering
 Copyright (c) 2013, 2014 Daniel Edler, Anton Holmgren, Martin Rosvall

 This file is part of the Infomap software package.
 See file LICENSE_GPLv3.txt for full license details.
 For more information, see <http://www.mapequation.org>
 ******************************************************************************/

#ifndef INFOMAPERROR_H_
#define INFOMAPERROR_H_

#include <cstdint>
#include <stdexcept>
#include <string>

namespace infomap {

enum class ExitCode : std::uint8_t {
  Success = 0,
  InvalidArguments = 1,
  InputError = 2,
  OutputError = 3,
  // 4 (MemoryEstimateExceeded) and 5 (NoValidTrials) are reserved but not yet
  // wired to a code path. 6 is reserved for checkpoint support. Keep these
  // values stable when they are implemented.
  // MemoryEstimateExceeded = 4,
  // NoValidTrials = 5,
  InternalError = 7,
  // Ctrl-C / cooperative cancellation; 128 + SIGINT(2), the shell convention (#412).
  Interrupted = 130,
};

class InfomapError : public std::runtime_error {
public:
  InfomapError(ExitCode code, const std::string& message)
      : std::runtime_error(message),
        m_code(code) {}

  ExitCode code() const { return m_code; }

private:
  ExitCode m_code;
};

inline int exitCodeValue(ExitCode code)
{
  return static_cast<int>(code);
}

} // namespace infomap

#endif // INFOMAPERROR_H_
