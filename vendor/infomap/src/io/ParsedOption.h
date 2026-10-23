/*******************************************************************************
 Infomap software package for multi-level network clustering
 Copyright (c) 2013, 2014 Daniel Edler, Anton Holmgren, Martin Rosvall

 This file is part of the Infomap software package.
 See file LICENSE_GPLv3.txt for full license details.
 For more information, see <http://www.mapequation.org>
 ******************************************************************************/

#ifndef PARSED_OPTION_H_
#define PARSED_OPTION_H_

#include <ostream>
#include <string>

namespace infomap {

struct ParsedOption {
  char shortName = '\0';
  std::string longName;
  std::string description;
  std::string group;
  bool isAdvanced = false;
  bool requireArgument = false;
  bool incrementalArgument = false;
  std::string argumentName;
  bool negated = false;
  std::string value;
};

inline std::ostream& operator<<(std::ostream& out, const ParsedOption& option)
{
  if (option.negated)
    out << "no ";
  out << option.longName;
  if (option.requireArgument)
    out << " = " << option.value;
  return out;
}

} // namespace infomap

#endif // PARSED_OPTION_H_
