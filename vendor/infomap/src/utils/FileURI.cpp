/*******************************************************************************
 Infomap software package for multi-level network clustering
 Copyright (c) 2013, 2014 Daniel Edler, Anton Holmgren, Martin Rosvall

 This file is part of the Infomap software package.
 See file LICENSE_GPLv3.txt for full license details.
 For more information, see <http://www.mapequation.org>
 ******************************************************************************/

#include "FileURI.h"
#include "convert.h"
#include <stdexcept>
#include <utility>

using std::string;

namespace infomap {

FileURI::FileURI(string filename, bool requireExtension)
    : m_filename(std::move(filename)), m_requireExtension(requireExtension)
{
  auto getErrorMessage = [](const auto& name, auto requireExt) {
    string s = io::Str() << "Filename '" << name << "' must match the pattern \"[dir/]name" << (requireExt ? ".extension\"" : "[.extension]\"");
    return s;
  };

  auto name = m_filename;
  auto pos = m_filename.find_last_of('/');
  if (pos != string::npos) {
    if (pos == m_filename.length()) // File could not end with slash
      throw std::invalid_argument(getErrorMessage(m_filename, m_requireExtension));
    m_directory = m_filename.substr(0, pos + 1); // Include the last slash in the directory
    name = m_filename.substr(pos + 1); // No slash in the name
  } else {
    m_directory = "";
  }

  pos = name.find_last_of('.');
  if (pos == string::npos || pos == 0 || pos == name.length() - 1) {
    if (pos != string::npos || m_requireExtension)
      throw std::invalid_argument(getErrorMessage(m_filename, m_requireExtension));
    m_name = name;
    m_extension = "";
  } else {
    m_name = name.substr(0, pos);
    m_extension = name.substr(pos + 1);
  }
}

} // namespace infomap
