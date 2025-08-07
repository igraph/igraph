/*******************************************************************************
 Infomap software package for multi-level network clustering
 Copyright (c) 2013, 2014 Daniel Edler, Anton Holmgren, Martin Rosvall

 This file is part of the Infomap software package.
 See file LICENSE_GPLv3.txt for full license details.
 For more information, see <http://www.mapequation.org>
 ******************************************************************************/

#ifndef FILEURI_H_
#define FILEURI_H_

#include <iostream>
#include <string>

namespace infomap {

/**
 * Filename class to simplify handling of parts of a filename.
 * If a path is path/to/file.ext, the member methods of this class give these parts:
 * getFilename -> "path/to/file.ext"
 * getName -> "file"
 * getExtension -> "ext"
 * Can throw std::invalid_argument on creation.
 */
class FileURI {
public:
  FileURI() = default;

  explicit FileURI(std::string filename, bool requireExtension = false);

  const std::string& getFilename() const { return m_filename; }

  const std::string& getName() const { return m_name; }

  const std::string& getExtension() const { return m_extension; }

  friend std::ostream& operator<<(std::ostream& out, const FileURI& file)
  {
    return out << file.getFilename();
  }

private:
  std::string m_filename;
  bool m_requireExtension = false;

  std::string m_directory;
  std::string m_name;
  std::string m_extension;
};

} // namespace infomap

#endif // FILEURI_H_
