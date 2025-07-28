/*******************************************************************************
 Infomap software package for multi-level network clustering
 Copyright (c) 2013, 2014 Daniel Edler, Anton Holmgren, Martin Rosvall

 This file is part of the Infomap software package.
 See file LICENSE_GPLv3.txt for full license details.
 For more information, see <http://www.mapequation.org>
 ******************************************************************************/

#ifndef SAFEFILE_H_
#define SAFEFILE_H_

#include "../utils/convert.h"
#include <iostream>
#include <fstream>
#include <ios>
#include <cstdio>
#include <stdexcept>

namespace infomap {

using std::ifstream;
using std::ofstream;

/**
 * A wrapper for the C++ file stream class that automatically closes
 * the file stream when the destructor is called. Allocate it on the
 * stack to have it automatically closed when going out of scope.
 *
 * Note:
 * In C++, the only code that can be guaranteed to be executed after an
 * exception is thrown are the destructors of objects residing on the stack.
 *
 * You can exploit that fact to avoid resource leaks by tying all resources
 * to the lifespan of an object allocated on the stack. This technique is
 * called Resource Acquisition Is Initialization (RAII).
 *
 */
class SafeInFile : public ifstream {
public:
  SafeInFile(const std::string& filename, ios_base::openmode mode = ios_base::in)
      : ifstream(filename, mode)
  {
    if (fail())
      throw std::runtime_error(io::Str() << "Error opening file '" << filename << "'. Check that the path points to a file and that you have read permissions.");
  }

  ~SafeInFile() override
  {
    if (is_open())
      close();
  }
};

class SafeOutFile : public ofstream {
public:
  SafeOutFile(const std::string& filename, ios_base::openmode mode = ios_base::out)
      : ofstream(filename, mode)
  {
    if (fail())
      throw std::runtime_error(io::Str() << "Error opening file '" << filename << "'. Check that the directory you are writing to exists and that you have write permissions.");
  }

  ~SafeOutFile() override
  {
    if (is_open())
      close();
  }
};

inline bool isDirectoryWritable(const std::string& dir)
{
  std::string path = io::Str() << dir << "_1nf0m4p_.tmp";
  bool ok = true;
  try {
    SafeOutFile out(path);
  } catch (const std::runtime_error&) {
    ok = false;
  }
  if (ok)
    std::remove(path.c_str());
  return ok;
}

} // namespace infomap

#endif // SAFEFILE_H_
