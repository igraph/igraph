#[=======================================================================[.rst:
FindGLPK
--------

Finds the GLPK library.

Result Variables
^^^^^^^^^^^^^^^^

This will define the following variables:

``GLPK_FOUND``
  True if the system has the GLPK library.
``GLPK_VERSION``
  The version of the GLPK library which was found.
``GLPK_INCLUDE_DIRS``
  Include directories needed to use Foo.
``GLPK_LIBRARIES``
  Libraries needed to link to Foo.

Cache Variables
^^^^^^^^^^^^^^^

The following cache variables may also be set:

``GLPK_INCLUDE_DIR``
  The directory containing ``glpk.h``.
``GLPK_LIBRARY``
  The path to the GLPK library.

#]=======================================================================]

find_path(GLPK_INCLUDE_DIR
  NAMES glpk.h
)

find_library(GLPK_LIBRARY
  NAMES glpk
)

# parse version from header
if(GLPK_INCLUDE_DIR)
  set(GLPK_VERSION_FILE ${GLPK_INCLUDE_DIR}/glpk.h)
  file(READ ${GLPK_VERSION_FILE} GLPK_VERSION_FILE_CONTENTS)

  string(REGEX MATCH "#define[ ]+GLP_MAJOR_VERSION[ ]+[0-9]+"
    GLPK_VERSION_MAJOR "${GLPK_VERSION_FILE_CONTENTS}")
  string(REGEX REPLACE "#define[ ]+GLP_MAJOR_VERSION[ ]+([0-9]+)" "\\1"
    GLPK_VERSION_MAJOR "${GLPK_VERSION_MAJOR}")

  string(REGEX MATCH "#define[ ]+GLP_MINOR_VERSION[ ]+[0-9]+"
    GLPK_VERSION_MINOR "${GLPK_VERSION_FILE_CONTENTS}")
  string(REGEX REPLACE "#define[ ]+GLP_MINOR_VERSION[ ]+([0-9]+)" "\\1"
    GLPK_VERSION_MINOR "${GLPK_VERSION_MINOR}")

  set(GLPK_VERSION "${GLPK_VERSION_MAJOR}.${GLPK_VERSION_MINOR}")

  # compatibility variables
  set(GLPK_VERSION_STRING "${GLPK_VERSION}")
endif()

# behave like a CMake module is supposed to behave
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GLPK
  FOUND_VAR GLPK_FOUND
  REQUIRED_VARS
    GLPK_LIBRARY
    GLPK_INCLUDE_DIR
  VERSION_VAR GLPK_VERSION
)

# hide the introduced cmake cached variables in cmake GUIs
mark_as_advanced(
  GLPK_INCLUDE_DIR
  GLPK_LIBRARY
)

if(GLPK_FOUND)
  set(GLPK_LIBRARIES ${GLPK_LIBRARY})
  set(GLPK_INCLUDE_DIRS ${GLPK_INCLUDE_DIR})
endif()
