set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "igraph library")
set(CPACK_PACKAGE_HOMEPAGE_URL "https://igraph.org")
set(CPACK_PACKAGE_VENDOR "The igraph development team")

set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/COPYING")

# Alias "dist" to "package_source"
add_custom_target(dist
  COMMAND "${CMAKE_COMMAND}"
    --build "${CMAKE_BINARY_DIR}"
    --target package_source
  VERBATIM
  USES_TERMINAL
)

#############################################################################
## Configuration of the source package
#############################################################################

# Set source package name and format
set(CPACK_SOURCE_PACKAGE_FILE_NAME "igraph-${CMAKE_PROJECT_VERSION}")
set(CPACK_SOURCE_GENERATOR "TGZ")

# Declare what to include in the source tarball. Unfortunately we can only
# declare full directories here, not individual files.
set(
	CPACK_SOURCE_INSTALLED_DIRECTORIES
	"${CMAKE_SOURCE_DIR}/examples;/examples"
	"${CMAKE_SOURCE_DIR}/include;/include"
	"${CMAKE_SOURCE_DIR}/optional;/optional/glpk"
	"${CMAKE_SOURCE_DIR}/src;/src"
	"${CMAKE_SOURCE_DIR}/tests;/tests"
)

# To add: igraph.pc.cmake.in, *.md, all documentation stuff from the root folder,
# all CMakeLists.txt files

# Ignore the build and all hidden folders. Also ignore obsolete autoconf-related
# stuff. The latter won't be needed once we fully transitioned to CMake.
set(
	CPACK_SOURCE_IGNORE_FILES
	"\\\\..*/"
	"${CMAKE_SOURCE_DIR}/build"
	"Makefile.am"
	"Makefile.in"
	"configure.ac"
	"atlocal.in"
	"\.at$"
)

#############################################################################
## Now we can include CPack
#############################################################################

include(CPack)

