# Custom CPack install script that allows us to whitelist files to be copied
# to the tarball from the root directory, instead of copying the entire root
# directory recursively

if(CPACK_SOURCE_INSTALLED_DIRECTORIES)
	list(GET CPACK_BUILD_SOURCE_DIRS 0 SOURCE_DIR)
    # This branch runs only if CPack generates the source package, and within
	# this branch, CMAKE_CURRENT_BINARY_DIR refers to the root of the staging
	# area where the tarball is assembled
	file(GLOB FILES_TO_COPY "${SOURCE_DIR}/*.md")
	file(
		INSTALL ${FILES_TO_COPY}
		DESTINATION "${CMAKE_CURRENT_BINARY_DIR}"
	)
	file(
		INSTALL
		"${SOURCE_DIR}/AUTHORS"
		"${SOURCE_DIR}/CMakeLists.txt"
		"${SOURCE_DIR}/COPYING"
		"${SOURCE_DIR}/ChangeLog"
		"${SOURCE_DIR}/IGRAPH_VERSION"
		"${SOURCE_DIR}/INSTALL"
		"${SOURCE_DIR}/NEWS"
		"${SOURCE_DIR}/ONEWS"
		"${SOURCE_DIR}/src/config.h.in"
		"${SOURCE_DIR}/igraph.pc.cmake.in"
		DESTINATION "${CMAKE_CURRENT_BINARY_DIR}"
	)
endif()

