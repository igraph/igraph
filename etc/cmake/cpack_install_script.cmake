# Custom CPack install script that allows us to whitelist files to be copied
# to the tarball from the root directory, instead of copying the entire root
# directory recursively

if(CPACK_SOURCE_INSTALLED_DIRECTORIES)
    # Make sure that the parser sources are built
    execute_process(
        COMMAND "${CMAKE_COMMAND}"
        --build "${CPACK_PACKAGE_DIRECTORY}"
        --target parsersources
        RESULT_VARIABLE EXIT_CODE
    )
    if(NOT EXIT_CODE EQUAL 0)
        message(FATAL_ERROR "Failed to build the parser sources.")
    endif()

    # Generate a version file in the build folder if we don't have one in the
    # source folder
    if(EXISTS "${SOURCE_DIR}/IGRAPH_VERSION")
        set(IGRAPH_VERSION_FILE "${SOURCE_DIR}/IGRAPH_VERSION")
    else()
        execute_process(
            COMMAND "${CMAKE_COMMAND}"
            --build "${CPACK_PACKAGE_DIRECTORY}"
            --target versionfile
            RESULT_VARIABLE EXIT_CODE
        )
        if(NOT EXIT_CODE EQUAL 0)
            message(FATAL_ERROR "Failed to determine the version number of igraph that is being packaged.")
        endif()
        set(IGRAPH_VERSION_FILE "${CPACK_PACKAGE_DIRECTORY}/IGRAPH_VERSION")
    endif()

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
        "${SOURCE_DIR}/INSTALL"
        "${SOURCE_DIR}/NEWS"
        "${SOURCE_DIR}/ONEWS"
        "${SOURCE_DIR}/igraph.pc.in"
        "${IGRAPH_VERSION_FILE}"
        DESTINATION "${CMAKE_CURRENT_BINARY_DIR}"
    )
    file(
        INSTALL
        "${SOURCE_DIR}/src/config.h.in"
        DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/src"
    )
    file(
		INSTALL
		"${CPACK_PACKAGE_DIRECTORY}/src/io/parsers"
		DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/src/io"
	)
    file(
		INSTALL
		"${CPACK_PACKAGE_DIRECTORY}/doc/html"
		DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/doc"
	)
endif()

