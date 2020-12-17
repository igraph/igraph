include(CMakeParseArguments)

find_program(DIFF_TOOL diff)
if(NOT DIFF_TOOL)
  find_program(FC_TOOL fc)
endif()

function(add_legacy_test FOLDER NAME)
  add_executable(test_${NAME} EXCLUDE_FROM_ALL ${CMAKE_SOURCE_DIR}/${FOLDER}/${NAME})
  add_dependencies(build_tests test_${NAME})
  target_link_libraries(test_${NAME} PRIVATE igraph)

  if (WIN32 AND NOT BUILD_SHARED_LIBS)
    # Add a compiler definition required to compile igraph in static mode on Windows
    target_compile_definitions(test_${NAME} PRIVATE IGRAPH_STATIC)
  endif()

  # Some tests depend on internal igraph headers so we also have to add src/
  # to the include path even though it's not part of the public API
  target_include_directories(
    test_${NAME} PRIVATE ${CMAKE_SOURCE_DIR}/src ${CMAKE_SOURCE_DIR}/vendor ${CMAKE_BINARY_DIR}/src
  )

  if (MSVC)
    # Add MSVC-specific include path for some headers that are missing on Windows
    target_include_directories(test_${NAME} PRIVATE ${CMAKE_SOURCE_DIR}/msvc/include)
  endif()

  set(EXPECTED_OUTPUT_FILE ${CMAKE_SOURCE_DIR}/${FOLDER}/${NAME}.out)
  set(OBSERVED_OUTPUT_FILE ${CMAKE_CURRENT_BINARY_DIR}/test_${NAME}.out)
  set(DIFF_FILE ${CMAKE_CURRENT_BINARY_DIR}/test_${NAME}.diff)
  get_filename_component(WORK_DIR ${EXPECTED_OUTPUT_FILE} DIRECTORY)

  if(EXISTS ${EXPECTED_OUTPUT_FILE})
    add_test(
      NAME ${NAME}
      COMMAND ${CMAKE_COMMAND}
        -DTEST_EXECUTABLE=$<TARGET_FILE:test_${NAME}>
        -DEXPECTED_OUTPUT_FILE=${EXPECTED_OUTPUT_FILE}
        -DOBSERVED_OUTPUT_FILE=${OBSERVED_OUTPUT_FILE}
        -DDIFF_FILE=${DIFF_FILE}
        -DDIFF_TOOL=${DIFF_TOOL}
        -DFC_TOOL=${FC_TOOL}
        -DIGRAPH_VERSION=${PACKAGE_VERSION}
        -P ${CMAKE_SOURCE_DIR}/etc/cmake/run_legacy_test.cmake
    )
    set_property(TEST ${NAME} PROPERTY SKIP_REGULAR_EXPRESSION "Test skipped")
  else()
    add_test(
      NAME ${NAME}
      COMMAND test_${NAME}
      WORKING_DIRECTORY ${WORK_DIR}
    )
    set_property(TEST ${NAME} PROPERTY SKIP_RETURN_CODE 77)
  endif()
  if (WIN32 AND BUILD_SHARED_LIBS)
    # On Windows the built igraph.dll is not automatically found by the tests. We therefore
    # add the dir that contains the built igraph.dll to the path environment variable
    # so that igraph.dll is found when running the tests.
    SET(IGRAPH_LIBDIR $<TARGET_FILE_DIR:igraph>)

    # The next line is necessitated by MinGW on Windows. MinGW uses forward slashes in
    # IGRAPH_LIBDIR, but we need to supply CTest with backslashes because CTest is executed
    # in a cmd.exe shell. So we simply replace forward slashes with backslases in
    # IGRAPH_LIBDIR.
    string(REPLACE "/" "\\" IGRAPH_LIBDIR ${IGRAPH_LIBDIR})

    # Semicolons are used as list separators in CMake so we need to escape them in the PATH,
    # otherwise the PATH envvar gets split by CMake before it passes the PATH on to CTest.
    string(JOIN "\;" CORRECT_PATH $ENV{PATH})
    SET_TESTS_PROPERTIES( ${NAME} PROPERTIES ENVIRONMENT "PATH=${IGRAPH_LIBDIR}\;${CORRECT_PATH}" )
  endif()
endfunction()

function(add_legacy_tests)
  cmake_parse_arguments(
    PARSED "" "FOLDER" "NAMES;LIBRARIES" ${ARGN}
  )
  foreach(NAME ${PARSED_NAMES})
    add_legacy_test(${PARSED_FOLDER} ${NAME})
    if(PARSED_LIBRARIES)
      target_link_libraries(test_${NAME} PRIVATE ${PARSED_LIBRARIES})
    endif()
  endforeach()
endfunction()
