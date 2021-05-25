include(CMakeParseArguments)

find_program(DIFF_TOOL diff)
if(NOT DIFF_TOOL)
  find_program(FC_TOOL fc)
endif()

function(add_legacy_test FOLDER NAME NAMESPACE)
  set(TARGET_NAME ${NAMESPACE}_${NAME})
  set(TEST_NAME "${NAMESPACE}::${NAME}")

  add_executable(${TARGET_NAME} EXCLUDE_FROM_ALL ${PROJECT_SOURCE_DIR}/${FOLDER}/${NAME}.c)
  add_dependencies(build_tests ${TARGET_NAME})
  target_link_libraries(${TARGET_NAME} PRIVATE igraph)

  if (NOT BUILD_SHARED_LIBS)
    # Add a compiler definition required to compile igraph in static mode
    target_compile_definitions(${TARGET_NAME} PRIVATE IGRAPH_STATIC)
  endif()

  # Some tests depend on internal igraph headers so we also have to add src/
  # to the include path even though it's not part of the public API
  target_include_directories(
    ${TARGET_NAME} PRIVATE ${CMAKE_SOURCE_DIR}/src ${CMAKE_SOURCE_DIR}/vendor ${CMAKE_BINARY_DIR}/src
  )

  if (MSVC)
    # Add MSVC-specific include path for some headers that are missing on Windows
    target_include_directories(${TARGET_NAME} PRIVATE ${CMAKE_SOURCE_DIR}/msvc/include)
  endif()

  set(EXPECTED_OUTPUT_FILE ${CMAKE_SOURCE_DIR}/${FOLDER}/${NAME}.out)
  set(OBSERVED_OUTPUT_FILE ${CMAKE_CURRENT_BINARY_DIR}/${TARGET_NAME}.out)
  set(DIFF_FILE ${CMAKE_CURRENT_BINARY_DIR}/${TARGET_NAME}.diff)
  get_filename_component(WORK_DIR ${EXPECTED_OUTPUT_FILE} DIRECTORY)

  if(EXISTS ${EXPECTED_OUTPUT_FILE})
    add_test(
      NAME ${TEST_NAME}
      COMMAND ${CMAKE_COMMAND}
        -DTEST_EXECUTABLE=$<TARGET_FILE:${TARGET_NAME}>
        -DEXPECTED_OUTPUT_FILE=${EXPECTED_OUTPUT_FILE}
        -DOBSERVED_OUTPUT_FILE=${OBSERVED_OUTPUT_FILE}
        -DDIFF_FILE=${DIFF_FILE}
        -DDIFF_TOOL=${DIFF_TOOL}
        -DFC_TOOL=${FC_TOOL}
        -DIGRAPH_VERSION=${PACKAGE_VERSION}
        -P ${CMAKE_SOURCE_DIR}/etc/cmake/run_legacy_test.cmake
    )
    set_property(TEST ${TEST_NAME} PROPERTY SKIP_REGULAR_EXPRESSION "Test skipped")
  else()
    add_test(
      NAME ${TEST_NAME}
      COMMAND ${TARGET_NAME}
      WORKING_DIRECTORY ${WORK_DIR}
    )
    set_property(TEST ${TEST_NAME} PROPERTY SKIP_RETURN_CODE 77)
  endif()
  if (WIN32 AND BUILD_SHARED_LIBS)
    # On Windows the built igraph.dll is not automatically found by the tests. We therefore
    # add the dir that contains the built igraph.dll to the path environment variable
    # so that igraph.dll is found when running the tests.
    SET(IGRAPH_LIBDIR $<TARGET_FILE_DIR:igraph>)

    # The next line is necessitated by MinGW on Windows. MinGW uses forward slashes in
    # IGRAPH_LIBDIR, but we need to supply CTest with backslashes because CTest is executed
    # in a cmd.exe shell. We therefore explicitly ensure that that path is transformed to a
    # native path.
    file(TO_NATIVE_PATH "${IGRAPH_LIBDIR}" IGRAPH_LIBDIR)

    # Semicolons are used as list separators in CMake so we need to escape them in the PATH,
    # otherwise the PATH envvar gets split by CMake before it passes the PATH on to CTest.
    # We process each path separately to ensure it is a proper path. In particular, we need
    # to ensure that a trailing backslash is not incorrectly interpreted as an escape
    # character. Presumably, with cmake 3.20, this can be changed to using TO_NATIVE_PATH_LIST.
    SET(TEST_PATHS)
    foreach (PATH $ENV{PATH})
      file(TO_NATIVE_PATH "${PATH}" CORRECT_PATH)
      # Remove trailing backslash
      STRING(REGEX REPLACE "\\$" "" CORRECT_PATH ${CORRECT_PATH})
      list(APPEND TEST_PATHS ${CORRECT_PATH})
    endforeach()

    # Join all paths in a single string, separated by an escaped semi-colon.
    string(JOIN "\;" CORRECT_PATHS ${TEST_PATHS})
    SET_TESTS_PROPERTIES(${TEST_NAME} PROPERTIES ENVIRONMENT "PATH=${IGRAPH_LIBDIR}\;${CORRECT_PATHS}" )
  endif()
endfunction()

function(add_legacy_tests)
  cmake_parse_arguments(
    PARSED "" "FOLDER" "NAMES;LIBRARIES" ${ARGN}
  )
  foreach(NAME ${PARSED_NAMES})
    add_legacy_test(${PARSED_FOLDER} ${NAME} test)
    if(PARSED_LIBRARIES)
      target_link_libraries(test_${NAME} PRIVATE ${PARSED_LIBRARIES})
    endif()
  endforeach()
endfunction()

function(add_examples)
  cmake_parse_arguments(
    PARSED "" "FOLDER" "NAMES;LIBRARIES" ${ARGN}
  )
  foreach(NAME ${PARSED_NAMES})
    add_legacy_test(${PARSED_FOLDER} ${NAME} example)
    if(PARSED_LIBRARIES)
      target_link_libraries(example_${NAME} PRIVATE ${PARSED_LIBRARIES})
    endif()
  endforeach()
endfunction()
