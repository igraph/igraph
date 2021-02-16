# Runs a legacy autotools-based test with a file containing the expected output
#
# Parameters of the script:
#
# - TEST_EXECUTABLE: full path of the compiled test executable
# - EXPECTED_OUTPUT_FILE: full path of the file containing the expected output
# - OBSERVED_OUTPUT_FILE: full path of the file where the observed output
#   can be written
# - DIFF_FILE: full path of the file where the differences between the expectd
#   and the observed output should be written
# - DIFF_TOOL: full path to a "diff" tool on the system of the user, if present
# - FC_TOOL: full path to a "fc" tool on the system of the user, if present
# - IGRAPH_VERSION: version string of igraph that should be replaced in
#   expected outputs

function(print_file FILENAME)
  # Replacement of "cmake -E cat" for older CMake versions. cat was added in
  # CMake 3.18
  file(TO_NATIVE_PATH "${FILENAME}" FILENAME_NATIVE)
  if(UNIX OR APPLE)
    # Most likely Linux or macOS
    execute_process(COMMAND "/bin/sh" "-c" "cat ${FILENAME_NATIVE}")
  elseif(WIN32)
    # Most likely Windows
    execute_process(COMMAND "cmd" "/c" "type" "${FILENAME_NATIVE}")
  endif()
endfunction()

get_filename_component(WORK_DIR ${EXPECTED_OUTPUT_FILE} DIRECTORY)

execute_process(
  COMMAND ${TEST_EXECUTABLE}
  WORKING_DIRECTORY ${WORK_DIR}
  RESULT_VARIABLE ERROR_CODE
  OUTPUT_VARIABLE OBSERVED_OUTPUT
)

if(ERROR_CODE EQUAL 77)
  message(STATUS "Test skipped")
elseif(ERROR_CODE)
  set(MESSAGE "Test exited abnormally with error: ${ERROR_CODE}")
  file(WRITE ${OBSERVED_OUTPUT_FILE} "${MESSAGE}\n=========================================\n${OBSERVED_OUTPUT}")
  print_file("${OBSERVED_OUTPUT_FILE}")
  file(REMOVE ${DIFF_FILE})
  message(FATAL_ERROR "Exiting test.")
else()
  string(REPLACE ${IGRAPH_VERSION} "\@VERSION\@" OBSERVED_OUTPUT "${OBSERVED_OUTPUT}")
  file(WRITE ${OBSERVED_OUTPUT_FILE} "${OBSERVED_OUTPUT}")

  execute_process(
    COMMAND ${CMAKE_COMMAND} -E compare_files --ignore-eol
    ${EXPECTED_OUTPUT_FILE} ${OBSERVED_OUTPUT_FILE}
    RESULT_VARIABLE ARE_DIFFERENT
  )

  if(ARE_DIFFERENT)
    if(DIFF_TOOL)
      execute_process(
        COMMAND ${DIFF_TOOL} -u ${EXPECTED_OUTPUT_FILE} ${OBSERVED_OUTPUT_FILE}
        OUTPUT_FILE ${DIFF_FILE}
      )
    elseif(FC_TOOL)
      file(TO_NATIVE_PATH "${EXPECTED_OUTPUT_FILE}" REAL_EXPECTED_OUTPUT_FILE)
      file(TO_NATIVE_PATH "${OBSERVED_OUTPUT_FILE}" REAL_OBSERVED_OUTPUT_FILE)
      execute_process(
        COMMAND ${FC_TOOL} /A ${REAL_EXPECTED_OUTPUT_FILE} ${REAL_OBSERVED_OUTPUT_FILE}
        OUTPUT_FILE ${DIFF_FILE}
      )
    endif()

    message(STATUS "Test case output differs from the expected output")
    if(EXISTS ${DIFF_FILE})
      message(STATUS "See diff below:")
      message(STATUS "-------------------------------------------------------")
      print_file("${DIFF_FILE}")
      message(STATUS "-------------------------------------------------------")
    else()
      message(STATUS "Diff omitted; no diff tool was installed.")
    endif()
    message(FATAL_ERROR "Exiting test.")
  else()
    file(REMOVE ${DIFF_FILE})
  endif()

  file(REMOVE ${OBSERVED_OUTPUT_FILE})
endif()
