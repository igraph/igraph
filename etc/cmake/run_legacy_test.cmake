# Runs a legacy autotools-based test with a file containing the expected output
#
# Parameters of the script:
#
# - TEST_EXECUTABLE: full path of the compiled test executable
# - EXPECTED_OUTPUT_FILE: full path of the file containing the expected output
# - OBSERVED_OUTPUT_FILE: full path of the file where the observed output
#   can be written
# - IGRAPH_VERSION: version string of igraph that should be replaced in
#   expected outputs

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
  message(FATAL_ERROR ${MESSAGE})
else()
  string(REPLACE ${IGRAPH_VERSION} "\@VERSION\@" OBSERVED_OUTPUT "${OBSERVED_OUTPUT}")
  file(WRITE ${OBSERVED_OUTPUT_FILE} "${OBSERVED_OUTPUT}")

  execute_process(
    COMMAND ${CMAKE_COMMAND} -E compare_files --ignore-eol
    ${EXPECTED_OUTPUT_FILE} ${OBSERVED_OUTPUT_FILE}
    RESULT_VARIABLE ARE_DIFFERENT
  )

  if(ARE_DIFFERENT)
    message(FATAL_ERROR "Test case output differs from the expected output")
  endif()

  file(REMOVE ${OBSERVED_OUTPUT_FILE})
endif()
