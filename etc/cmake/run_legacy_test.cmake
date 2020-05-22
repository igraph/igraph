# Runs a legacy autotools-based test with a file containing the expected output
#
# Parameters of the script:
#
# - TEST_EXECUTABLE: full path of the compiled test executable
# - EXPECTED_OUTPUT_FILE: full path of the file containing the expected output
# - OBSERVED_OUTPUT_FILE: full path of the fhile where the observed output
#   can be written

execute_process(
  COMMAND ${TEST_EXECUTABLE}
  RESULT_VARIABLE ERROR_CODE
  OUTPUT_VARIABLE OBSERVED_OUTPUT
)

if(ERROR_CODE)
  message(FATAL_ERROR "Test exited abnormally with error: ${ERROR_CODE}")
endif()

file(WRITE ${OBSERVED_OUTPUT_FILE} ${OBSERVED_OUTPUT})

execute_process(
  COMMAND ${CMAKE_COMMAND} -E compare_files
  ${EXPECTED_OUTPUT_FILE} ${OBSERVED_OUTPUT_FILE}
  RESULT_VARIABLE ARE_DIFFERENT
)

if(ARE_DIFFERENT)
  message(FATAL_ERROR "Test case output differs from the expected output")
endif()

file(REMOVE ${OBSERVED_OUTPUT_FILE})
