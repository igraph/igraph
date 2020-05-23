function(add_simple_test NAME)
  add_executable(test_${NAME} ${CMAKE_SOURCE_DIR}/examples/simple/${NAME})
  target_link_libraries(test_${NAME} PRIVATE igraph)

  # Some tests depend on internal igraph headers so we also have to add src/
  # to the include path even though it's not part of the public API
  target_include_directories(
    test_${NAME} PRIVATE ${CMAKE_SOURCE_DIR}/src ${CMAKE_BINARY_DIR}/src
  )

  set(EXPECTED_OUTPUT_FILE ${CMAKE_SOURCE_DIR}/examples/simple/${NAME}.out)
  set(OBSERVED_OUTPUT_FILE ${CMAKE_CURRENT_BINARY_DIR}/test_${NAME}.out)
  get_filename_component(WORK_DIR ${EXPECTED_OUTPUT_FILE} DIRECTORY)

  if(EXISTS ${EXPECTED_OUTPUT_FILE})
    # TODO: this is not entirely correct; tests that are skipped because of their
    # return codes are marked as "successful"
    add_test(
      NAME ${NAME}
      COMMAND ${CMAKE_COMMAND}
        -DTEST_EXECUTABLE=$<TARGET_FILE:test_${NAME}>
        -DEXPECTED_OUTPUT_FILE=${EXPECTED_OUTPUT_FILE}
        -DOBSERVED_OUTPUT_FILE=${OBSERVED_OUTPUT_FILE}
        -DIGRAPH_VERSION=${PACKAGE_VERSION}
        -P ${CMAKE_SOURCE_DIR}/etc/cmake/run_legacy_test.cmake
    )
  else()
    add_test(
      NAME ${NAME}
      COMMAND test_${NAME}
      WORKING_DIRECTORY ${WORK_DIR}
    )
  endif()

  set_tests_properties(${NAME} PROPERTIES SKIP_RETURN_CODE 77)
endfunction()

function(add_simple_tests)
  foreach(NAME ${ARGV})
    add_simple_test(${NAME})
  endforeach()
endfunction()
