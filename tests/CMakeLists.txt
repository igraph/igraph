function(add_simple_test NAME)
  add_executable(test_${NAME} ${CMAKE_SOURCE_DIR}/examples/simple/${NAME}.c)
  target_link_libraries(test_${NAME} PRIVATE igraph)

  # Some tests depend on internal igraph headers so we also have to add src/
  # to the include path even though it's not part of the public API
  target_include_directories(
    test_${NAME} PRIVATE ${CMAKE_SOURCE_DIR}/src ${CMAKE_BINARY_DIR}/src
  )

  set(EXPECTED_OUTPUT_FILE ${CMAKE_SOURCE_DIR}/examples/simple/${NAME}.out)
  set(OBSERVED_OUTPUT_FILE ${CMAKE_CURRENT_BINARY_DIR}/test_${NAME}.out)

  if(EXISTS ${EXPECTED_OUTPUT_FILE})
    add_test(
      NAME ${NAME}
      COMMAND ${CMAKE_COMMAND}
        -DTEST_EXECUTABLE=$<TARGET_FILE:test_${NAME}>
        -DEXPECTED_OUTPUT_FILE=${EXPECTED_OUTPUT_FILE}
        -DOBSERVED_OUTPUT_FILE=${OBSERVED_OUTPUT_FILE}
        -P ${CMAKE_SOURCE_DIR}/etc/cmake/run_legacy_test.cmake
    )
  else()
    add_test(NAME ${NAME} COMMAND test_${NAME})
  endif()
endfunction()

function(add_simple_tests)
  foreach(NAME ${ARGV})
    add_simple_test(${NAME})
  endforeach()
endfunction()

# version.at
add_simple_test(igraph_version)

# types.at
add_simple_tests(
  vector vector2 vector3
  matrix matrix2 matrix3
  dqueue
  vector_ptr
  stack
  heap indheap d_indheap
  igraph_strvector
  igraph_trie
  igraph_psumtree
  igraph_array
  igraph_hashtable
  igraph_i_cutheap
  igraph_set
  2wheap
  igraph_sparsemat igraph_sparsemat2 igraph_sparsemat3 igraph_sparsemat4
  igraph_sparsemat5 igraph_sparsemat6 igraph_sparsemat7 igraph_sparsemat8
  igraph_sparsemat9 igraph_sparsemat_is_symmetric igraph_sparsemat_minmax
  igraph_sparsemat_which_minmax
  spmatrix
  biguint
  igraph_marked_queue
  igraph_complex
)

# basic.at
add_simple_tests(
  igraph_empty
  igraph_copy
  igraph_add_edges
  igraph_add_vertices
  igraph_delete_edges
  igraph_delete_vertices
)
