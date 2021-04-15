include(CMakeParseArguments)

function(add_benchmark NAME NAMESPACE)
  set(TARGET_NAME ${NAMESPACE}_${NAME})

  add_executable(${TARGET_NAME} EXCLUDE_FROM_ALL ${PROJECT_SOURCE_DIR}/tests/benchmarks/${NAME}.c)
  add_dependencies(build_benchmarks ${TARGET_NAME})
  target_link_libraries(${TARGET_NAME} PRIVATE igraph)

  if (NOT BUILD_SHARED_LIBS)
    # Add a compiler definition required to compile igraph in static mode
    target_compile_definitions(${TARGET_NAME} PRIVATE IGRAPH_STATIC)
  endif()

  if (MSVC)
    # Add MSVC-specific include path for some headers that are missing on Windows
    target_include_directories(${TARGET_NAME} PRIVATE ${CMAKE_SOURCE_DIR}/msvc/include)
  endif()

  add_custom_command(
    TARGET benchmark
    POST_BUILD
    COMMAND ${TARGET_NAME}
    COMMENT "Running benchmark: ${NAME}"
  )
endfunction()

function(add_benchmarks)
  cmake_parse_arguments(
    PARSED "" "" "NAMES;LIBRARIES" ${ARGN}
  )
  foreach(NAME ${PARSED_NAMES})
    add_benchmark(${NAME} benchmark)
    if(PARSED_LIBRARIES)
      target_link_libraries(benchmark_${NAME} PRIVATE ${PARSED_LIBRARIES})
    endif()
  endforeach()
endfunction()
