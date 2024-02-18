
function(add_fuzzer NAME)
  set(TARGET_NAME fuzzer_${NAME})

  add_executable(${TARGET_NAME} EXCLUDE_FROM_ALL ${PROJECT_SOURCE_DIR}/fuzzing/${NAME}.cpp)

  add_dependencies(build_fuzzers ${TARGET_NAME})

  target_link_libraries(${TARGET_NAME} PRIVATE igraph)

  # The -fsanitize=fuzzer-no-link is already added by the top-level CMakeLists.txt
  # for general fuzzer instrumentation. Additionally, we need -fsanitize=fuzzer
  # for the fuzz targets, which do not contain a main() function, to link in the
  # fuzz driver. See https://llvm.org/docs/LibFuzzer.html
  target_link_options(${TARGET_NAME} PRIVATE -fsanitize=fuzzer)
endfunction()
