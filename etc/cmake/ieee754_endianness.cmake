include(CheckCSourceRuns)

cmake_push_check_state(RESET)

# Check whether IEEE754 doubles are laid out in little-endian order. We do this
# only when not cross-compiling; during cross-compilation, the host architecture
# might have different endianness conventions than the target, and we are running
# the test on the host here
if(CMAKE_CROSSCOMPILING AND NOT CMAKE_CROSSCOMPILING_EMULATOR)
    # If we are cross-compiling and we have no emulator, let's just assume that
    # IEEE754 doubles use the same endianness as uint64_t
    set(IEEE754_DOUBLE_ENDIANNESS_MATCHES YES)
    message(WARNING "\
igraph is being cross-compiled, therefore we cannot validate whether the \
endianness of IEEE754 doubles is the same as the endianness of uint64_t. \
Most likely it is, unless you are compiling for some esoteric platform, \
in which case you need make sure that this is the case on your own.\
")
else()
    if(NOT DEFINED CACHE{IEEE754_DOUBLE_ENDIANNESS_MATCHES})
        try_run(
            IEEE754_DOUBLE_ENDIANNESS_TEST_EXIT_CODE
            IEEE754_DOUBLE_ENDIANNESS_TEST_COMPILES
            ${CMAKE_BINARY_DIR}
            ${PROJECT_SOURCE_DIR}/etc/cmake/ieee754_endianness_check.c
            RUN_OUTPUT_VARIABLE IEEE754_DOUBLE_ENDIANNESS_TEST_RESULT
        )
        # Strip trailing newline, which is necessary on some platforms (such as node.js)
        # to complete printing the output.
        string(STRIP "${IEEE754_DOUBLE_ENDIANNESS_TEST_RESULT}" IEEE754_DOUBLE_ENDIANNESS_TEST_RESULT)
        if(IEEE754_DOUBLE_ENDIANNESS_TEST_EXIT_CODE EQUAL 0)
            if(IEEE754_DOUBLE_ENDIANNESS_TEST_RESULT STREQUAL "OK")
                set(TEST_RESULT YES)
            else()
                set(TEST_RESULT NO)
            endif()
        else()
            message(FATAL_ERROR "IEEE754 double endianness test terminated abnormally.")
        endif()

        set(
            IEEE754_DOUBLE_ENDIANNESS_MATCHES ${TEST_RESULT} CACHE BOOL
            "Specifies whether the endianness of IEEE754 doubles is the same as the endianness of uint64_t."
            FORCE
        )
        mark_as_advanced(IEEE754_DOUBLE_ENDIANNESS_MATCHES)
    endif()
endif()

cmake_pop_check_state()

if(NOT IEEE754_DOUBLE_ENDIANNESS_MATCHES)
    message(FATAL_ERROR "igraph only supports platforms where IEEE754 doubles have the same endianness as uint64_t.")
endif()
