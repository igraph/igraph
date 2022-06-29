include(CheckCSourceRuns)

cmake_push_check_state(RESET)

# Check whether IEEE754 doubles are laid out in little-endian order. We do this
# only when not cross-compiling; during cross-compilation, the host architecture
# might have different endianness conventions than the target, and we are running
# the test on the host here
if(NOT CMAKE_CROSSCOMPILING)
    if(NOT DEFINED CACHE{IEEE754_DOUBLE_ENDIANNESS_MATCHES})
        try_run(
            IEEE754_DOUBLE_ENDIANNESS_TEST_EXIT_CODE
            IEEE754_DOUBLE_ENDIANNESS_TEST_COMPILES
            ${CMAKE_BINARY_DIR}
            ${PROJECT_SOURCE_DIR}/etc/cmake/ieee754_endianness_check.c
        )
        if(IEEE754_DOUBLE_ENDIANNESS_TEST_EXIT_CODE EQUAL 0)
            set(TEST_RESULT YES)
        else()
            set(TEST_RESULT NO)
        endif()

		set(
			IEEE754_DOUBLE_ENDIANNESS_MATCHES ${TEST_RESULT} CACHE BOOL
			"Specifies whether the endianness of IEEE754 doubles is the same as the endianness of uint64_t."
			FORCE
		)
		mark_as_advanced(IEEE754_DOUBLE_ENDIANNESS_MATCHES)
	endif()
else()
    # If we are cross-compiling, let's just assume that IEEE754 doubles use the
    # same endianness as uint64_t
    set(IEEE754_DOUBLE_ENDIANNESS_MATCHES YES)
endif()

cmake_pop_check_state()

if(NOT IEEE754_DOUBLE_ENDIANNESS_MATCHES)
    message(FATAL_ERROR "igraph only supports platforms where IEEE754 doubles have the same endianness as uint64_t")
endif()
