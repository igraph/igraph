include(CheckCXXSourceCompiles)
include(CheckTypeSize)

cmake_push_check_state(RESET)

# Check whether the compiler supports the _umul128() intrinsic
check_cxx_source_compiles("
    #include <intrin.h>

    int main(void) {
        unsigned long long a = 0, b = 0;
        unsigned long long c;
        volatile unsigned long long d;
        d = _umul128(a, b, &c);
        return 0;
    }
    "
    HAVE__UMUL128
)

# Check whether the compiler supports the __umulh() intrinsic
check_cxx_source_compiles("
    #include <intrin.h>

    int main(void) {
        unsigned long long a = 0, b = 0;
        volatile unsigned long long c;
        c = __umulh(a, b);
        return 0;
    }
    "
    HAVE___UMULH
)

# Check whether the compiler has __uint128_t
check_type_size("__uint128_t" UINT128 LANGUAGE CXX)
if(UINT128 EQUAL 16)
    set(HAVE___UINT128_T ON)
else()
    set(HAVE___UINT128_T OFF)
endif()

cmake_pop_check_state()
