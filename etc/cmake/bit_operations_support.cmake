include(CheckCXXSourceCompiles)
include(CheckTypeSize)

cmake_push_check_state(RESET)

# Check whether the compiler supports the __popcnt64() intrinsic
check_cxx_source_compiles("
    #include <intrin.h>

    int main(void) {
        unsigned long long a = 0xDEADBEEF;
        volatile unsigned long long b;
        b = __popcnt64(a);
        return 0;
    }
    "
    HAVE__POPCNT64
)

# Check whether the compiler supports the __popcnt() intrinsic
check_cxx_source_compiles("
    #include <intrin.h>

    int main(void) {
        unsigned long a = 0xDEADBEEF;
        volatile unsigned long long b;
        b = __popcnt(a);
        return 0;
    }
    "
    HAVE__POPCNT
)

# Check whether the compiler supports the _BitScanForward64() intrinsic
check_cxx_source_compiles("
    #include <intrin.h>

    int main(void) {
        unsigned long long a = 0xDEADBEEF;
        unsinged long b;
        volatile unsigned long c;
        c = _BitScanForward64(&b, a) ? b : 64;
        return 0;
    }
    "
    HAVE__BITSCANFORWARD64
)

# Check whether the compiler supports the _BitScanForward() intrinsic
check_cxx_source_compiles("
    #include <intrin.h>

    int main(void) {
        unsigned long a = 0xDEADBEEF b;
        volatile unsigned long c;
        c = _BitScanForward(&b, a) ? b : 64;
        return 0;
    }
    "
    HAVE__BITSCANFORWARD
)

# Check whether the compiler supports the _BitScanReverse64() intrinsic
check_cxx_source_compiles("
    #include <intrin.h>

    int main(void) {
        unsigned long long a = 0xDEADBEEF;
        unsigned long b;
        volatile unsigned long c;
        c = _BitScanReverse64(&b, a) ? b : 64;
        return 0;
    }
    "
    HAVE__BITSCANREVERSE64
)

# Check whether the compiler supports the _BitScanReverse() intrinsic
check_cxx_source_compiles("
    #include <intrin.h>

    int main(void) {
        unsigned long a = 0xDEADBEEF, b;
        volatile unsigned long long c;
        c = _BitScanReverse(&b, a) ? b : 64;
        return 0;
    }
    "
    HAVE__BITSCANREVERSE
)

cmake_pop_check_state()
