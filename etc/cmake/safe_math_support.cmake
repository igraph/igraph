include(CheckCXXSourceCompiles)

# Check whether the compiler supports the __builtin_add_overflow() and __builtin_mul_overflow()
# builtins. These are present in recent GCC-compatible compilers.
cmake_push_check_state(RESET)

check_cxx_source_compiles("
    int main(void) {
        long long a=1, b=2, c;
        __builtin_add_overflow(a, b, &c);
        __builtin_mul_overflow(a, b, &c);
        return 0;
    }
    "
    HAVE_BUILTIN_OVERFLOW
)

cmake_pop_check_state()
