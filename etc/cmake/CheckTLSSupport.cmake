include(CheckCSourceCompiles)
include(CMakePushCheckState)

macro(check_tls_support VAR)
  if(NOT DEFINED "${VAR}")
    cmake_push_check_state()
    set(CMAKE_REQUIRED_QUIET 1)

    check_c_source_compiles("
    __thread int tls;

    int main(void) {
        return 0;
    }" HAVE_GCC_TLS)

    if(HAVE_GCC_TLS)
      message(STATUS "Thread-local storage: supported (__thread)")
      set(${VAR} "__thread" CACHE INTERNAL "Thread-local storage support keyword in compiler")
    else()
      check_c_source_compiles("
      __declspec(thread) int tls;

      int main(void) {
          return 0;
      }" HAVE_MSVC_TLS)
      if(HAVE_MSVC_TLS)
        message(STATUS "Thread-local storage: supported (__declspec(thread))")
        set(${VAR} "__declspec(thread)" CACHE INTERNAL "Thread-local storage keyword in compiler")
      else()
        message(STATUS "Thread-local storage: not supported")
        set(${VAR} "" CACHE INTERNAL "Thread-local storage keyword in compiler")
      endif()
    endif()
    cmake_pop_check_state()
  endif()
endmacro()
