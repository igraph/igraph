if(MSVC)
  add_compile_options(/FS)
endif()

macro(use_all_warnings TARGET_NAME)
  if(MSVC)
    target_compile_options(${TARGET_NAME} PRIVATE /W4)
  else()
    target_compile_options(${TARGET_NAME} PRIVATE 
      # GCC-style compilers:
      $<$<C_COMPILER_ID:GCC,Clang,AppleClang>:
        -Wall -Wextra -pedantic -Werror -Wno-unused-parameter -Wno-sign-compare -Wno-varargs
      >
      # Intel compiler:
      $<$<C_COMPILER_ID:Intel>:
        -Wall -Wextra -pedantic -Werror -Wno-unused-parameter -Wno-sign-compare
        # disable #279: controlling expression is constant; affecting assert(condition && "message")
        # disable #188: enumerated type mixed with another type; affecting IGRAPH_CHECK
        # disable #592: variable "var" is used before its value is set; affecting IGRAPH_UNUSED
        -wd279 -wd188 -wd592 -diag-disable=remark
      >
    )
  endif()
endmacro()

