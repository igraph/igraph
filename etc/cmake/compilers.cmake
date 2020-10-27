if(MSVC)
  add_compile_options(/FS)
else()
  add_compile_options(
    # Intel compiler:
    # Do not print remarks:
    $<$<C_COMPILER_ID:Intel>:-diag-disable=remark>
  )
endif()

macro(use_all_warnings TARGET_NAME)
  if(MSVC)
    target_compile_options(${TARGET_NAME} PRIVATE /W4)
  else()
    target_compile_options(${TARGET_NAME} PRIVATE 
      -Wall -Wextra -pedantic -Werror -Wno-unused-parameter -Wno-sign-compare
      $<$<NOT:$<C_COMPILER_ID:Intel>>:-Wno-varargs>
      # Intel compiler:
      # disable #279: controlling expression is constant; affecting assert(condition && "message")
      # disable #188: enumerated type mixed with another type; affecting IGRAPH_CHECK
      # disable #592: variable "var" is used before its value is set; affecting IGRAPH_UNUSED
      $<$<C_COMPILER_ID:Intel>:-wd279 -wd188 -wd592>
    )
  endif()
endmacro()

