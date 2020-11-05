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

# Helper function to add preprocesor definition of IGRAPH_FILE_BASENAME
# to pass the filename without directory path for debugging use.
#
# Example:
#
#   define_file_basename_for_sources(my_target)
#
# Will add -DIGRAPH_FILE_BASENAME="filename" for each source file depended
# on by my_target, where filename is the name of the file.
#
# Source: https://stackoverflow.com/a/27990434/156771
function(define_file_basename_for_sources targetname)
  get_target_property(source_files "${targetname}" SOURCES)
  foreach(sourcefile ${source_files})
    # Add the IGRAPH_FILE_BASENAME=filename compile definition to the list.
    get_filename_component(basename "${sourcefile}" NAME)
    # Set the updated compile definitions on the source file.
    set_property(
      SOURCE "${sourcefile}" APPEND
      PROPERTY COMPILE_DEFINITIONS "IGRAPH_FILE_BASENAME=\"${basename}\""
    )
  endforeach()
endfunction()

