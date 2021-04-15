include(CheckCCompilerFlag)

if(MSVC)
  add_compile_options(/FS)
endif()

if (NOT MSVC)
  check_c_compiler_flag("-Wno-varargs" COMPILER_SUPPORTS_NO_VARARGS_FLAG)
endif()

macro(use_all_warnings TARGET_NAME)
  if(MSVC)
    target_compile_options(${TARGET_NAME} PRIVATE 
      /W4 # enable most warnings, then disable:
      /wd4244 # 'conversion' conversion from 'type1' to 'type2', possible loss of data
      /wd4267 # 'var' : conversion from 'size_t' to 'type', possible loss of data 
      /wd4996 # deprecated functions, e.g. 'sprintf': This function or variable may be unsafe. Consider using sprintf_s instead.
      /wd4456 # declaration of 'identifier' hides previous local declaration
      /wd4800 # forcing value to 'true' or 'false' (performance warning)
      /wd4204 # nonstandard extension used: non-constant aggregate initializer
      /wd4701 # potentially uninitialized local variable
    )
  else()
    target_compile_options(${TARGET_NAME} PRIVATE 
      # GCC-style compilers:
      $<$<C_COMPILER_ID:GCC,Clang,AppleClang,Intel>:
        -Wall -Wextra -pedantic -Werror -Wno-unused-function -Wno-unused-parameter -Wno-sign-compare
      >
      $<$<BOOL:${COMPILER_SUPPORTS_NO_VARARGS_FLAG}>:-Wno-varargs>
      # Intel compiler:
      $<$<C_COMPILER_ID:Intel>:
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
  get_target_property(source_dir "${targetname}" SOURCE_DIR)
  foreach(sourcefile ${source_files})
    # Turn relative paths into absolute
    get_filename_component(source_full_path "${sourcefile}" ABSOLUTE BASE_DIR "${source_dir}")

    # Figure out whether the relative path from the source or the build folder
    # is shorter
    file(RELATIVE_PATH source_rel_path "${PROJECT_SOURCE_DIR}" "${source_full_path}")
    file(RELATIVE_PATH binary_rel_path "${PROJECT_BINARY_DIR}" "${source_full_path}")
    string(LENGTH "${source_rel_path}" source_rel_path_length)
    string(LENGTH "${binary_rel_path}" binary_rel_path_length)
    if(binary_rel_path_length LESS source_rel_path_length)
      set(basename "${binary_rel_path}")
    else()
      set(basename "${source_rel_path}")
    endif()

    # Add the IGRAPH_FILE_BASENAME=filename compile definition to the source file
    set_property(
      SOURCE "${sourcefile}" APPEND
      PROPERTY COMPILE_DEFINITIONS "IGRAPH_FILE_BASENAME=\"${basename}\""
    )
  endforeach()
endfunction()

