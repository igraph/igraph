include(CheckCCompilerFlag)

# Enable POSIX features. This needs to be set here instead of in source files so
# that it affects CMake-based feature tests.
#
# See:
#  - https://pubs.opengroup.org/onlinepubs/007904875/functions/xsh_chap02_02.html
#  - https://www.gnu.org/software/libc/manual/html_node/Feature-Test-Macros.html
add_compile_definitions(_POSIX_C_SOURCE=200809L)

if(MSVC)
  add_compile_options(/FS)
  add_compile_definitions(_CRT_SECURE_NO_WARNINGS) # necessary to compile for UWP
endif()

if(NOT MSVC)
  # Even though we will later use 'no-unknown-warning-option', we perform the test for
  # 'unknown-warning-option', without the 'no-' prefix. This is necessary because GCC
  # will accept any warning option starting with 'no-', and will not error, yet it still
  # prints a message about the unrecognized option.
  check_c_compiler_flag("-Wunknown-warning-option" COMPILER_SUPPORTS_UNKNOWN_WARNING_OPTION_FLAG)
endif()

set(
  IGRAPH_WARNINGS_AS_ERRORS ON CACHE BOOL
  "Treat warnings as errors with GCC-like compilers"
)

option(FORCE_COLORED_OUTPUT "Always produce ANSI-colored output (GNU/Clang only)." FALSE)
if(FORCE_COLORED_OUTPUT)
  if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
   add_compile_options(-fdiagnostics-color=always)
  elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
   add_compile_options(-fcolor-diagnostics)
  elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang")
   add_compile_options(-fcolor-diagnostics)
  endif()
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
      /wd4221 # nonstandard extension used: '...': cannot be initialized using address of automatic variable '...'
      /wd4127 # conditional expression is constant
      /wd4702 # unreachable code
    )
  else()
    # Notes:
    # GCC does not complain when encountering an unsupported "no"-prefixed wanring option such as -Wno-foo.
    # Clang does complain, but these complaints can be silenced with -Wno-unknown-warning-option.
    # Therefore it is generally safe to use -Wno-... options that are only supported by recent GCC/Clang.
    target_compile_options(${TARGET_NAME} PRIVATE
      # GCC-style compilers:
      $<$<C_COMPILER_ID:GCC,Clang,AppleClang,Intel,IntelLLVM>:
        $<$<BOOL:${IGRAPH_WARNINGS_AS_ERRORS}>:-Werror>
        -Wall -Wextra -pedantic
        -Wstrict-prototypes
        -Wno-unused-function -Wno-unused-parameter -Wno-unused-but-set-variable -Wno-sign-compare -Wno-constant-logical-operand
      >
      $<$<BOOL:${COMPILER_SUPPORTS_UNKNOWN_WARNING_OPTION_FLAG}>:-Wno-unknown-warning-option>
      # Intel compiler:
      $<$<C_COMPILER_ID:Intel>:
        # disable #279: controlling expression is constant; affecting assert(condition && "message")
        # disable #592: variable "var" is used before its value is set; affecting IGRAPH_UNUSED
        -wd279 -wd592 -diag-disable=remark
      >
      # Intel LLVM:
      $<$<C_COMPILER_ID:IntelLLVM>:
        -fp-model=precise # The default 'fast' mode is not compatible with igraph's extensive use of NaN/Inf
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
