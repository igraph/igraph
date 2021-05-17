# Declare the files needed to compile our vendored CXSparse copy
add_library(
  cxsparse_vendored
  OBJECT
  EXCLUDE_FROM_ALL
  cs_add.c
  cs_amd.c
  cs_chol.c
  cs_cholsol.c
  cs_compress.c
  cs_counts.c
  cs_cumsum.c
  cs_dfs.c
  cs_dmperm.c
  cs_droptol.c
  cs_dropzeros.c
  cs_dupl.c
  cs_entry.c
  cs_ereach.c
  cs_etree.c
  cs_fkeep.c
  cs_gaxpy.c
  cs_happly.c
  cs_house.c
  cs_ipvec.c
  cs_leaf.c
  cs_load.c
  cs_lsolve.c
  cs_ltsolve.c
  cs_lu.c
  cs_lusol.c
  cs_malloc.c
  cs_maxtrans.c
  cs_multiply.c
  cs_norm.c
  cs_permute.c
  cs_pinv.c
  cs_post.c
  cs_pvec.c
  cs_qr.c
  cs_qrsol.c
  cs_randperm.c
  cs_reach.c
  cs_scatter.c
  cs_scc.c
  cs_schol.c
  cs_spsolve.c
  cs_sqr.c
  cs_symperm.c
  cs_tdfs.c
  cs_transpose.c
  cs_updown.c
  cs_usolve.c
  cs_util.c
  cs_utsolve.c
  # the following files are not needed - they contain no symbols
  # cs_print.c
)

target_include_directories(
  cxsparse_vendored
  PRIVATE
  ${PROJECT_SOURCE_DIR}/include
  ${PROJECT_BINARY_DIR}/include
  PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}
)

if (BUILD_SHARED_LIBS)
  set_property(TARGET cxsparse_vendored PROPERTY POSITION_INDEPENDENT_CODE ON)
endif()

# Disable complex number support for CXSparse because:
#   - It is necessary to compile with MSVC
#   - igraph does not need complex number support from CXSparse on any platform
target_compile_definitions(cxsparse_vendored PUBLIC NCOMPLEX)

# Since these are included as object files, they should call the
# function as is (without a visibility specification)
target_compile_definitions(cxsparse_vendored PRIVATE IGRAPH_STATIC)

use_all_warnings(cxsparse_vendored)

if (MSVC)
  target_compile_options(cxsparse_vendored PRIVATE /wd4100) # disable unreferenced parameter warning
endif()
