# Declare the files needed to compile our vendored CXSparse copy
add_library(
  cxsparse_vendored
  STATIC
  cs_add.c       cs_happly.c    cs_pvec.c
  cs_amd.c       cs_house.c     cs_qr.c
  cs_chol.c      cs_ipvec.c     cs_qrsol.c
  cs_cholsol.c   cs_leaf.c      cs_randperm.c
  cs_compress.c  cs_load.c      cs_reach.c
  cs_counts.c    cs_lsolve.c    cs_scatter.c
  cs_cumsum.c    cs_ltsolve.c   cs_scc.c
  cs_dfs.c       cs_lu.c        cs_schol.c
  cs_dmperm.c    cs_lusol.c     cs_spsolve.c
  cs_droptol.c   cs_malloc.c    cs_sqr.c
  cs_dropzeros.c cs_maxtrans.c  cs_symperm.c
  cs_dupl.c      cs_multiply.c  cs_tdfs.c
  cs_entry.c     cs_norm.c      cs_transpose.c
  cs_ereach.c    cs_permute.c   cs_updown.c
  cs_etree.c     cs_pinv.c      cs_usolve.c
  cs_fkeep.c     cs_post.c      cs_util.c
  cs_gaxpy.c     cs_utsolve.c
  # the following files are not needed - they contain no symbols
  # cs_print.c
)

target_include_directories(
  cxsparse_vendored
  PRIVATE
  ${CMAKE_SOURCE_DIR}/include
)

# Suppress some warnings that occur in the output because we do not want to
# mess around with the source of f2c too much to fix these
use_all_warnings(cxsparse_vendored)
