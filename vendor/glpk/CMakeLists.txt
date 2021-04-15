
add_library(
  glpk_vendored
  OBJECT
  EXCLUDE_FROM_ALL
  glpapi01.c glpapi02.c glpapi03.c glpapi04.c glpapi05.c glpapi06.c glpapi07.c
  glpapi08.c glpapi09.c glpapi10.c glpapi11.c glpapi12.c glpapi13.c glpapi14.c
  glpapi15.c glpapi16.c glpapi17.c glpapi18.c glpapi19.c glpavl.c glpbfd.c
  glpbfx.c glpcpx.c glpdmp.c glpdmx.c glpenv01.c glpenv02.c glpenv03.c
  glpenv04.c glpenv05.c glpenv06.c glpenv07.c glpenv08.c glpfhv.c glpgmp.c
  glphbm.c glpini01.c glpini02.c glpios01.c glpios02.c glpios03.c glpios04.c
  glpios05.c glpios06.c glpios07.c glpios08.c glpios09.c glpios10.c glpios11.c
  glpios12.c glpipm.c glplib01.c glplib02.c glplib03.c glplpf.c
  glplpx01.c glplpx02.c glplpx03.c glpluf.c glpluf.h glplux.c glpmat.c
  glpmpl01.c glpmpl02.c glpmpl03.c glpmpl04.c glpmpl05.c glpmpl06.c glpmps.c
  glpnet01.c glpnet02.c glpnet03.c glpnet04.c glpnet05.c glpnet06.c glpnet07.c
  glpnet08.c glpnet09.c glpnpp.h glpnpp01.c glpnpp02.c glpnpp03.c glpnpp04.c
  glpnpp05.c glpqmd.c glpqmd.h glprgr.c glprng01.c glprng02.c glpscf.c
  glpscl.c glpsdf.c glpspm.c glpspm.h glpspx.h glpspx01.c glpspx02.c glpsql.c
  glpssx01.c glpssx02.c glpstd.h glptsp.c
  amd/amd_1.c amd/amd_2.c amd/amd_aat.c amd/amd_control.c amd/amd_defaults.c
  amd/amd_info.c amd/amd_order.c amd/amd_post_tree.c
  amd/amd_postorder.c amd/amd_preprocess.c amd/amd_valid.c
  colamd/colamd.c
  # The following files have no symbols
  amd/amd_dump.c
)

target_include_directories(
  glpk_vendored
  PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}
  PRIVATE
  ${PROJECT_SOURCE_DIR}/include
  ${PROJECT_BINARY_DIR}/include
  ${PROJECT_BINARY_DIR}/src # config.h for TLS
)

if (BUILD_SHARED_LIBS)
  set_property(TARGET glpk_vendored PROPERTY POSITION_INDEPENDENT_CODE ON)
endif()

# Since these are included as object files, they should call the
# function as is (without visibility specification)
target_compile_definitions(glpk_vendored PRIVATE IGRAPH_STATIC)

if (MSVC)
  target_compile_options(glpk_vendored PRIVATE
    /wd4068
  )
else()
  target_compile_options(glpk_vendored PRIVATE
    $<$<C_COMPILER_ID:Intel>:-wd161>
  )
endif()
