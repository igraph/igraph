
add_library(
  glpk_vendored
  OBJECT
  EXCLUDE_FROM_ALL

  amd/amd_1.c amd/amd_2.c amd/amd_aat.c amd/amd_control.c amd/amd_defaults.c
  amd/amd_info.c amd/amd_order.c amd/amd_post_tree.c
  amd/amd_postorder.c amd/amd_preprocess.c amd/amd_valid.c

  api/advbas.c api/asnhall.c api/asnlp.c api/asnokalg.c api/ckasn.c api/ckcnf.c
  api/cplex.c api/cpp.c api/cpxbas.c api/graph.c api/gridgen.c api/intfeas1.c
  api/maxffalg.c api/maxflp.c api/mcflp.c api/mcfokalg.c api/mcfrelax.c
  api/minisat1.c api/mpl.c api/mps.c api/netgen.c api/npp.c api/pript.c
  api/prmip.c api/prob1.c api/prob2.c api/prob3.c api/prob4.c api/prob5.c
  api/prrngs.c api/prsol.c api/rdasn.c api/rdcc.c api/rdcnf.c api/rdipt.c
  api/rdmaxf.c api/rdmcf.c api/rdmip.c api/rdprob.c api/rdsol.c api/rmfgen.c
  api/strong.c api/topsort.c api/wcliqex.c api/weak.c api/wrasn.c api/wrcc.c
  api/wrcnf.c api/wript.c api/wrmaxf.c api/wrmcf.c api/wrmip.c api/wrprob.c
  api/wrsol.c

  bflib/btf.c bflib/btfint.c bflib/fhv.c bflib/fhvint.c bflib/ifu.c bflib/luf.c
  bflib/lufint.c bflib/scf.c bflib/scfint.c bflib/sgf.c bflib/sva.c

  colamd/colamd.c

  draft/bfd.c draft/bfx.c draft/glpapi06.c draft/glpapi07.c draft/glpapi08.c
  draft/glpapi09.c draft/glpapi10.c draft/glpapi12.c draft/glpapi13.c
  draft/glpios01.c draft/glpios02.c draft/glpios03.c draft/glpios07.c
  draft/glpios09.c draft/glpios11.c draft/glpios12.c draft/glpipm.c
  draft/glpmat.c draft/glpscl.c draft/glpssx01.c draft/glpssx02.c draft/lux.c

  env/alloc.c env/dlsup.c env/env.c env/error.c env/stdc.c env/stdout.c
  env/stream.c env/time.c env/tls.c

  intopt/cfg.c intopt/cfg1.c intopt/cfg2.c intopt/clqcut.c intopt/covgen.c
  intopt/fpump.c intopt/gmicut.c intopt/gmigen.c intopt/mirgen.c intopt/spv.c

  minisat/minisat.c

  misc/avl.c misc/bignum.c misc/dimacs.c misc/dmp.c misc/ffalg.c misc/fp2rat.c
  misc/fvs.c misc/gcd.c misc/hbm.c misc/jd.c misc/keller.c misc/ks.c
  misc/mc13d.c misc/mc21a.c misc/mt1.c misc/mygmp.c misc/okalg.c misc/qmd.c
  misc/relax4.c misc/rgr.c misc/rng.c misc/rng1.c misc/round2n.c misc/spm.c
  misc/str2int.c misc/str2num.c misc/strspx.c misc/strtrim.c misc/triang.c
  misc/wclique.c misc/wclique1.c

  mpl/mpl1.c mpl/mpl2.c mpl/mpl3.c mpl/mpl4.c mpl/mpl5.c mpl/mpl6.c
  mpl/mplsql.c

  npp/npp1.c npp/npp2.c npp/npp3.c npp/npp4.c npp/npp5.c npp/npp6.c

  proxy/proxy.c proxy/proxy1.c

  simplex/spxat.c simplex/spxchuzc.c simplex/spxchuzr.c simplex/spxlp.c
  simplex/spxnt.c simplex/spxprim.c simplex/spxprob.c simplex/spychuzc.c
  simplex/spychuzr.c simplex/spydual.c

  # amd/amd_dump.c has no symbols
)

target_include_directories(
  glpk_vendored
  PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}
  PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/amd
  ${CMAKE_CURRENT_SOURCE_DIR}/api
  ${CMAKE_CURRENT_SOURCE_DIR}/bflib
  ${CMAKE_CURRENT_SOURCE_DIR}/colamd
  ${CMAKE_CURRENT_SOURCE_DIR}/draft
  ${CMAKE_CURRENT_SOURCE_DIR}/env
  ${CMAKE_CURRENT_SOURCE_DIR}/intopt
  ${CMAKE_CURRENT_SOURCE_DIR}/minisat
  ${CMAKE_CURRENT_SOURCE_DIR}/misc
  ${CMAKE_CURRENT_SOURCE_DIR}/mpl
  ${CMAKE_CURRENT_SOURCE_DIR}/npp
  ${CMAKE_CURRENT_SOURCE_DIR}/simplex
  ${PROJECT_SOURCE_DIR}/include
  ${PROJECT_BINARY_DIR}/include
  ${PROJECT_BINARY_DIR}/src # config.h for TLS
)

# We are using IGRAPH_FILE_BASENAME in glpk/env/env.h
define_file_basename_for_sources(glpk_vendored)

if (BUILD_SHARED_LIBS)
  set_property(TARGET glpk_vendored PROPERTY POSITION_INDEPENDENT_CODE ON)
endif()

# Since these are included as object files, they should call the
# function as is (without visibility specification)
target_compile_definitions(glpk_vendored PRIVATE IGRAPH_STATIC)

# GLPK requires __WOE__ to be defined when building for Windows,
# either with MSVC or with MinGW.
# See w64/config_VC in the original GLPK distribution
if (WIN32)
  target_compile_definitions(glpk_vendored PRIVATE __WOE__=1)
endif()

if (MSVC)
  target_compile_options(glpk_vendored PRIVATE
    /wd4068
  )
else()
  target_compile_options(glpk_vendored PRIVATE
    $<$<C_COMPILER_ID:Intel>:-wd161 -Wno-return-type>
    $<$<C_COMPILER_ID:GCC,Clang,AppleClang,IntelLLVM>:-Wno-return-type -Wno-unused-value -Wno-dangling-else -Wno-logical-op-parentheses>
  )
endif()
