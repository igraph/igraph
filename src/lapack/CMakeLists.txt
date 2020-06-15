# Declare the files needed to compile our vendored BLAS copy
add_library(
  blas_vendored
  STATIC
  EXCLUDE_FROM_ALL
  dscal.c dswap.c lsame.c dnrm2.c daxpy.c dgemv.c dger.c dgemm.c
  dcopy.c dtrmm.c dtrmv.c drot.c ddot.c dasum.c dsymv.c dsyr2k.c dsyr2.c
  dtrsm.c dsyrk.c dtrsv.c idamax.c xerbla.c len_trim.c
)

target_link_libraries(
  blas_vendored
  PUBLIC
  f2c_vendored
)

# Declare the files needed to compile our vendored LAPACK copy
add_library(
  lapack_vendored
  STATIC
  EXCLUDE_FROM_ALL
  dgeev.c dgebak.c dgebal.c disnan.c dlaisnan.c dgehrd.c dgehd2.c
  dlarf.c iladlc.c iladlr.c dlarfg.c dlapy2.c dlahr2.c dlacpy.c dlarfb.c
  ilaenv.c ieeeck.c iparmq.c dhseqr.c dlahqr.c dlabad.c dlanv2.c dlaqr0.c
  dlaqr3.c dlaqr4.c dlaqr2.c dlaset.c dormhr.c dormqr.c dlarft.c dorm2r.c
  dtrexc.c dlaexc.c dlange.c dlassq.c dlarfx.c dlartg.c dlasy2.c dlaqr5.c
  dlaqr1.c dlascl.c dorghr.c dorgqr.c dorg2r.c dtrevc.c dlaln2.c dladiv.c
  dsyevr.c dlansy.c dormtr.c dormql.c dorm2l.c dstebz.c dlaebz.c dstein.c
  dlagtf.c dlagts.c dlarnv.c dlaruv.c dstemr.c dlae2.c dlaev2.c dlanst.c
  dlarrc.c dlarre.c dlarra.c dlarrb.c dlaneg.c dlarrd.c dlarrk.c dlasq2.c
  dlasq3.c dlasq4.c dlasq5.c dlasq6.c dlasrt.c dlarrj.c dlarrr.c dlarrv.c
  dlar1v.c dlarrf.c dsterf.c dsytrd.c dlatrd.c dsytd2.c dlanhs.c dgeqr2.c
  dtrsen.c dlacn2.c dtrsyl.c dlasr.c dsteqr.c dgeevx.c dtrsna.c dlaqtr.c
  dgetrf.c dgetf2.c dlaswp.c dgetrs.c dgesv.c dpotrf.c dpotf2.c
  dlamch.c ../fortran_intrinsics.c
)

target_link_libraries(
  lapack_vendored
  PUBLIC
  blas_vendored
  f2c_vendored
)

# Declare the files needed to compile our vendored ARPACK copy
add_library(
  arpack_vendored
  STATIC
  EXCLUDE_FROM_ALL
  dnaupd.c dnaup2.c dgetv0.c dvout.c second.c dmout.c dnaitr.c ivout.c dnapps.c
  dnconv.c dneigh.c dlaqrb.c dngets.c dsortc.c dstatn.c dneupd.c dsaupd.c
  dsaup2.c dsaitr.c dsapps.c dsconv.c dseigt.c dstqrb.c dsgets.c dsortr.c
  dstats.c dseupd.c dsesrt.c
)

target_link_libraries(
  arpack_vendored
  PUBLIC
  f2c_vendored
)

# Suppress some warnings that occur in the output because we do not want to
# mess around with the source of f2c too much to fix these
if(NOT MSVC)
  target_compile_options(blas_vendored PRIVATE -Wno-logical-op-parentheses)
  target_compile_options(lapack_vendored PRIVATE -Wno-logical-op-parentheses -Wno-shift-op-parentheses)
endif()
