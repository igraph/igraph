# check_tls_support(TLS_KEYWORD)
set(TLS_KEYWORD "")  # when using a third-party ARPACK, TLS support cannot be guaranteed
if(TLS_KEYWORD)
  set(HAVE_TLS 1)
else()
  set(HAVE_TLS 0)
endif()

