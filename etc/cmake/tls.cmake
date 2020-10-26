option(IGRAPH_ENABLE_TLS "Enable thread-local storage for igraph global variables" OFF)

if(IGRAPH_ENABLE_TLS)
  include(CheckTLSSupport)
  check_tls_support(TLS_KEYWORD)

  if(NOT TLS_KEYWORD)
    message(FATAL_ERROR "Thread-local storage not supported on this compiler")
  endif()
else()
  set(TLS_KEYWORD "")
endif()
