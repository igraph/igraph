tristate(IGRAPH_ENABLE_TLS "Enable thread-local storage for igraph global variables" AUTO)

include(CheckTLSSupport)
check_tls_support(TLS_KEYWORD)

if(IGRAPH_ENABLE_TLS STREQUAL "AUTO")
  if(TLS_KEYWORD)
    set(IGRAPH_ENABLE_TLS ON)
  else()
    set(IGRAPH_ENABLE_TLS OFF)
  endif()
endif()

if(IGRAPH_ENABLE_TLS)
  if(NOT TLS_KEYWORD)
    message(FATAL_ERROR "Thread-local storage not supported on this compiler")
  endif()

  # TODO: we should probably set this only if we are building igraph with
  # internal-everything
  set(IGRAPH_THREAD_SAFE YES)
else()
  set(TLS_KEYWORD "")
  set(IGRAPH_THREAD_SAFE NO)
endif()
