option(IGRAPH_ENABLE_LTO "Enable link-time optimization" OFF)

include(CheckIPOSupported)

if(IGRAPH_ENABLE_LTO)
  check_ipo_supported(RESULT IPO_SUPPORTED OUTPUT IPO_NOT_SUPPORTED_REASON)
  if(IPO_SUPPORTED)
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
  else()
    message(FATAL_ERROR "Link-time optimization not supported on this compiler")
  endif()
endif()

