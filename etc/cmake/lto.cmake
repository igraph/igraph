include(helpers)

tristate(IGRAPH_ENABLE_LTO "Enable link-time optimization" AUTO)

include(CheckIPOSupported)

if(IGRAPH_ENABLE_LTO)
  # this matches both ON and AUTO
  check_ipo_supported(RESULT IPO_SUPPORTED OUTPUT IPO_NOT_SUPPORTED_REASON)
  if(IGRAPH_ENABLE_LTO STREQUAL "AUTO")
    set(IGRAPH_ENABLE_LTO ${IPO_SUPPORTED})
  endif()

  if(IPO_SUPPORTED)
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
  else()
    message(FATAL_ERROR "Link-time optimization not supported on this compiler")
  endif()
endif()
