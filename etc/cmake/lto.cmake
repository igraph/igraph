include(helpers)

tristate(IGRAPH_ENABLE_LTO "Enable link-time optimization" OFF)

include(CheckIPOSupported)

if(IGRAPH_ENABLE_LTO)
  # this matches both ON and AUTO
  check_ipo_supported(RESULT IPO_SUPPORTED OUTPUT IPO_NOT_SUPPORTED_REASON)
  if(IGRAPH_ENABLE_LTO STREQUAL "AUTO")
    # autodetection
    set(IGRAPH_ENABLE_LTO ${IPO_SUPPORTED})
  elseif(IPO_SUPPORTED)
    # user wanted LTO and the compiler supports it
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
  else()
    # user wanted LTO and the compiler does not support it
    message(FATAL_ERROR "Link-time optimization not supported on this compiler")
  endif()
endif()
