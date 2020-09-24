file(GET_RUNTIME_DEPENDENCIES
  LIBRARIES ${lib}
  RESOLVED_DEPENDENCIES_VAR found_deps
  UNRESOLVED_DEPENDENCIES_VAR unfound_deps
  PRE_EXCLUDE_REGEXES "^(api|ext)-ms-(win|onecore).*.dll"
  PRE_EXCLUDE_REGEXES "^vcruntime.*1d.dll"
  POST_EXCLUDE_REGEXES ^[A-Z]:[\\/][Ww][Ii][Nn][Dd][Oo][Ww][Ss]
  )

if (unfound_deps)
  message(WARNING "Could not find dependencies: ${unfound_deps}")
endif()

# Copy library itself
configure_file(${lib} ${dest_dir} COPYONLY)

# Copy each dependency to the destination
foreach (found_dep IN LISTS found_deps)
  configure_file(${found_dep} ${dest_dir} COPYONLY)
endforeach ()