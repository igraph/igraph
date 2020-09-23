file(GET_RUNTIME_DEPENDENCIES
  LIBRARIES ${lib}
  RESOLVED_DEPENDENCIES_VAR found_deps
  UNRESOLVED_DEPENDENCIES_VAR unfound_deps
  POST_EXCLUDE_REGEXES WINDOWS)

configure_file(${lib} ${dest_dir} COPYONLY)

# Copy each dependency to the destination
foreach (found_dep IN LISTS found_deps)
  configure_file(${found_dep} ${dest_dir} COPYONLY)
endforeach ()