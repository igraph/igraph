# Specify the list of .xml files that are used as-is
set(
  DOCBOOK_SOURCES
  fdl.xml
  gpl.xml
  igraph-docs.xml
  installation.xml
  introduction.xml
  licenses.xml
  pmt.xml
  tutorial.xml
)

# Specify the list of .xxml files that have to be piped through doxrox to
# obtain the final set of .xml files that serve as an input to DocBook
set(
  DOXROX_SOURCES
  adjlist.xxml
  arpack.xxml
  attributes.xxml
  basicigraph.xxml
  bipartite.xxml
  cliques.xxml
  coloring.xxml
  community.xxml
  dqueue.xxml
  embedding.xxml
  error.xxml
  flows.xxml
  foreign.xxml
  generators.xxml
  graphlets.xxml
  heap.xxml
  hrg.xxml
  isomorphism.xxml
  iterators.xxml
  layout.xxml
  matrix.xxml
  memory.xxml
  motifs.xxml
  nongraph.xxml
  operators.xxml
  progress.xxml
  random.xxml
  scg.xxml
  separators.xxml
  sparsemat.xxml
  sparsematrix.xxml
  spatialgames.xxml
  stack.xxml
  status.xxml
  structural.xxml
  strvector.xxml
  threading.xxml
  vector.xxml
  visitors.xxml
)

# Specify the igraph source files that may contain documentation chunks
file(
  GLOB IGRAPH_SOURCES_FOR_DOXROX
  LIST_DIRECTORIES FALSE
  ${CMAKE_SOURCE_DIR}/include/*.h
  ${CMAKE_SOURCE_DIR}/include/*.h.in
  ${CMAKE_SOURCE_DIR}/src/*.c
  ${CMAKE_SOURCE_DIR}/src/*.cc
  ${CMAKE_SOURCE_DIR}/src/*.cpp
  ${CMAKE_SOURCE_DIR}/src/*.h
  ${CMAKE_SOURCE_DIR}/src/*.pmt
)

# Specify the igraph source files that are used as examples in the
# documentation
file(
  GLOB DOCBOOK_EXAMPLES
  LIST_DIRECTORIES FALSE
  RELATIVE ${CMAKE_SOURCE_DIR}
  ${CMAKE_SOURCE_DIR}/examples/simple/*.c
  ${CMAKE_SOURCE_DIR}/examples/tests/*.c
)

# You should not need to change anything below this line if you are simply
# trying to add new files to produce documentation from

# Documentation build requires Python and DocBook
find_package(PYTHON)
find_program(SOURCE_HIGHLIGHT source-highlight)
find_program(XMLTO xmlto)

set(DOC_BUILD_SUPPORTED Python_FOUND AND XMLTO AND SOURCE_HIGHLIGHT)

if(DOC_BUILD_SUPPORTED)
  set(DOXROX_COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/doxrox.py)
  set(DOXROX_RULES ${CMAKE_CURRENT_SOURCE_DIR}/c-docbook.re)
  set(DOXROX_CHUNKS ${CMAKE_CURRENT_BINARY_DIR}/chunks.pickle)

  set(HTML_STAMP ${CMAKE_CURRENT_BINARY_DIR}/html/stamp)

  set(DOCBOOK_INPUTS "")
  set(DOXROX_OUTPUTS "")

  # Specify that each DocBook .xml file is to be copied to the build folder
  # TODO(ntamas): currently this works with out-of-tree builds only
  set(IGRAPH_VERSION ${PACKAGE_VERSION})   # for replacement in igraph-docs.xml
  foreach(DOCBOOK_SOURCE ${DOCBOOK_SOURCES})
    set(DOCBOOK_INPUT "${CMAKE_CURRENT_BINARY_DIR}/${DOCBOOK_SOURCE}")
    list(APPEND DOCBOOK_INPUTS "${DOCBOOK_INPUT}")
    configure_file(${DOCBOOK_SOURCE} ${DOCBOOK_INPUT})
  endforeach()

  # Specify that each example source file is to be piped through source-higlight
  # to produce an .xml representation that can be used in the DocBook
  # documentation
  foreach(DOCBOOK_EXAMPLE_SOURCE ${DOCBOOK_EXAMPLES})
    string(REGEX REPLACE "[.]c$" ".c.xml" DOCBOOK_EXAMPLE_OUTPUT ${DOCBOOK_EXAMPLE_SOURCE})
    set(COMMENT "Highlighting source code in ${DOCBOOK_EXAMPLE_SOURCE}")

    set(DOCBOOK_EXAMPLE_OUTPUT "${CMAKE_BINARY_DIR}/${DOCBOOK_EXAMPLE_SOURCE}.xml")
    list(APPEND DOCBOOK_INPUTS "${DOCBOOK_EXAMPLE_OUTPUT}")

    get_filename_component(DOCBOOK_EXAMPLE_OUTPUT_DIR "${DOCBOOK_EXAMPLE_OUTPUT}" DIRECTORY)

    add_custom_command(
      OUTPUT ${DOCBOOK_EXAMPLE_OUTPUT}
      COMMAND ${CMAKE_COMMAND} -E make_directory ${DOCBOOK_EXAMPLE_OUTPUT_DIR}
      COMMAND ${SOURCE_HIGHLIGHT}
      ARGS
      --src-lang c
      --out-format docbook
      --input ${CMAKE_SOURCE_DIR}/${DOCBOOK_EXAMPLE_SOURCE}
      --output ${DOCBOOK_EXAMPLE_OUTPUT}
      # --gen-references inline
      # --ctags=""
      --outlang-def ${CMAKE_SOURCE_DIR}/doc/docbook.outlang
      MAIN_DEPENDENCY ${CMAKE_SOURCE_DIR}/${DOCBOOK_EXAMPLE_SOURCE}
      COMMENT ${COMMENT}
    )
  endforeach()

  foreach(DOXROX_SOURCE ${DOXROX_SOURCES})
    string(REGEX REPLACE "[.]xxml$" ".xml" DOXROX_OUTPUT ${DOXROX_SOURCE})
    set(COMMENT "Generating ${DOXROX_OUTPUT} from ${DOXROX_SOURCE}")

    string(PREPEND DOXROX_OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/")
    list(APPEND DOXROX_OUTPUTS "${DOXROX_OUTPUT}")

    add_custom_command(
      OUTPUT ${DOXROX_OUTPUT}
      COMMAND ${DOXROX_COMMAND}
      ARGS
      -t ${CMAKE_CURRENT_SOURCE_DIR}/${DOXROX_SOURCE}
      --chunks ${DOXROX_CHUNKS}
      -o ${DOXROX_OUTPUT}
      MAIN_DEPENDENCY ${CMAKE_CURRENT_SOURCE_DIR}/${DOXROX_SOURCE}
      DEPENDS ${DOXROX_CHUNKS}
      COMMENT ${COMMENT}
    )
  endforeach()

  add_custom_command(
    OUTPUT ${DOXROX_CHUNKS}
    COMMAND ${DOXROX_COMMAND}
    ARGS
    -e ${DOXROX_RULES}
    -o ${DOXROX_CHUNKS}
    ${IGRAPH_SOURCES_FOR_DOXROX}
    MAIN_DEPENDENCY ${DOXROX_RULES}
    DEPENDS ${IGRAPH_SOURCES_FOR_DOXROX}
    COMMENT "Parsing documentation chunks from source code"
  )

  add_custom_command(
    OUTPUT ${HTML_STAMP}
    COMMAND ${XMLTO} -x ${CMAKE_CURRENT_SOURCE_DIR}/gtk-doc.xsl -o html xhtml igraph-docs.xml
    COMMAND ${CMAKE_COMMAND} -E touch ${HTML_STAMP}
    MAIN_DEPENDENCY igraph-docs.xml
    DEPENDS ${DOXROX_OUTPUTS} ${DOCBOOK_INPUTS}
    COMMENT "Generating HTML documentation with xmlto"
  )

  add_custom_target(
    doc
    DEPENDS ${HTML_STAMP}
    COMMENT
    "Built DocBook XML files for the documentation"
  )
endif()
