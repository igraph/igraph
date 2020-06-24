find_package(PYTHON)

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
set(DOXROX_OUTPUTS "")

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
if(Python_FOUND)
  set(DOXROX_COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/doxrox.py)
  set(DOXROX_RULES ${CMAKE_CURRENT_SOURCE_DIR}/c-docbook.re)
  set(DOXROX_CHUNKS ${CMAKE_CURRENT_BINARY_DIR}/chunks.pickle)

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

  add_custom_target(
    doc
    DEPENDS
    ${DOXROX_OUTPUTS}
    COMMENT
    "Built DocBook XML files for the documentation"
  )
endif(Python_FOUND)
