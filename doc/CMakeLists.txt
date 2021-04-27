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
  cycles.xxml
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
  psumtree.xxml
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
  GLOB_RECURSE IGRAPH_SOURCES_FOR_DOXROX
  LIST_DIRECTORIES FALSE
  ${CMAKE_SOURCE_DIR}/include/*.h
  ${CMAKE_SOURCE_DIR}/include/*.h.in
  ${CMAKE_BINARY_DIR}/include/*.h
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
  ${CMAKE_SOURCE_DIR}/examples/tutorial/*.c
)

# You should not need to change anything below this line if you are simply
# trying to add new files to produce documentation from

# Documentation build requires Python and source-highlight
find_package(Python3)
find_program(SOURCE_HIGHLIGHT_COMMAND source-highlight)

# HTML documentation additionally requires xmlto from DocBook
find_program(XMLTO_COMMAND xmlto)

# PDF documentation additionally requires xsltproc, xmllint and Apache FOP
find_program(FOP_COMMAND fop)
find_program(XMLLINT_COMMAND xmllint)
find_program(XSLTPROC_COMMAND xsltproc)

if(Python3_FOUND AND SOURCE_HIGHLIGHT_COMMAND)
  set(DOC_BUILD_SUPPORTED TRUE)
else()
  set(DOC_BUILD_SUPPORTED FALSE)
endif()

if(DOC_BUILD_SUPPORTED AND XMLTO_COMMAND)
  set(HTML_DOC_BUILD_SUPPORTED TRUE)
else()
  set(HTML_DOC_BUILD_SUPPORTED FALSE)
endif()

if(DOC_BUILD_SUPPORTED AND XMLLINT_COMMAND AND XSLTPROC_COMMAND AND FOP_COMMAND)
  set(PDF_DOC_BUILD_SUPPORTED TRUE)
else()
  set(PDF_DOC_BUILD_SUPPORTED FALSE)
endif()

if(DOC_BUILD_SUPPORTED)
  set(DOXROX_COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/doxrox.py)
  set(DOXROX_RULES ${CMAKE_CURRENT_SOURCE_DIR}/c-docbook.re)
  set(DOXROX_CHUNKS ${CMAKE_CURRENT_BINARY_DIR}/chunks.pickle)
  set(DOXROX_CACHE ${CMAKE_CURRENT_BINARY_DIR}/doxrox.cache)

  set(DOCBOOK_INPUTS "")
  set(DOCBOOK_GENERATED_INPUTS "")

  # Specify that each DocBook .xml file is to be copied to the build folder
  # TODO(ntamas): currently this works with out-of-tree builds only
  set(IGRAPH_VERSION ${PACKAGE_VERSION})   # for replacement in igraph-docs.xml
  foreach(DOCBOOK_SOURCE ${DOCBOOK_SOURCES})
    set(DOCBOOK_INPUT "${CMAKE_CURRENT_BINARY_DIR}/${DOCBOOK_SOURCE}")
    list(APPEND DOCBOOK_INPUTS "${DOCBOOK_INPUT}")
    configure_file(${DOCBOOK_SOURCE} ${DOCBOOK_INPUT})
  endforeach()

  # Specify that .xxml files should be piped through doxrox.py to get a
  # DocBook-compatible .xml file. This step inserts the documentation chunks
  # extracted from the igraph source to the DocBook sources
  foreach(DOXROX_SOURCE ${DOXROX_SOURCES})
    string(REGEX REPLACE "[.]xxml$" ".xml" DOXROX_OUTPUT ${DOXROX_SOURCE})
    set(COMMENT "Generating ${DOXROX_OUTPUT} from ${DOXROX_SOURCE}")

    string(PREPEND DOXROX_OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/")
    list(APPEND DOCBOOK_INPUTS "${DOXROX_OUTPUT}")
    list(APPEND DOCBOOK_GENERATED_INPUTS "${DOXROX_OUTPUT}")

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

  # When all .xxml and .xml files have been processed, we have to send them
  # through a custom Python script that extracts the ID references and produces
  # a ctags-compatible "tags" file. This will then be used later by
  # source-highlight to cross-reference the known tokens from the source code
  # of the examples
  list(JOIN DOCBOOK_GENERATED_INPUTS ";" DOCBOOK_GENERATED_INPUTS_AS_STRING)
  add_custom_command(
    OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/tags"
    COMMAND ${CMAKE_COMMAND}
    ARGS
    -DINPUT_FILES="${DOCBOOK_GENERATED_INPUTS_AS_STRING}"
    -DOUTPUT_FILE=${CMAKE_CURRENT_BINARY_DIR}/tags
    -P ${CMAKE_SOURCE_DIR}/etc/cmake/generate_tags_file.cmake
    DEPENDS ${DOCBOOK_GENERATED_INPUTS}
    COMMENT "Creating tags file from DocBook xmls"
  )

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
      COMMAND ${SOURCE_HIGHLIGHT_COMMAND}
      ARGS
      --src-lang c
      --out-format docbook
      --input ${CMAKE_SOURCE_DIR}/${DOCBOOK_EXAMPLE_SOURCE}
      --output ${DOCBOOK_EXAMPLE_OUTPUT}
      --gen-references inline
      --ctags=""
      --outlang-def ${CMAKE_SOURCE_DIR}/doc/docbook.outlang
      MAIN_DEPENDENCY ${CMAKE_SOURCE_DIR}/${DOCBOOK_EXAMPLE_SOURCE}
      DEPENDS tags
      COMMENT ${COMMENT}
    )
  endforeach()

  add_custom_command(
    OUTPUT ${DOXROX_CHUNKS} ${DOXROX_CACHE}
    COMMAND ${DOXROX_COMMAND}
    ARGS
    -e ${DOXROX_RULES}
    -o ${DOXROX_CHUNKS}
    --cache ${DOXROX_CACHE}
    ${IGRAPH_SOURCES_FOR_DOXROX}
    MAIN_DEPENDENCY ${DOXROX_RULES}
    DEPENDS ${IGRAPH_SOURCES_FOR_DOXROX}
    COMMENT "Parsing documentation chunks from source code"
  )

  if(HTML_DOC_BUILD_SUPPORTED)
    set(HTML_STAMP ${CMAKE_CURRENT_BINARY_DIR}/html/stamp)

    add_custom_command(
      OUTPUT ${HTML_STAMP}
      COMMAND ${XMLTO_COMMAND} -x ${CMAKE_CURRENT_SOURCE_DIR}/gtk-doc.xsl -o html xhtml igraph-docs.xml
      COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/html/*.css ${CMAKE_CURRENT_BINARY_DIR}/html
      COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/html/*.js ${CMAKE_CURRENT_BINARY_DIR}/html
      COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/html/*.png ${CMAKE_CURRENT_BINARY_DIR}/html
      COMMAND ${CMAKE_COMMAND} -E touch ${HTML_STAMP}
      MAIN_DEPENDENCY igraph-docs.xml
      DEPENDS ${DOCBOOK_INPUTS}
      COMMENT "Generating HTML documentation with xmlto"
    )

    add_custom_target(html DEPENDS ${HTML_STAMP})
    set(HTML_TARGET html)
  endif()

  if(PDF_DOC_BUILD_SUPPORTED)
    add_custom_command(
      OUTPUT igraph-docs-with-resolved-includes.xml
      COMMAND ${XMLLINT_COMMAND}
      ARGS
      --xinclude
      --output igraph-docs-with-resolved-includes-tmp.xml
      igraph-docs.xml
      COMMAND ${Python3_EXECUTABLE}
      ARGS
      ${CMAKE_SOURCE_DIR}/tools/removeexamples.py
      igraph-docs-with-resolved-includes-tmp.xml
      igraph-docs-with-resolved-includes.xml
      COMMAND ${CMAKE_COMMAND}
      ARGS
      -E remove igraph-docs-with-resolved-includes-tmp.xml
      MAIN_DEPENDENCY igraph-docs.xml
      DEPENDS ${DOCBOOK_INPUTS}
      COMMENT "Resolving includes in DocBook XML source"
    )

    add_custom_command(
      OUTPUT igraph-docs.fo
      COMMAND ${XSLTPROC_COMMAND}
      ARGS
      --output igraph-docs.fo
      --stringparam paper.type A4
      http://docbook.sourceforge.net/release/xsl/current/fo/docbook.xsl
      igraph-docs-with-resolved-includes.xml
      MAIN_DEPENDENCY igraph-docs-with-resolved-includes.xml
      COMMENT "Converting DocBook XML to Apache FOP format"
    )

    add_custom_command(
      OUTPUT igraph-docs.pdf
      COMMAND ${FOP_COMMAND}
      ARGS -fo igraph-docs.fo -pdf igraph-docs.pdf
      MAIN_DEPENDENCY igraph-docs.fo
      COMMENT "Generating PDF documentation with Apache FOP"
    )

    add_custom_target(pdf DEPENDS igraph-docs.pdf)
    set(PDF_TARGET pdf)
  endif()

  add_custom_target(doc DEPENDS ${HTML_TARGET} ${PDF_TARGET})
endif()

set(HTML_DOC_BUILD_SUPPORTED ${HTML_DOC_BUILD_SUPPORTED} PARENT_SCOPE)
set(PDF_DOC_BUILD_SUPPORTED ${PDF_DOC_BUILD_SUPPORTED} PARENT_SCOPE)
