
#include <igraph.h>
#include <stdio.h>

#define FILENAME "mybool.graphml.xml"

int main() {
  igraph_t graph;

  FILE* graphFile = fopen("cattr_bool_bug2.graphml", "r");
  igraph_i_set_attribute_table(&igraph_cattribute_table);

  if (!graphFile) {
    printf("Cannot open input file");
    return 1;
  }

  igraph_read_graph_graphml(&graph, graphFile, 0);

  fclose(graphFile);

  if (igraph_cattribute_has_attr(&graph, IGRAPH_ATTRIBUTE_GRAPH, "mybool")) {
      igraph_bool_t value = igraph_cattribute_GAB(&graph, "mybool");
      printf("found boolean value %d\n", value);
  } else {
      return 2;
  }

  igraph_destroy(&graph);
  
  return 0;
}
