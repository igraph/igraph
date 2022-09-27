
#include <igraph.h>

#include <locale.h>
#include <stdio.h>
#include <string.h>

int main(void) {
    const char *filename = "weighted.gml";
    igraph_t graph;
    igraph_safelocale_t loc;

    /* Attempt to set a locale that uses a decimal comma. Locale names
     * differ between platforms, and not all locales are available,
     * so the locale change may not be successful. */
    const char *locname = setlocale(LC_ALL, "de_DE");
    struct lconv *lc = localeconv();
    if (strcmp(lc->decimal_point, ",")) {
        /* If decimal point is not a comma, presumably because the requested
         * locale was not available, report locale information. */
        fprintf(stderr, "setlocale() returned '%s', decimal point is '%s'\n",
                locname ? locname : "NULL",
                lc->decimal_point);
    }

    FILE *file = fopen(filename, "r");
    if (! file) {
        fprintf(stderr, "Cannot open %s file.\n", filename);
        exit(1);
    }

    /* An attribute table is needed to read graph attributes. */
    igraph_set_attribute_table(&igraph_cattribute_table);

    /* At this point, the current locale may use decimal commas.
     * We temporarily set a C locale using enter_safelocale() to
     * allow the GML reader and writer to work correctly.*/
    igraph_enter_safelocale(&loc);
    if (igraph_read_graph_gml(&graph, file) != IGRAPH_SUCCESS) {
        fprintf(stderr, "Reading %s failed.\n", filename);
        abort();
    }
    igraph_write_graph_gml(&graph, stdout, IGRAPH_WRITE_GML_DEFAULT_SW, NULL, "");
    igraph_exit_safelocale(&loc);

    igraph_destroy(&graph);

    return 0;
}
