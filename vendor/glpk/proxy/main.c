/* Last update: 08-May-2013 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "glpk.h"
#include "proxy.h"

/**********************************************************************/
int main(int argc, char **argv)
/**********************************************************************/
{
    glp_prob *lp;
    int ncols, status;
    double *initsol, zstar, *xstar;

    /* check arguments */
    if ( (argc == 1) || (argc > 3) ) {
        printf("ERROR: Usage: ts <instance> <(possibly) xml initsols>\n"
              );
        exit(1);
    }

    /* creating the problem */
    lp = glp_create_prob();
    glp_set_prob_name(lp, "Proxy");

    /* reading the problem */
    glp_term_out(GLP_OFF);
#if 0 /* by mao */
    status = glp_read_lp(lp, NULL, argv[1]);
#else
    status = glp_read_mps(lp, GLP_MPS_FILE, NULL, argv[1]);
#endif
    glp_term_out(GLP_ON);
    if ( status ) {
        printf("Problem %s does not exist!!!, status %d\n",
               argv[1], status);
        exit(1);
    }

    ncols = glp_get_num_cols(lp);

    initsol = (double *) calloc(ncols+1, sizeof(double));

    if (argc == 3) {
        FILE *fp=fopen(argv[2],"r");
        char  tmp[256]={0x0};
        int counter = 1;
        while(fp!=NULL && fgets(tmp, sizeof(tmp),fp)!=NULL)
        {
            char *valini = strstr(tmp, "value");
            if (valini!=NULL){
                int num;
                double dnum;
                valini +=7;
                sscanf(valini, "%d%*s",&num);
                dnum = (double)num;
                initsol[counter] = dnum;
                counter++;
            }
        }
        fclose(fp);
    }

    xstar = (double *) calloc(ncols+1, sizeof(double));

    if (argc == 3) {
        status = proxy(lp, &zstar, xstar, initsol, 0.0, 0, 1);
    }
    else {
        status = proxy(lp, &zstar, xstar, NULL, 0.0, 0, 1);
    }

    printf("Status = %d; ZSTAR = %f\n",status,zstar);
    /*
    int i;
    for (i=1; i< ncols+1; i++) {
        printf("XSTAR[%d] = %f\n",i, xstar[i]);
    }
     */

    glp_delete_prob(lp);

    return 0;
}
