
/*
 * This file contains the graph handling routines.
 *
 * Copyright (C) 2002 Sampo Niskanen, Patric Östergård.
 * Licensed under the GNU GPL, read the file LICENSE for details.
 */


#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include "graph.h"


/*
static graph_t *graph_read_dimacs_binary(FILE *fp,char *firstline);
static graph_t *graph_read_dimacs_ascii(FILE *fp,char *firstline);
*/


/*
 * graph_new()
 *
 * Returns a newly allocated graph with n vertices all with weight 1,
 * and no edges.
 */
graph_t *graph_new(int n) {
	graph_t *g;
	int i;

	ASSERT((sizeof(setelement)*8)==ELEMENTSIZE);
	ASSERT(n>0);

	g=malloc(sizeof(graph_t));
	g->n=n;
	g->edges=malloc(g->n * sizeof(set_t));
	g->weights=malloc(g->n * sizeof(int));
	for (i=0; i < g->n; i++) {
		g->edges[i]=set_new(n);
		g->weights[i]=1;
	}
	return g;
}

/*
 * graph_free()
 *
 * Frees the memory associated with the graph g.
 */
void graph_free(graph_t *g) {
	int i;

	ASSERT((sizeof(setelement)*8)==ELEMENTSIZE);
	ASSERT(g!=NULL);
	ASSERT(g->n > 0);

	for (i=0; i < g->n; i++) {
		set_free(g->edges[i]);
	}
	free(g->weights);
	free(g->edges);
	free(g);
	return;
}


/*
 * graph_resize()
 *
 * Resizes graph g to given size.  If size > g->n, the new vertices are
 * not connected to any others and their weights are set to 1.
 * If size < g->n, the last g->n - size vertices are removed.
 */
void graph_resize(graph_t *g, int size) {
	int i;

	ASSERT(g!=NULL);
	ASSERT(g->n > 0);
	ASSERT(size > 0);

	if (g->n == size)
		return;

	/* Free/alloc extra edge-sets */
	for (i=size; i < g->n; i++)
		set_free(g->edges[i]);
	g->edges=realloc(g->edges, size * sizeof(set_t));
	for (i=g->n; i < size; i++)
		g->edges[i]=set_new(size);

	/* Resize original sets */
	for (i=0; i < MIN(g->n,size); i++) {
		g->edges[i]=set_resize(g->edges[i],size);
	}

	/* Weights */
	g->weights=realloc(g->weights,size * sizeof(int));
	for (i=g->n; i<size; i++)
		g->weights[i]=1;

	g->n=size;
	return;
}

/*
 * graph_crop()
 *
 * Resizes the graph so as to remove all highest-valued isolated vertices.
 */
void graph_crop(graph_t *g) {
	int i;

	for (i=g->n-1; i>=1; i--)
		if (set_size(g->edges[i])>0)
			break;
	graph_resize(g,i+1);
	return;
}


/*
 * graph_weighted()
 *
 * Returns TRUE if all vertex weights of graph g are all the same.
 *
 * Note: Does NOT require weights to be 1.
 */
boolean graph_weighted(graph_t *g) {
	int i,w;

	w=g->weights[0];
	for (i=1; i < g->n; i++)
		if (g->weights[i] != w)
			return TRUE;
	return FALSE;
}

/*
 * graph_edge_count()
 *
 * Returns the number of edges in graph g.
 */
int graph_edge_count(graph_t *g) {
	int i;
	int count=0;

	for (i=0; i < g->n; i++) {
		count += set_size(g->edges[i]);
	}
	return count/2;
}


#if 0
/*
 * graph_write_dimacs_ascii_file()
 *
 * Writes an ASCII dimacs-format file of graph g, with comment, to
 * given file.
 *
 * Returns TRUE if successful, FALSE if an error occurred.
 */
boolean graph_write_dimacs_ascii_file(graph_t *g, char *comment, char *file) {
	FILE *fp;

	ASSERT((sizeof(setelement)*8)==ELEMENTSIZE);
	ASSERT(file!=NULL);

	if ((fp=fopen(file,"wb"))==NULL)
		return FALSE;
	if (!graph_write_dimacs_ascii(g,comment,fp)) {
		fclose(fp);
		return FALSE;
	}
	fclose(fp);
	return TRUE;
}

/*
 * graph_write_dimacs_ascii()
 *
 * Writes an ASCII dimacs-format file of graph g, with comment, to the
 * file stream fp.
 *
 * Returns TRUE if successful, FALSE if an error occurred.
 */
boolean graph_write_dimacs_ascii(graph_t *g, char *comment, FILE *fp) {
	int i,j;

	ASSERT((sizeof(setelement)*8)==ELEMENTSIZE);
	ASSERT(graph_test(g,NULL));
	ASSERT(fp!=NULL);

	if (comment)
		fprintf(fp,"c %s\n",comment);
	fprintf(fp,"p edge %d %d\n",g->n,graph_edge_count(g));
	for (i=0; i < g->n; i++)
		if (g->weights[i]!=1)
			fprintf(fp,"n %d %d\n",i+1,g->weights[i]);
	for (i=0; i < g->n; i++)
		for (j=0; j<i; j++)
			if (GRAPH_IS_EDGE_FAST(g,i,j))
				fprintf(fp,"e %d %d\n",i+1,j+1);
	return TRUE;
}

/*
 * graph_write_dimacs_binary_file()
 *
 * Writes a binary dimacs-format file of graph g, with comment, to
 * given file.
 *
 * Returns TRUE if successful, FALSE if an error occurred.
 */
boolean graph_write_dimacs_binary_file(graph_t *g, char *comment, char *file) {
	FILE *fp;

	ASSERT((sizeof(setelement)*8)==ELEMENTSIZE);
	ASSERT(file!=NULL);

	if ((fp=fopen(file,"wb"))==NULL)
		return FALSE;
	if (!graph_write_dimacs_binary(g,comment,fp)) {
		fclose(fp);
		return FALSE;
	}
	fclose(fp);
	return TRUE;
}

/*
 * graph_write_dimacs_binary()
 *
 * Writes a binary dimacs-format file of graph g, with comment, to the
 * file stream fp.
 *
 * Returns TRUE if successful, FALSE if an error occurred.
 */

#define STR_APPEND(s) \
if (headerlength+strlen(s) >= headersize) {  \
	headersize+=1024;                    \
	header=realloc(header,headersize);   \
}                                            \
strncat(header,s,1000);                      \
headerlength+=strlen(s);

boolean graph_write_dimacs_binary(graph_t *g, char *comment,FILE *fp) {
	char *buf;
	char *header=NULL;
	int headersize=0;
	int headerlength=0;
	int i,j;

	ASSERT((sizeof(setelement)*8)==ELEMENTSIZE);
	ASSERT(graph_test(g,NULL));
	ASSERT(fp!=NULL);

	buf=malloc(MAX(1024,g->n/8+1));
	header=malloc(1024);
	header[0]=0;
	headersize=1024;
	if (comment) {
		strcpy(buf,"c ");
		strncat(buf,comment,1000);
		strcat(buf,"\n");
		STR_APPEND(buf);
	}
	sprintf(buf,"p edge %d %d\n",g->n,graph_edge_count(g));
	STR_APPEND(buf);
	for (i=0; i < g->n; i++) {
		if (g->weights[i]!=1) {
			sprintf(buf,"n %d %d\n",i+1,g->weights[i]);
			STR_APPEND(buf);
		}
	}

	fprintf(fp,"%d\n",(int)strlen(header));
	fprintf(fp,"%s",header);
	free(header);

	for (i=0; i < g->n; i++) {
		memset(buf,0,i/8+1);
		for (j=0; j<i; j++) {
			if (GRAPH_IS_EDGE_FAST(g,i,j)) {
				buf[j/8] |= SET_BIT_MASK(7-j%8);
			}
		}
		fwrite(buf,1,i/8+1,fp);
	}
	free(buf);
	return TRUE;
}



/*
 * graph_read_dimacs_file()
 *
 * Reads a dimacs-format (ASCII or binary) file from the given file.
 *
 * Returns a newly allocated graph, or NULL if an error occurred, and an
 * error message is printed to stderr.
 */
graph_t *graph_read_dimacs_file(char *file) {
	FILE *fp;
	graph_t *g;

	ASSERT((sizeof(setelement)*8)==ELEMENTSIZE);
	ASSERT(file!=NULL);

	if ((fp=fopen(file,"rb"))==NULL) {
		perror(file);
		return NULL;
	}
	g=graph_read_dimacs(fp);
	fclose(fp);
	return g;
}


/*
 * graph_read_dimacs()
 *
 * Reads a dimacs-format (ASCII or binary) file from the file stream fp.
 *
 * Returns a newly allocated graph, or NULL if an error occurred, and an
 * error message is printed to stderr.
 */
graph_t *graph_read_dimacs(FILE *fp) {
	char buffer[1024];
	graph_t *g;
	char tmp[10];
	int n;

	ASSERT((sizeof(setelement)*8)==ELEMENTSIZE);
	ASSERT(fp!=NULL);

	if (fgets(buffer,1023,fp)==NULL) {
		fprintf(stderr,"Input does not contain any data.\n");
		return NULL;
	}
	if (sscanf(buffer," %d %2s",&n,tmp)!=1) {
		g=graph_read_dimacs_ascii(fp,buffer);
	} else {
		g=graph_read_dimacs_binary(fp,buffer);
	}
	return g;
}


/*
 * parse_input()
 *
 * Parses the string str for ASCII-format dimacs commands, and modifies
 * the graph g accordingly.
 *
 * Returns TRUE if successful, FALSE if a bad command was encountered.
 *
 * Note: Ignores all unknown commands.  The 'd', 'v' and 'x' commands
 *       (mainly generator-specific information) are ignored silently,
 *       for all others a warning message is printed to stderr.
 */
static boolean parse_input(char *str,graph_t *g) {
	int i,j,w;
	char tmp[16];

	for (i=0; i<strlen(str); i++) {
		if (!isspace((int)str[i]))
			break;
	}
	if (i>=strlen(str))  /* blank line */
		return TRUE;
	if (str[i+1]!=0 && !isspace(str[i+1]))  /* not 1-char field */
		return FALSE;

	switch (str[i]) {
	case 'c':
		return TRUE;
	case 'p':
		if (g->n != 0)
			return FALSE;
		if (sscanf(str," p %15s %d %d %2s",tmp,&(g->n),&i,tmp)!=3)
			return FALSE;
		if (g->n <= 0)
			return FALSE;
		g->edges=calloc(g->n,sizeof(set_t));
		for (i=0; i<g->n; i++)
			g->edges[i]=set_new(g->n);
		g->weights=calloc(g->n,sizeof(int));
		for (i=0; i<g->n; i++)
			g->weights[i]=1;
		return TRUE;
	case 'n':
		if ((g->n <= 0) || (g->weights == NULL))
			return FALSE;
		if (sscanf(str," n %d %d %2s",&i,&w,tmp)!=2)
			return FALSE;
		if (i<1 || i>g->n)
			return FALSE;
		if (w<=0)
			return FALSE;
		g->weights[i-1]=w;
		return TRUE;
	case 'e':
		if ((g->n <= 0) || (g->edges == NULL))
			return FALSE;
		if (sscanf(str," e %d %d %2s",&i,&j,tmp)!=2)
			return FALSE;
		if (i<1 || j<1 || i>g->n || j>g->n)
			return FALSE;
		if (i==j)   /* We want antireflexive graphs. */
			return TRUE;
		GRAPH_ADD_EDGE(g,i-1,j-1);
		return TRUE;
	case 'd':
	case 'v':
	case 'x':
		return TRUE;
	default:
		fprintf(stderr,"Warning: ignoring field '%c' in "
			"input.\n",str[i]);
		return TRUE;
	}
}


/*
 * graph_read_dimacs_binary()
 *
 * Reads a dimacs-format binary file from file stream fp with the first
 * line being firstline.
 *
 * Returns the newly-allocated graph or NULL if an error occurred.
 *
 * TODO: This function leaks memory when reading erroneous files.
 */
static graph_t *graph_read_dimacs_binary(FILE *fp,char *firstline) {
	int length=0;
	graph_t *g;
	int i,j;
	char *buffer;
	char *start;
	char *end;
	char **buf;
	char tmp[10];

	if (sscanf(firstline," %d %2s",&length,tmp)!=1)
		return NULL;
	if (length<=0) {
		fprintf(stderr,"Malformed preamble: preamble size < 0.\n");
		return NULL;
	}
	buffer=malloc(length+2);
	if (fread(buffer,1,length,fp)<length) {
		fprintf(stderr,"Malformed preamble: unexpected "
			"end of file.\n");
		free(buffer);
		return NULL;
	}

	g=calloc(1,sizeof(graph_t));
	start=buffer;
	while (start < buffer+length) {
		end=strchr(start,'\n');
		if (end==NULL)
			end=buffer+length;
		end[0]=0;
		if (!parse_input(start,g)) {
			fprintf(stderr,"Malformed preamble: %s\n",start);
			free (buffer);
			return NULL;
		}
		start=end+1;
	}

	free(buffer);
	if (g->n <= 0) {
		fprintf(stderr,"Malformed preamble: number of "
			"vertices <= 0\n");
		free(g);
		return NULL;
	}

	/* Binary part. */
	buf=calloc(g->n,sizeof(char*));
	for (i=0; i < g->n; i++) {
		buf[i]=calloc(g->n,1);
		if (fread(buf[i],1,i/8+1,fp) < (i/8+1)) {
			fprintf(stderr,"Unexpected end of file when "
				"reading graph.\n");
			return NULL;
		}
	}

	for (i=0; i < g->n; i++) {
		for (j=0; j<i; j++) {
			if (buf[i][j/8]&(1<<(7-(j%8)))) {
				GRAPH_ADD_EDGE(g,i,j);
			}
		}
		free(buf[i]);
	}
	free(buf);

	return g;
}


/*
 * graph_read_dimacs_ascii()
 *
 * Reads a dimacs-format ASCII file from file stream fp with the first
 * line being firstline.
 *
 * Returns the newly-allocated graph or NULL if an error occurred.
 *
 * TODO:  This function leaks memory when reading erroneous files.
 */
static graph_t *graph_read_dimacs_ascii(FILE *fp, char *firstline) {
	graph_t *g;
	char buffer[1024];

	g=calloc(1,sizeof(graph_t));

	if (!parse_input(firstline,g)) {
		fprintf(stderr,"Malformed input: %s",firstline);
		free(g);
		return NULL;
	}
	while (fgets(buffer,1023,fp)) {
		if (!parse_input(buffer,g)) {
			fprintf(stderr,"Malformed input: %s",buffer);
			return NULL;
		}
	}
	if (g->n <= 0) {
		free(g);
		fprintf(stderr,"Unexpected end of file when reading graph.\n");
		return NULL;
	}

	return g;
}
#endif


#if 0
/*
 * graph_print()
 *
 * Prints a representation of the graph g to stdout (along with any errors
 * noticed).  Mainly useful for debugging purposes and trivial output.
 *
 * The output consists of a first line describing the dimensions and then
 * one line per vertex containing the vertex number (numbered 0,...,n-1),
 * the vertex weight (if the graph is weighted), "->" and then a list
 * of all vertices it is adjacent to.
 */
void graph_print(graph_t *g) {
	int i,j;
	int asymm=0;
	int refl=0;
	int nonpos=0;
	int extra=0;
	unsigned int weight=0;
	boolean weighted;

	ASSERT((sizeof(setelement)*8)==ELEMENTSIZE);

	if (g==NULL) {
		printf("   WARNING: Graph pointer is NULL!\n");
		return;
	}
	if (g->n <= 0) {
		printf("   WARNING: Graph has %d vertices "
		       "(should be positive)!\n",g->n);
		return;
	}

	weighted=graph_weighted(g);

	printf("%s graph has %d vertices, %d edges (density %.2f).\n",
	       weighted?"Weighted":((g->weights[0]==1)?
				    "Unweighted":"Semi-weighted"),
	       g->n,graph_edge_count(g),
	       (float)graph_edge_count(g)/((float)(g->n - 1)*(g->n)/2));

	for (i=0; i < g->n; i++) {
		printf("%2d",i);
		if (weighted) {
			printf(" w=%d",g->weights[i]);
			if (g->weights[i] <= 0) {
				printf("*NON-POSITIVE*");
				nonpos++;
			}
		}
		if (weight < INT_MAX)
			weight+=g->weights[i];
		printf(" ->");
		for (j=0; j < g->n; j++) {
			if (SET_CONTAINS_FAST(g->edges[i],j)) {
				printf(" %d",j);
				if (i==j) {
					printf("*REFLEXIVE*");
					refl++;
				}
				if (!SET_CONTAINS_FAST(g->edges[j],i)) {
					printf("*ASYMMERTIC*");
					asymm++;
				}
			}
		}
		for (j=g->n; j < SET_ARRAY_LENGTH(g->edges[i])*ELEMENTSIZE;
		     j++) {
			if (SET_CONTAINS_FAST(g->edges[i],j)) {
				printf(" %d*NON-EXISTENT*",j);
				extra++;
			}
		}
		printf("\n");
	}

	if (asymm)
		printf("   WARNING: Graph contained %d asymmetric edges!\n",
		       asymm);
	if (refl)
		printf("   WARNING: Graph contained %d reflexive edges!\n",
		       refl);
	if (nonpos)
		printf("   WARNING: Graph contained %d non-positive vertex "
		       "weights!\n",nonpos);
	if (extra)
		printf("   WARNING: Graph contained %d edges to "
		       "non-existent vertices!\n",extra);
	if (weight>=INT_MAX)
		printf("   WARNING: Total graph weight >= INT_MAX!\n");
	return;
}


/*
 * graph_test()
 *
 * Tests graph g to be valid.  Checks that g is non-NULL, the edges are
 * symmetric and anti-reflexive, and that all vertex weights are positive.
 * If output is non-NULL, prints a few lines telling the status of the graph
 * to file descriptor output.
 *
 * Returns TRUE if the graph is valid, FALSE otherwise.
 */
boolean graph_test(graph_t *g,FILE *output) {
	int i,j;
	int edges=0;
	int asymm=0;
	int nonpos=0;
	int refl=0;
	int extra=0;
	unsigned int weight=0;
	boolean weighted;

	ASSERT((sizeof(setelement)*8)==ELEMENTSIZE);

	if (g==NULL) {
		if (output)
			fprintf(output,"   WARNING: Graph pointer is NULL!\n");
		return FALSE;
	}

	weighted=graph_weighted(g);

	for (i=0; i < g->n; i++) {
		if (g->edges[i]==NULL) {
			if (output)
				fprintf(output,"   WARNING: Graph edge set "
					"NULL!\n"
					"   (further warning suppressed)\n");
			return FALSE;
		}
		if (SET_MAX_SIZE(g->edges[i]) < g->n) {
			if (output)
				fprintf(output,"   WARNING: Graph edge set "
					"too small!\n"
					"   (further warnings suppressed)\n");
			return FALSE;
		}
		for (j=0; j < g->n; j++) {
			if (SET_CONTAINS_FAST(g->edges[i],j)) {
				edges++;
				if (i==j) {
					refl++;
				}
				if (!SET_CONTAINS_FAST(g->edges[j],i)) {
					asymm++;
				}
			}
		}
		for (j=g->n; j < SET_ARRAY_LENGTH(g->edges[i])*ELEMENTSIZE;
		     j++) {
			if (SET_CONTAINS_FAST(g->edges[i],j))
				extra++;
		}
		if (g->weights[i] <= 0)
			nonpos++;
		if (weight<INT_MAX)
			weight += g->weights[i];
	}

	edges/=2;  /* Each is counted twice. */

	if (output) {
		/* Semi-weighted means all weights are equal, but not 1. */
		fprintf(output,"%s graph has %d vertices, %d edges "
			"(density %.2f).\n",
			weighted?"Weighted":
			((g->weights[0]==1)?"Unweighted":"Semi-weighted"),
			g->n,edges,(float)edges/((float)(g->n - 1)*(g->n)/2));

		if (asymm)
			fprintf(output,"   WARNING: Graph contained %d "
				"asymmetric edges!\n",asymm);
		if (refl)
			fprintf(output,"   WARNING: Graph contained %d "
				"reflexive edges!\n",refl);
		if (nonpos)
			fprintf(output,"   WARNING: Graph contained %d "
				"non-positive vertex weights!\n",nonpos);
		if (extra)
			fprintf(output,"   WARNING: Graph contained %d edges "
				"to non-existent vertices!\n",extra);
		if (weight>=INT_MAX)
			fprintf(output,"   WARNING: Total graph weight >= "
				"INT_MAX!\n");
		if (asymm==0 && refl==0 && nonpos==0 && extra==0 &&
		    weight<INT_MAX)
			fprintf(output,"Graph OK.\n");
	}

	if (asymm || refl || nonpos || extra || weight>=INT_MAX)
		return FALSE;

	return TRUE;
}


/*
 * graph_test_regular()
 *
 * Returns the vertex degree for regular graphs, or -1 if the graph is
 * not regular.
 */
int graph_test_regular(graph_t *g) {
	int i,n;

	n=set_size(g->edges[0]);

	for (i=1; i < g->n; i++) {
		if (set_size(g->edges[i]) != n)
			return -1;
	}
	return n;
}

#endif
