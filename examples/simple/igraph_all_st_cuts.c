/* -*- mode: C -*-  */
/* 
   IGraph library.
   Copyright (C) 2010  Gabor Csardi <csardi.gabor@gmail.com>
   Rue de l'Industrie 5, Lausanne 1005, Switzerland
   
   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
   
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
   
   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 
   02110-1301 USA

*/

#include <igraph.h>
#include <igraph_marked_queue.h>

int igraph_i_all_st_cuts_pivot(const igraph_t *graph,
			       const igraph_marked_queue_t *S,
			       const igraph_stack_t *T,
			       const igraph_vector_bool_t *TV,
			       long int source,
			       long int target,
			       long int *v,
			       igraph_vector_t *Isv);

int main() {
  igraph_t g;
  igraph_vector_ptr_t cuts, partition1s;
  long int i, n;

  igraph_marked_queue_t S;
  igraph_stack_t T;
  igraph_vector_bool_t TV;
  long int v;
  igraph_vector_t Isv;

  /* ----------------------------------------------------------- */
  /* This is the example from the Provan-Shier paper, 
     for calculating the dominator tree and finding the right pivot 
     element */
  
  igraph_small(&g, 12, IGRAPH_DIRECTED,
  	       /* a->b */ 0,1,
  	       /* b->t */ 1,11,
  	       /* c->b */ 2,1,  /* c->d */ 2,3,
  	       /* d->e */ 3,4,  /* d->i */ 3,8,
  	       /* e->c */ 4,2,
  	       /* f->c */ 5,2,  /* f->e */ 5,4,
  	       /* g->d */ 6,3,  /* g->e */ 6,4,  /* g->f */ 6,5,
  	                        /* g->j */ 6,9,
  	       /* h->g */ 7,6,  /* h->t */ 7,11,
  	       /* i->a */ 8,0,
  	       /* j->i */ 9,8,
  	       /* s->a */ 10,0, /* s->c */ 10,2, /* s->h */ 10,7,
  	       -1);
  
  /* S={s,a} */
  igraph_marked_queue_init(&S, igraph_vcount(&g));
  igraph_marked_queue_start_batch(&S);
  igraph_marked_queue_push(&S, 10);
  igraph_marked_queue_push(&S, 0);
  
  /* T={t} */
  igraph_stack_init(&T, 1);
  igraph_stack_push(&T, 11);
  igraph_vector_bool_init(&TV, igraph_vcount(&g));
  VECTOR(TV)[11] = 1;

  igraph_vector_init(&Isv, 0);
  igraph_i_all_st_cuts_pivot(&g, &S, &T, &TV,
  				/*source=*/ 10, /*target=*/ 11,
  				&v, &Isv);

  /* Expected result: v=c, Isv={c,d,e,i} */
  printf("%li; ", v);
  igraph_vector_print(&Isv);
  
  igraph_vector_destroy(&Isv);
  igraph_vector_bool_destroy(&TV);
  igraph_stack_destroy(&T);
  igraph_marked_queue_destroy(&S);
  igraph_destroy(&g);

  /* ----------------------------------------------------------- */

  igraph_small(&g, 3, IGRAPH_DIRECTED,
  	       0,1, 1,2,
  	       -1);
  
  /* S={}, T={} */
  igraph_marked_queue_init(&S, igraph_vcount(&g));
  igraph_stack_init(&T, 3);
  igraph_vector_bool_init(&TV, igraph_vcount(&g));

  igraph_vector_init(&Isv, 0);
  igraph_i_all_st_cuts_pivot(&g, &S, &T, &TV,
  				/*source=*/ 0, /*target=*/ 2,
  				&v, &Isv);
  printf("%li; ", v);
  igraph_vector_print(&Isv);

  igraph_vector_destroy(&Isv);
  igraph_vector_bool_destroy(&TV);
  igraph_stack_destroy(&T);
  igraph_marked_queue_destroy(&S);
  igraph_destroy(&g);

  /* ----------------------------------------------------------- */

  igraph_small(&g, 3, IGRAPH_DIRECTED,
  	       0,1, 1,2,
  	       -1);
  
  /* S={}, T={0} */
  igraph_marked_queue_init(&S, igraph_vcount(&g));

  igraph_stack_init(&T, 3);
  igraph_vector_bool_init(&TV, igraph_vcount(&g));
  igraph_stack_push(&T, 0);
  VECTOR(TV)[0]=1;

  igraph_vector_init(&Isv, 0);
  igraph_i_all_st_cuts_pivot(&g, &S, &T, &TV,
  				/*source=*/ 0, /*target=*/ 2,
  				&v, &Isv);
  printf("%li; ", v);
  igraph_vector_print(&Isv);

  igraph_vector_destroy(&Isv);
  igraph_vector_bool_destroy(&TV);
  igraph_stack_destroy(&T);
  igraph_marked_queue_destroy(&S);
  igraph_destroy(&g);

  /* ----------------------------------------------------------- */

  igraph_small(&g, 3, IGRAPH_DIRECTED,
  	       0,1, 1,2,
  	       -1);
  
  /* S={0}, T={} */
  igraph_marked_queue_init(&S, igraph_vcount(&g));
  igraph_marked_queue_push(&S, 0);

  igraph_stack_init(&T, 3);
  igraph_vector_bool_init(&TV, igraph_vcount(&g));

  igraph_vector_init(&Isv, 0);
  igraph_i_all_st_cuts_pivot(&g, &S, &T, &TV,
  				/*source=*/ 0, /*target=*/ 2,
  				&v, &Isv);
  printf("%li; ", v);
  igraph_vector_print(&Isv);

  igraph_vector_destroy(&Isv);
  igraph_vector_bool_destroy(&TV);
  igraph_stack_destroy(&T);
  igraph_marked_queue_destroy(&S);
  igraph_destroy(&g);

  /* ----------------------------------------------------------- */

  igraph_small(&g, 3, IGRAPH_DIRECTED,
  	       0,1, 1,2,
  	       -1);
  
  /* S={0}, T={1} */
  igraph_marked_queue_init(&S, igraph_vcount(&g));
  igraph_marked_queue_push(&S, 0);

  igraph_stack_init(&T, 3);
  igraph_vector_bool_init(&TV, igraph_vcount(&g));
  igraph_stack_push(&T, 1);
  VECTOR(TV)[1] = 1;

  igraph_vector_init(&Isv, 0);
  igraph_i_all_st_cuts_pivot(&g, &S, &T, &TV,
  				/*source=*/ 0, /*target=*/ 2,
  				&v, &Isv);
  printf("%li; ", v);
  igraph_vector_print(&Isv);

  igraph_vector_destroy(&Isv);
  igraph_vector_bool_destroy(&TV);
  igraph_stack_destroy(&T);
  igraph_marked_queue_destroy(&S);
  igraph_destroy(&g);

  /* ----------------------------------------------------------- */

  igraph_small(&g, 3, IGRAPH_DIRECTED,
  	       0,1, 1,2,
  	       -1);
  
  /* S={0,1}, T={} */
  igraph_marked_queue_init(&S, igraph_vcount(&g));
  igraph_marked_queue_push(&S, 0);
  igraph_marked_queue_push(&S, 1);

  igraph_stack_init(&T, 3);
  igraph_vector_bool_init(&TV, igraph_vcount(&g));

  igraph_vector_init(&Isv, 0);
  igraph_i_all_st_cuts_pivot(&g, &S, &T, &TV,
  				/*source=*/ 0, /*target=*/ 2,
  				&v, &Isv);
  printf("%li; ", v);
  igraph_vector_print(&Isv);

  igraph_vector_destroy(&Isv);
  igraph_vector_bool_destroy(&TV);
  igraph_stack_destroy(&T);
  igraph_marked_queue_destroy(&S);
  igraph_destroy(&g);

  /* ----------------------------------------------------------- */

  igraph_small(&g, 3, IGRAPH_DIRECTED,
  	       0,1, 1,2,
  	       -1);

  igraph_vector_ptr_init(&cuts, 0);
  igraph_vector_ptr_init(&partition1s, 0);
  igraph_all_st_cuts(&g, /*cuts=*/ 0, &partition1s,
		     /*source=*/ 0, /*target=*/ 2);

  n=igraph_vector_ptr_size(&partition1s);
  for (i=0; i<n; i++) {
    igraph_vector_t *v=VECTOR(partition1s)[i];
    igraph_vector_print(v);
  }
  
  igraph_destroy(&g);

  /* ----------------------------------------------------------- */

  igraph_small(&g, 5, IGRAPH_DIRECTED,
  	       0,1, 1,2, 1,3, 2,4, 3,4,
  	       -1);

  igraph_vector_ptr_init(&cuts, 0);
  igraph_vector_ptr_init(&partition1s, 0);
  igraph_all_st_cuts(&g, /*cuts=*/ 0, &partition1s,
		     /*source=*/ 0, /*target=*/ 4);

  n=igraph_vector_ptr_size(&partition1s);
  for (i=0; i<n; i++) {
    igraph_vector_t *v=VECTOR(partition1s)[i];
    igraph_vector_print(v);
  }
  
  igraph_destroy(&g);  

  /* ----------------------------------------------------------- */

  igraph_small(&g, 6, IGRAPH_DIRECTED,
  	       0,1, 1,2, 1,3, 2,4, 3,4, 1,5, 5,4,
  	       -1);

  igraph_vector_ptr_init(&cuts, 0);
  igraph_vector_ptr_init(&partition1s, 0);
  igraph_all_st_cuts(&g, &cuts, &partition1s,
		     /*source=*/ 0, /*target=*/ 4);

  n=igraph_vector_ptr_size(&partition1s);
  printf("Partitions and cuts:\n");
  for (i=0; i<n; i++) {
    igraph_vector_t *v=VECTOR(partition1s)[i];
    igraph_vector_t *v2=VECTOR(cuts)[i];
    printf("P: ");
    igraph_vector_print(v);
    printf("C: ");
    igraph_vector_print(v2);
  }
  
  igraph_destroy(&g);  

  return 0;
}
