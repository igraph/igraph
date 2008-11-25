/*
 *  SCGlib : A C library for the spectral coarse graining of matrices
 *	as described in the paper: Shrinking Matrices while preserving their
 *	eigenpairs with Application to the Spectral Coarse Graining of Graphs.
 *	Preprint available at <http://people.epfl.ch/david.morton>
 *  
 *	Copyright (C) 2008 David Morton de Lachapelle <david.morton@a3.epfl.ch>
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 
 *  02110-1301 USA
 *
 *  DESCRIPTION
 *	-----------
 *    The kmeans_Lloyd function is adapted from the R-stats package.
 *	  It perfoms Lloyd's k-means clustering on a p x n data matrix
 *	  stored row-wise in a vector 'x'. 'cen' contains k initial centers.
 *	  The group label to which each object belongs is stored in 'cl'.
 *	  Labels are positive consecutive integers starting from 0.
 *	  See also Section 5.3.3 of the above reference.
 */
 
#include "scg_headers.h"

INT kmeans_Lloyd(const REAL *x, const UINT n, const UINT p, REAL *cen,
						const UINT k, INT *cl, const UINT maxiter)
{
    UINT iter, i, j, c, it, inew = 0;
    REAL best, dd, tmp;
    UINT updated;
	UINT *nc = uint_vector(k);

    for(i = 0; i < n; i++) cl[i] = -1;
    for(iter = 0; iter < maxiter; iter++) {
		updated = 0;
		for(i = 0; i < n; i++) {
	//find nearest centre for each point
			best = LDBL_MAX;
			for(j = 0; j < k; j++) {
				dd = 0.0;
				for(c = 0; c < p; c++) {
					tmp = x[i+n*c] - cen[j+k*c];
					dd += tmp * tmp;
				}
				if(dd < best) {
					best = dd;
					inew = j+1;
				}
			}
			if(cl[i] != inew) {
				updated = 1;
				cl[i] = inew;
			}
		}
		if(!updated) break;
	//update each centre
		for(j = 0; j < k*p; j++) cen[j] = 0.0;
		for(j = 0; j < k; j++) nc[j] = 0;
		for(i = 0; i < n; i++) {
			it = cl[i] - 1;
			nc[it]++;
			for(c = 0; c < p; c++) cen[it+c*k] += x[i+c*n];
			}
		for(j = 0; j < k*p; j++) cen[j] /= nc[j % k];
	}
	free_uint_vector(nc);
	//returns 1 if converged else -1
	if(iter<maxiter-1)
		return 1;
	else
		return -1;
}

/*another kmeans algorithm that might one day be useful (needs refactoring)*/
/*void kmeans_MacQueen(double *x, int *pn, int *pp, double *cen, int *pk,
		     int *cl, int *pmaxiter, int *nc, double *wss)
{
    int n = *pn, k = *pk, p = *pp, maxiter = *pmaxiter;
    int iter, i, j, c, it, inew = 0, iold;
    double best, dd, tmp;
    Rboolean updated;

    /* first assign each point to the nearest cluster centre 
    for(i = 0; i < n; i++) {
	best = R_PosInf;
	for(j = 0; j < k; j++) {
	    dd = 0.0;
	    for(c = 0; c < p; c++) {
		tmp = x[i+n*c] - cen[j+k*c];
		dd += tmp * tmp;
	    }
	    if(dd < best) {
		best = dd;
		inew = j+1;
	    }
	}
	if(cl[i] != inew) cl[i] = inew;
    }
   /* and recompute centres as centroids 
    for(j = 0; j < k*p; j++) cen[j] = 0.0;
    for(j = 0; j < k; j++) nc[j] = 0;
    for(i = 0; i < n; i++) {
	it = cl[i] - 1; nc[it]++;
	for(c = 0; c < p; c++) cen[it+c*k] += x[i+c*n];
    }
    for(j = 0; j < k*p; j++) cen[j] /= nc[j % k];

    for(iter = 0; iter < maxiter; iter++) {
	updated = FALSE;
	for(i = 0; i < n; i++) {
	    best = R_PosInf;
	    for(j = 0; j < k; j++) {
		dd = 0.0;
		for(c = 0; c < p; c++) {
		    tmp = x[i+n*c] - cen[j+k*c];
		    dd += tmp * tmp;
		}
		if(dd < best) {
		    best = dd;
		    inew = j;
		}
	    }
	    if((iold = cl[i] - 1) != inew) {
		updated = TRUE;
		cl[i] = inew + 1;
		nc[iold]--; nc[inew]++;
		/* update old and new cluster centres 
		for(c = 0; c < p; c++) {
		    cen[iold+k*c] += (cen[iold+k*c] - x[i+n*c])/nc[iold];
		    cen[inew+k*c] += (x[i+n*c] - cen[inew+k*c])/nc[inew];
		}
	    }
	}
	if(!updated) break;
    }

    *pmaxiter = iter + 1;
    for(j = 0; j < k; j++) wss[j] = 0.0;
    for(i = 0; i < n; i++) {
	it = cl[i] - 1;
	for(c = 0; c < p; c++) {
	    tmp = x[i+n*c] - cen[it+k*c];
	    wss[it] += tmp * tmp;
	}
    }
}*/
