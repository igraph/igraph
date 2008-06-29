// This file contains the member definitions of the DensityGrid.h class
// This code is modified from the original code by B.N. Wylie

#include <string>
#include <deque>
#include <iostream>
#include <cmath>
#include <cstdlib>

using namespace std;

#include "drl_Node.h"
#include "DensityGrid.h"

namespace drl {

//*******************************************************
// Density Grid Destructor -- deallocates memory used
// for Density matrix, fall_off matrix, and node deque.

DensityGrid::~DensityGrid ()
{
	delete[] Density;
	delete[] fall_off;
	delete[] Bins;
}

/*********************************************
* Function: Density_Grid::Reset		         *
* Description: Reset the density grid		 *
*********************************************/
// changed from reset to init since we will only
// call this once in the parallel version of layout

void DensityGrid::Init() 
{
  
  try
    {
      Density = new float[GRID_SIZE][GRID_SIZE];
      fall_off = new float[RADIUS*2+1][RADIUS*2+1];
      Bins = new deque<Node>[GRID_SIZE][GRID_SIZE];
    }
  catch (bad_alloc errora)
    {
      cout << "Error: Out of memory! Program stopped." << endl;
	  #ifdef MUSE_MPI
        MPI_Abort ( MPI_COMM_WORLD, 1 );
	  #else
	    exit (1);
	  #endif
    }
	
  // Clear Grid
  int i;
  for (i=0; i< GRID_SIZE; i++) 
    for (int j=0; j< GRID_SIZE; j++) {
      Density[i][j] = 0;
      Bins[i][j].erase(Bins[i][j].begin(),Bins[i][j].end());
    }
  
  // Compute fall off
  for(i=-RADIUS; i<=RADIUS; i++)
    for(int j=-RADIUS; j<=RADIUS; j++) {
      fall_off[i+RADIUS][j+RADIUS] = (float)((RADIUS-fabs((float)i))/RADIUS) * 
	(float)((RADIUS-fabs((float)j))/RADIUS);
    }
	
}


/***************************************************
 * Function: DensityGrid::GetDensity               *
 * Description: Get_Density from density grid      *
 **************************************************/
float DensityGrid::GetDensity(float Nx, float Ny, bool fineDensity) 
{
	deque<Node>::iterator BI;
	int x_grid, y_grid;
	float x_dist, y_dist, distance, density=0;
	int boundary=10;	// boundary around plane


	/* Where to look */
	x_grid = (int)((Nx+HALF_VIEW+.5)*VIEW_TO_GRID);
	y_grid = (int)((Ny+HALF_VIEW+.5)*VIEW_TO_GRID);

	// Check for edges of density grid (10000 is arbitrary high density)
	if (x_grid > GRID_SIZE-boundary || x_grid < boundary) return 10000;
	if (y_grid > GRID_SIZE-boundary || y_grid < boundary) return 10000;

	// Fine density?
	if (fineDensity) {

		// Go through nearest bins
		for(int i=y_grid-1; i<=y_grid+1; i++)
			for(int j=x_grid-1; j<=x_grid+1; j++) {

			// Look through bin and add fine repulsions
			for(BI = Bins[i][j].begin(); BI < Bins[i][j].end(); ++BI) {
				x_dist =  Nx-(BI->x);
				y_dist =  Ny-(BI->y);
				distance = x_dist*x_dist+y_dist*y_dist;
				density += 1e-4/(distance + 1e-50);
		 }
		}

	// Course density
	} else {

		// Add rough estimate
		density = Density[y_grid][x_grid];
		density *= density;
	}

	return density;
}

/// Wrapper functions for the Add and subtract methods
/// Nodes should all be passed by constant ref

void DensityGrid::Add(Node &n, bool fineDensity)
{
  if(fineDensity)
    fineAdd(n);
  else
    Add(n);
}

void DensityGrid::Subtract( Node &n, bool first_add,
							bool fine_first_add, bool fineDensity)
{
  if ( fineDensity && !fine_first_add ) fineSubtract (n);
  else if ( !first_add ) Subtract(n);
}

			
/***************************************************
 * Function: DensityGrid::Subtract                *
 * Description: Subtract a node from density grid  *
 **************************************************/
void DensityGrid::Subtract(Node &N) 
{
  int x_grid, y_grid, diam;
  float *den_ptr, *fall_ptr;
	
  /* Where to subtract */
  x_grid = (int)((N.sub_x+HALF_VIEW+.5)*VIEW_TO_GRID);
  y_grid = (int)((N.sub_y+HALF_VIEW+.5)*VIEW_TO_GRID);
  x_grid -= RADIUS;
  y_grid -= RADIUS;
  diam = 2*RADIUS;

  /* Subtract density values */
  den_ptr = &Density[y_grid][x_grid];
  fall_ptr = &fall_off[0][0];
  for(int i = 0; i <= diam; i++) {
    for(int j = 0; j <= diam; j++)
	 *den_ptr++ -= *fall_ptr++;
    den_ptr += GRID_SIZE - (diam+1);
  }
}

/***************************************************
 * Function: DensityGrid::Add                     *
 * Description: Add a node to the density grid     *
 **************************************************/
void DensityGrid::Add(Node &N) 
{

  int x_grid, y_grid, diam;
  float *den_ptr, *fall_ptr;


  /* Where to add */
  x_grid = (int)((N.x+HALF_VIEW+.5)*VIEW_TO_GRID);
  y_grid = (int)((N.y+HALF_VIEW+.5)*VIEW_TO_GRID);
 
  N.sub_x = N.x;
  N.sub_y = N.y;
  
  x_grid -= RADIUS;
  y_grid -= RADIUS;
  diam = 2*RADIUS;

  // check to see that we are inside grid
  if ( (x_grid >= GRID_SIZE) || (x_grid < 0) ||
       (y_grid >= GRID_SIZE) || (y_grid < 0) )
    {
      cout << endl << "Error: Exceeded density grid with x_grid = " << x_grid 
	       << " and y_grid = " << y_grid << ".  Program stopped." << endl;
      #ifdef MUSE_MPI
 	    MPI_Abort ( MPI_COMM_WORLD, 1 );
	  #else
	    exit (1);
	  #endif
    }    

  /* Add density values */
  den_ptr = &Density[y_grid][x_grid];
  fall_ptr = &fall_off[0][0];
  for(int i = 0; i <= diam; i++) {
    for(int j = 0; j <= diam; j++)
	 *den_ptr++ += *fall_ptr++;
    den_ptr += GRID_SIZE - (diam+1);
  }
  
}

/***************************************************
 * Function: DensityGrid::fineSubtract             *
 * Description: Subtract a node from bins		   *
 **************************************************/
void DensityGrid::fineSubtract(Node &N) 
{
  int x_grid, y_grid;

  /* Where to subtract */
  x_grid = (int)((N.sub_x+HALF_VIEW+.5)*VIEW_TO_GRID);
  y_grid = (int)((N.sub_y+HALF_VIEW+.5)*VIEW_TO_GRID);
  Bins[y_grid][x_grid].pop_front();
}

/***************************************************
 * Function: DensityGrid::fineAdd                  *
 * Description: Add a node to the bins			   *
 **************************************************/
void DensityGrid::fineAdd(Node &N) 
{
  int x_grid, y_grid;

  /* Where to add */
  x_grid = (int)((N.x+HALF_VIEW+.5)*VIEW_TO_GRID);
  y_grid = (int)((N.y+HALF_VIEW+.5)*VIEW_TO_GRID);
  N.sub_x = N.x;
  N.sub_y = N.y;
  Bins[y_grid][x_grid].push_back(N);
}

} // namespace drl
