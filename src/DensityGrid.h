#ifndef __DENSITY_GRID_H__
#define __DENSITY_GRID_H__


// Compile time adjustable parameters


#include <deque>

using namespace std;

#include "drl_layout.h"
#include "drl_Node.h"
#ifdef MUSE_MPI
  #include <mpi.h>
#endif

namespace drl {

class DensityGrid {

public:
  
	  // Methods
	  void Init();
	  void Subtract(Node &n, bool first_add, bool fine_first_add, bool fineDensity);
	  void Add(Node &n, bool fineDensity );
	  float GetDensity(float Nx, float Ny, bool fineDensity);

	  // Contructor/Destructor
	  DensityGrid() {};
	  ~DensityGrid();

private:

	  // Private Members
	  void Subtract( Node &N );
	  void Add( Node &N );
	  void fineSubtract( Node &N );
	  void fineAdd( Node &N );

	  // new dynamic variables -- SBM
	  float (*fall_off)[RADIUS*2+1];
	  float (*Density)[GRID_SIZE];
	  deque<Node> (*Bins)[GRID_SIZE];

	  // old static variables
	  //float fall_off[RADIUS*2+1][RADIUS*2+1];
	  //float Density[GRID_SIZE][GRID_SIZE];
	  //deque<Node *> Bins[GRID_SIZE][GRID_SIZE];
};

} // namespace drl

#endif // __DENSITY_GRID_H__

