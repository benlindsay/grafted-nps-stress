#define MAIN
#include "globals.h"
#include "r2.h"
double simulate(void);

// The main routine is essentially a wrapper for the simulation routine
int main(int argc, char** argv) {

#ifdef PAR
  // MPI initialization stuff
  MPI_Init( &argc , &argv );
  MPI_Comm_rank( MPI_COMM_WORLD, &myrank);
  MPI_Comm_size( MPI_COMM_WORLD, &nprocs);
  fftw_mpi_init();
#else
  myrank = 0;
#endif

  // Other initialization stuff
  idum = -long(time(0)) * (myrank + 1);
  read_input();
  initialize_1();
  allocate();

  // Single simulation if not doing Brent's method
  simulate();

#ifdef PAR
  MPI_Finalize();
#endif

  return 0;
}
