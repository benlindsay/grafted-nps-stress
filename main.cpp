#define MAIN
#include "globals.h"
#include "r2.h"
void write_outputs(void);
void accumulate_all_averages(void);
void update_1s(void);
void update_Euler(void);
void write_data(char*, complex<double>*);
double brent_method(double, double, double);
void change_L0(double);
double simulate(void);

// The main routine is essentially a wrapper that allocates memory, then either
// runs a single simuation or implements Brent's method to find the optimum box
// size before calculating the structure and energy at that box size.
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

  // Flag that's only marked true for first simulation init
  first_sim = 1;

  if (do_brent) {
    // Brent's method to optimize box size
    double L_brent = brent_method(L_low, L_high, brent_tol);
    if (myrank == 0) {
      printf("L_brent = %lf\n\n", L_brent);
    }
    // Simulate at L_brent if that wasn't just done
    if (L[0] != L_brent) {
      change_L0(L_brent);
      simulate();
    }
    if (myrank == 0) printf("L_ideal = %lf\n\n", L_ideal);
    // Run a final simulation at the global optimum L
    if (L[0] != L_ideal) {
      change_L0(L_ideal);
      simulate();
    }
  }
  else
    // Single simulation if not doing Brent's method
    simulate();

#ifdef PAR
  MPI_Finalize();
#endif

  return 0;
}
