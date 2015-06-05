#define MAIN
#include "globals.h"
#include "r2.h"
void numerical_F_test( int , int ) ;
complex<double> compute_dumb_term(void);
void write_outputs( void );
void accumulate_all_averages(void );
void remove_particle(int);
void add_particle(int);
void update_1s( void ) ;
void update_Euler( void ) ;
double get_r( int , double* ) ;
void numerical_calc_P( void ) ;
void write_data( char*, complex<double>* );
double brent_method(double, double, double);
double get_slope(double);
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
#endif

  // Other initialization stuff
  idum = -long(time(0)) * (myrank + 1); // -25*myrank ;//
  read_input();
  initialize_1();
  allocate();

  // Flag that's only marked true for first simulation init
  first_sim = 1;

  if (do_brent) {
    // Brent's method to optimize box size
    double L_ideal = brent_method(L_low, L_high, brent_tol);
    if (myrank == 0)
      cout << "L_ideal=" << L_ideal << endl;
  }
  else
    // Single simulation if not doing Brent's method
    simulate();

#ifdef PAR
  MPI_Finalize();
#endif

  return 0;
}

// This implementation of Brent's method optimizes box size in the x-direction
// (L[0]) and maintains the same x-y ratio. The z dimension doesn't change.
double brent_method(double lowerLimit, double upperLimit, double errorTol) {
  double a = lowerLimit;
  double b = upperLimit;
  double c = 0;
  double d = upperLimit * 1.25;

  // Open file that will store final results for each simulation run.
  // One line is added to this file at the end of each pass through the
  // simulate method
  brent_otp = fopen("brent.dat", "w");

  double fa = get_slope(a);
  double fb = get_slope(b);

  double fc = 0;
  double s = 0;
  double fs = 0;

  // if f(a) f(b) >= 0 then error-exit
  if (fa * fb >= 0 ) {
    if (myrank == 0) cout << "Bad length range for Brent's method" << endl;
    if (fa < fb)
      return a;
    else
      return b;
  }

  // if |f(a)| < |f(b)| then swap (a,b) end if
  if (abs(fa) < abs(fb)) {
    double tmp = a; a = b; b = tmp; tmp = fa; fa = fb; fb = tmp;
  }

  c = a;
  fc = fa;
  bool mflag = true;

  // Iterate max 5 times so this doesn't go on forever
  for (int j=0; j<5; j++) {
    if ((fa != fc) && (fb != fc)) {
      // Inverse quadratic interpolation
      s = a*fb*fc / (fa-fb) / (fa-fc) +
          b*fa*fc / (fb-fa) / (fb-fc) +
          c*fa*fb / (fc-fa) / (fc-fb);
    }
    else
      // Secant method (linear interpolation)
      s = b - fb * (b-a) / (fb-fa);

    // Use bisection method if any of the following conditions is true:
    // Condition 1: s is not between (3a+b)/4 and b
    double tmp2 = (3 * a + b) / 4;
    bool cond1 = (s<tmp2 && s<b) || (s>tmp2 && s>b);
    // Condition 2: mflag is set and |s-b|>=|b-c|/2
    bool cond2 = mflag && (abs(s-b) >= abs(b-c)/2);
    // Condition 3: mflag is cleared and |s-b|>=|c-d|/2
    bool cond3 = !mflag && (abs(s-b) >= abs(c-d)/2);
    // Condition 4: mflag is set and |b-c|<|delta|
    bool cond4 = mflag && (abs(b-c) < errorTol);
    // Condition 5: mflag is cleared and |c-d|<|delta|
    bool cond5 = !mflag && (abs(c-d) < errorTol);

    if (cond1 || cond2 || cond3 || cond4 || cond5) {
      s = (a + b) / 2;
      mflag = true;
    }
    else
      mflag = false;

    fs = get_slope(s);
    d = c;
    c = b;
    fc = fb;
    if (fa * fs < 0) { b = s; fb = fs; }
    else { a = s; fa = fs; }

    // if |f(a)| < |f(b)| then swap (a,b) end if
    if (abs(fa) < abs(fb)) {
      double tmp = a; a = b; b = tmp; tmp = fa; fa = fb; fb = tmp;
    }

    // exit loop if f(s)=0 or |b-a| is small enough
    if (fs==0.0 || abs(b-a) < errorTol) break;
  }

  // Close the brent.dat file
  fclose(brent_otp);

  return s;
}

// Returns the slope of the H/V vs L graph at the L-value input to get_slope.
// This is a part of Brent's method where H/V is calculated at 2 nearby L 
// values by running separate simulations, then the slope is calculated
double get_slope(double l) {
  // Change L[0] to l-L_step/2 and keep same L[1] / L[0] ratio for 3D
  if (Dim==3) L[1] *= (l-0.5*L_step) / L[0];
  L[0] = l - 0.5 * L_step;
  double H_over_V_1 = simulate();

  // Add L_step to L[0] and keep same L[1] / L[0] ratio for 3D
  if (Dim==3) L[1] += L_step * L[1] / L[0];
  L[0] += L_step;
  double H_over_V_2 = simulate();

  // Return slope of H/V vs L graph
  return (H_over_V_2 - H_over_V_1) / L_step;
}

// This is the routine that is essentially the main routine in a code that
// doesn't use Brent's method. The equilibrium structure and hamiltonian are
// calculated. The real part of the hamiltonian is returned.
double simulate() {
  complex<double> Hcur, Ho, H;
  double error;
  FILE *otp;
  otp = fopen("data.dat", "w");

  if (myrank == 0)
    printf("Starting simulation for L[0]=%lf, L[1]=%lf\n", L[0], L[1]);

  // Initialize variables and fields
  initialize_1();
  initialize_2();
 
#ifdef PAR
  MPI_Barrier( MPI_COMM_WORLD );
#endif
  calc_poly_density();

  if (myrank == 0) {
    printf("Initial densities calculated!\n");
  }

  printf("Segment counts:\n");
  cout << "nD * N = " << nD * N << " integ(rhoda + rhodb) = "
    << integ_trapPBC( rhoda ) + integ_trapPBC( rhodb ) << endl; 
  cout << "nAH * Nah = " << nAH * Nah << " integ(rhoha) = "
    << integ_trapPBC( rhoha ) << endl;
  fflush(stdout);

  Ho = calc_H();
  if (myrank == 0) 
    cout << "Starting H: " << Ho << endl;

  write_outputs();

  if (myrank == 0) { 
    printf("Entering main loop!\n\n"); 
  }

  ///////////////
  // MAIN LOOP //
  ///////////////

  for (iter; iter<=itermax; iter++) {
    if (update_scheme == 0)
      update_Euler();
    else
      update_1s();

    if (do_CL && iter >= sample_wait && iter % sample_freq == 0) 
      accumulate_all_averages();

    ////////////
    // OUTPUT //
    ////////////
    if (iter % print_freq == 0) {

      Ho = Hcur;
      H = Hcur = calc_H();

      if (Hcur != Hcur) {
        printf("Crashed! iteration %d\n", iter);
        printf("H: %lf Qd: %lf\n", real(H), real(Qd));
        write_data_bin("crashed.wpl", wpl);
        write_data_bin("crashed.wabp", wabp);
        write_data_bin("crashed.wabm", wabm);
        exit(1);
      }

      if (myrank == 0) {
        printf( "Iteration: %d, H=%lf -log(Qd): %lf -log(Qha): %lf ",
                iter, real(H), real(-log(Qd)), real(-log(Qha)) );
        if (do_CL)
          printf(" + i%lf", imag(H));
        printf("\n");
      }

      error = abs( H - Ho ) / V ;

      if ( iter > 0 && myrank == 0 ) {
        fprintf(otp,"%d %5.6lf %1.3e  %5.6lf %5.6lf " , 
            iter , real(H) , imag(H) , real(-log(Qd)) , real(-log(Qha) ) ) ;

        fprintf(otp , "\n");
        fflush(otp);
      }
      
      write_outputs();

    }// output


    if (!do_CL && iter > 25 && error < 1.0E-10) {     
      if (myrank == 0) {
        cout << "Tolerance reach! Error: " << error << endl;
      }
      break;
    }

  }// for ( iter ; iter < itermax //

  // Close output stream
  fclose(otp);

  double H_over_V = real(H) / V;
  if (first_sim || H_over_V < min_H_over_V) {
    min_H_over_V = H_over_V;
    min_L = L[0];
  }
  if (myrank == 0) {
    // Output length, H, and H/V data to brent.dat
    for (int i=0; i<Dim; i++) fprintf(brent_otp, "%5.6lf ", L[i]);
    fprintf(brent_otp, "%5.6lf %5.6lf\n", real(H), H_over_V);
    fflush(brent_otp);
    // Output results to standard output
    printf("For L[0]=%lf, H=%lf and H/V=%lf\n", L[0], real(H), H_over_V);
    printf("Global minimum: L=%lf, H/V=%lf\n", L[0], min_H_over_V);
    printf("Completed L[0]=%lf simulation\n", L[0]);
  }

  if (first_sim) first_sim = 0;

  return H_over_V;
}
