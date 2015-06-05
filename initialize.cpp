#include "globals.h"
void calc_P_constants(int);
void read_resume_files(void);
void init_particles(void);
void calc_gaa(complex<double>* , double);
void calc_gbb(complex<double>* , double);
void calc_gab(complex<double>* , double);
double get_r(int , double*) ;
double dot_prod(double[Dim], double[Dim]);
double cross_prod(double[Dim], double[Dim]);
void write_data(char*, complex<double>*); // Added for debugging
void init_negative_k(void) ;
void init_fields(void);
void explicit_nanorod(double, double, double, double[Dim], double[Dim]);
void explicit_nanosphere(double, double, double[Dim]);

void initialize_1() {

  I = complex<double>(0.0,1.0);
  M = 1; V = 1.0;
  for (int i=0; i<Dim; i++) {
    M *= Nx[i];
    V *= L[i];
    dx[i] = L[i] / double(Nx[i]);
  }
  
  N = Nda + Ndb ;
  fD = double(Nda) / double(N) ;

  if (a == -1.0) a_squared = 1.0 / double(N-1);
  else a_squared = a * a;

  if (myrank == 0) {
    cout << "V: " << V << endl;
    cout << "fD: " << fD << endl;
  }

  n_samples = 0.0;

  if (myrank == 0) {
    printf("First initialization complete\n");
    fflush(stdout);
  }
}

// Initializes fields and other variables
void initialize_2() {
  int i;
  double k2 , kv[Dim];
  double mdr2, mdr , r1[Dim], r2[Dim], dr[Dim] , x[Dim];

  // Initialize hhat
  for (i=0; i<ML; i++) {
    k2 = get_k(i , kv);
    hhat[i] = exp(-k2 * a_squared / 2.0) ;
  }

  // Initialize confining wall //
  if (do_film) {
    for (i=0; i<ML; i++) {

      // Distance from current position to z=0 plane
      r1[Dim-1] = 0.0;
      get_r(i, x);
      for (int j=0; j<Dim-1; j++)
        r1[j] = x[j];
      mdr2 = pbc_mdr2(x, r1, dr);
      mdr = sqrt(mdr2);

      rho_surf[i] = 0.5 * (1.0 - erf((fabs(mdr) - wallT) / wallXi));

    }
    write_data_bin( "rho_surf" ,  rho_surf ) ;
    fft_fwd_wrapper( rho_surf , surfH ) ;
  }
  else 
    for ( i=0 ; i<ML ; i++ )
      rho_surf[i] = surfH[i] = 0.0;

  // Initialize explicit nanorod density
  for ( i=0; i<ML ; i++ )
    rho_exp_nr[i] = exp_nrH[i] = 0.0;
  // Loop over nanorods and add nanorod density to rho_exp_nr
  for ( i=0; i<n_exp_nr; i++ ) {
    if (L_nr == R_nr)
      explicit_nanosphere(R_nr, xi_nr, exp_nr_c[i]);
    else
      explicit_nanorod(L_nr, R_nr, xi_nr, exp_nr_c[i], exp_nr_u[i]);
  }
  // If there are nanorods (and/or nanospheres) write data and take FT
  if (n_exp_nr>0) {
    write_data_bin( "rho_exp_nr" ,  rho_exp_nr ) ;
    fft_fwd_wrapper( rho_exp_nr , exp_nrH ) ;
  }
    
  // Initialize fields
  if (keep_fields && !first_sim) {
    // Do nothing if we're keeping the fields and it's not the first
    // simulation. The final fields from the previous simulation carry on to
    // the next one. (first_sim gets marked false at end of simulate() method)
    if (myrank == 0) printf("Keeping fields, not reinitializing\n");
  }
  else {
    init_fields();
    if (myrank == 0) printf("Reinitialized fields\n");
  }

  // Define the "free" volume
  Vf = V - real(integ_trapPBC( rho_surf ) + integ_trapPBC( rho_exp_nr ) );

  // Number of molecules of each component
  nD = C * Vf * ( 1.0 - phiH ) ;
  nAH = C * Vf * phiH * double( N ) / double ( Nah ) ;

  rho0 = ( nD * N + nAH * Nah ) / Vf ;

  if ( myrank == 0 ) {
    cout << "Total V_segment actual: " << nD * N + nAH * Nah << endl;
    cout << "Total V_segment theoretical: " << C * Vf * double(N)  << endl;
    cout << "V - Vf: " << V - Vf << endl;
    cout << "rho0 = " << rho0 << endl;
  }

  // Initialize Debye functions
  calc_gaa( gaa , fD );
  calc_gbb( gbb , fD );
  calc_gab( gab , fD );
  
  calc_gd( gd, double( Nah-1 ) / double( N-1 ) ) ;

  ///////////////////////////////
  // Set up bonding potentials //
  ///////////////////////////////
  if ( myrank == 0 ) 
    cout << "Setting up Gaussian bonds!" << endl; 
  for ( i=0 ; i<ML ; i++ ) {
    k2 = get_k( i , kv ) ;
    poly_bond_fft[i] = exp( -k2 / double(N-1) ) ;
  }

  iter = 0;

  // Zero all densities //
  for ( i=0 ; i<ML ; i++ )
    rhoda[i] = rhodb[i] = rhoha[i] = 0.0 ;

  if ( myrank == 0 ) {
    printf("Second initialization complete\n") ;
    fflush(stdout) ;
  }
}

// Initializes fields based on pattern specified in input file, then
// overwrites with restart file values if present.
void init_fields() {
  int sincos_dir;
  double x[Dim];
  // Initialize fields based on custom pattern specified in input file
  for (int i=0; i<ML; i++) {
    wpl[i] = wabp[i] = 0.0;
    get_r(i, x);

    // Sine and cosine pattern
    if (ic_flag[0] < -1) {
      if (real(rho_surf[i] + rho_exp_nr[i]) > 0.01)
        // -1 if there's some wall or nanorod present
        wabm[i] = -1;
      else {
        wabm[i] = 1.0;
        for (int j=0; j<2; j++) { // initial condition flag sets
          sincos_dir = ic_dir[j];
          if (ic_flag[j] == -2) // cosine flag
            wabm[i] *= ic_pre[j]*tanh( cos(2.0 * PI * ic_period[j]
                  * x[sincos_dir] / L[sincos_dir])/0.2 );
          else if (ic_flag[j] == -3) // sin flag
            wabm[i] *= ic_pre[j]*tanh( sin(2.0 * PI * ic_period[j]
                  * x[sincos_dir] / L[sincos_dir])/0.2 );
          else if (ic_flag[j] < -4) { // -4 flag is just a factor of 1
            cout << "Initial condition flags smaller than -4 not supported"
              << endl;
            exit(1);
          }
        } // j (initial condition flag sets)
      }
    }
    // Random fields
    else if (ic_flag[0] == -1) 
      wabm[i] = ic_pre[0] * ran2() ;
    // Zero fields
    else 
      wabm[i] = 0.0 ;
  } 

  // Overwrite fields with restart files if available
  read_resume_files();
}

// Initialize 1 explicit nanorod
void explicit_nanorod(double L, double R, double xi,
                      double center[Dim], double u[Dim]) {
  int i_global;
  int nn[Dim];
  double u_dot_r, u_cross_r;
  double dr[Dim], x[Dim];

  for (int i=0; i<ML; i++) {
    // Get distance from nanorod center to current position (dr)
    i_global = unstack_stack( i );
    unstack( i_global, nn);
    for (int j=0; j<Dim; j++)
      x[j] = double(nn[j])*dx[j];
    pbc_mdr2( x, center, dr );

    // Compute explicit nanorod density
    u_dot_r = dot_prod(u, dr);
    u_cross_r = cross_prod(u, dr);
    rho_exp_nr[i] += 0.25 * erfc((fabs(u_dot_r)-0.5*L)/xi)
      * erfc((fabs(u_cross_r)-R)/xi);
  }
}

// Initialize 1 explicit nanosphere
void explicit_nanosphere(double R, double xi, double center[Dim]) {
  int i_global;
  int nn[Dim];
  double dr2, dr_abs, dr[Dim], x[Dim];

  for (int i=0; i<ML; i++) {
    // Get distance from nanorod center to current position (dr)
    i_global = unstack_stack( i );
    unstack( i_global, nn);
    for (int j=0; j<Dim; j++)
      x[j] = double(nn[j])*dx[j];
    dr2 = pbc_mdr2( x, center, dr );
    dr_abs = sqrt(dr2);

    // Compute explicit nanorod density
    rho_exp_nr[i] += 0.5 * erfc((dr_abs-R)/xi);
  }
}
