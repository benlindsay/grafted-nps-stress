#include "globals.h"
/*
 * Note: "grafts" is expected to be the *single particle*, *normalized*
 * grafted density *in k space*. Its real-space integral should be unity. 
 * The grafted particle field will be constructed as:
 * 
 * WP_tot(r) = WP(r) - n_g / n_P * ( grafts * log( q[Ng-1] )(r)
 * 
 * where grafts is convolved with log(q).
 *
 */
complex<double> shift_field( complex<double>* ) ;

complex<double> grafted_nanoparticles( 
    complex<double> *WP,  // Potential felt by particle cores
    complex<double> *WG,  // Potential felt by grafts
    complex<double> *graft_pts,  // Graft distribution for one particle
    complex<double> *gammaNP,  // nanoparticle shape function
    complex<double> **q, complex<double> **qdag, // Propagators for grafts
    complex<double> *rp , complex<double> *smrp, // Nanoparticle densities
    complex<double> shift_WP, // Shift of chx pot field of the NPs
    complex<double> *rhg, // Grafted chain densities
    complex<double> npar,    // number of particles
    double ngrafts_per_np, // Number of grafted chains 
    int Ng // Length of grafts
    )

{

  int i,j ;

  complex<double> Qgp ;

  // Distribution of bare particles //
  if ( sigma == 0.0 ) {

    shift_WP = 0.0 ; //shift_field( WP  ) ;

    for ( i=0 ; i<ML ; i++ )
      tmp[i] = exp( -WP[i] ) ;

    Qgp = integ_trapPBC( tmp ) / V ;
 
    for ( i=0 ; i<ML ; i++ ) 
      rp[i] = tmp[i] * npar / ( Qgp * V ) ;

    fft_fwd_wrapper( rp , smrp ) ;
    for ( i=0 ; i<ML ; i++ ) 
      smrp[i] = gammaNP[i] * smrp[i] ;

    fft_bck_wrapper( smrp , smrp  ) ;

    return Qgp ;

  }

  // Grafted nanoparticles //
  else {

    // First, get the propagator from the free end
    for ( i=0 ; i<ML ; i++ ) 
      q[0][i] = exp( -WG[i] ) ;

    for ( j=1 ; j<Ng ; j++ ) {

      fft_fwd_wrapper( q[j-1] , q[j]  ) ;

      for ( i=0 ; i<ML ; i++ ) 
        q[j][i] *= poly_bond_fft[i] ;

      fft_bck_wrapper( q[j] , q[j] ) ;

      for ( i=0 ; i<ML ; i++ )
        q[j][i] *= exp( -WG[i] ) ;

    }

    // Construct the full particle field 
    for ( i=0 ; i<ML ; i++ ) 
      tmp[i] = ngrafts_per_np * log( q[Ng-1][i] ) ;

    fft_fwd_wrapper( tmp, tmp ) ;
    for ( i=0 ; i<ML ; i++ ) 
      tmp[i] *= graft_pts[i] ;
    fft_bck_wrapper( tmp, tmp ) ;

    for ( i=0 ; i<ML ; i++ ) 
      WP[i] = WP[i] - tmp[i] ;

    shift_WP = 0.0 ; //shift_field( WP ) ;

    // Calculate the particle partition function
    for ( i=0 ; i<ML ; i++ )
      tmp[i] = exp( - WP[i] ) ;

    Qgp = integ_trapPBC( tmp ) / V ;

    // Calculate the particle core densities //
    for ( i=0 ; i<ML ; i++ ) 
      rp[i] = tmp[i] * npar / ( Qgp * V ) ;

    fft_fwd_wrapper( rp , tmp2 ) ;
    for ( i=0 ; i<ML ; i++ ) 
      smrp[i] = gammaNP[i] * tmp2[i] ;

    fft_bck_wrapper( smrp , smrp ) ;

    // Calculate the total graft density
    for ( i=0 ; i<ML ; i++ ) 
      tmp2[i] = graft_pts[i] * tmp2[i] / npar ;

    fft_bck_wrapper( tmp2 , tmp2 ) ;

    // Include contribution from explicit particle if needed //
    if ( explicit_np ) {
      for ( i=0 ; i<ML ; i++ ) 
        tmp2[i] += expl_grafts[i] ;
    }

    // Initial condition for the complimentary graft propagator
    for ( i=0 ; i<ML ; i++ ) 
      qdag[0][i] = tmp2[i] * exp( -WG[i] ) / q[Ng-1][i] ;

    // Rest of the graft propagator
    for ( j=1 ; j<Ng ; j++ ) {

      fft_fwd_wrapper( qdag[j-1] , qdag[j] ) ;

      for ( i=0 ; i<ML ; i++ ) 
        qdag[j][i] *= poly_bond_fft[i] ;

      fft_bck_wrapper( qdag[j] , qdag[j] ) ;

      for ( i=0 ; i<ML ; i++ )
        qdag[j][i] *= exp( -WG[i] ) ;

    }

    // Finally, calculate the grafted chain density //
    for ( i=0 ; i<ML ; i++ ) {
      rhg[i] = 0.0 ;
      for ( j=0 ; j<Ng ; j++ ) 
        rhg[i] += q[j][i] * qdag[Ng-j-1][i] ;

      rhg[i] *= ngrafts_per_np * npar * exp( +WG[i] ) ;
    }

    return Qgp ;

  }// if (sigma > 0)

}


complex<double> shift_field( complex<double> *w  ) {

  double wmin, wtmp ;
  int i;

  wtmp = 12342.0 ;

  for ( i=0 ; i<ML ; i++ )
    if ( real( w[i] ) < wtmp )
      wtmp = real( w[i] ) ;

  MPI_Allreduce( &wtmp, &wmin, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD ) ;

  for ( i=0 ; i<ML ; i++ )
    w[i] -= wmin ;

  return complex<double>(wmin, 0.0) ;

}


void initialize_gNP() {

  int i, j;
  double k2, kv[Dim], r1[Dim], r2[Dim], dr[Dim], mdr2, mdr ;

  for ( j=0 ; j<Dim ; j++ ) 
    r1[j] = 0.0 ;

  for ( i=0 ; i<ML; i++ ) {
    get_r( i , r2  ) ;
    mdr2 = pbc_mdr2( r1, r2, dr ) ;

    mdr = sqrt( mdr2 ) ;
    
    double arg = (mdr - Rp) / Xi ;
    gammap[i] = double(N) * C * 0.5 * ( 1.0 - erf( arg ) ) ;
    
    gamma_vir[i] = -gammap[i] / V + exp( -arg*arg ) * mdr * C * double(N) /
                                    (sqrt(PI) * double(Dim) * Xi * V) ;

    double exp_arg = ( mdr - Rp ) / Xi ;
    //double exp_arg = ( mdr - Rp - Xi ) / Xi ;
    grafts[i] = exp( -exp_arg * exp_arg ) ;

  }

  if ( explicit_np ) {
    for ( j=0 ; j<Dim ; j++ )
      r1[j] = expl_x[j] ;
    for ( i=0 ; i<ML ; i++ ) {
      get_r( i , r2 ) ;
      mdr2 = pbc_mdr2( r1, r2, dr ) ;

      mdr = sqrt( mdr2 ) ;

      double arg = (mdr - Rp ) / Xi ;

      expl_rhop[i] = double(N) * C * 0.5 * ( 1.0 - erf( arg ) ) ;

      double exp_arg = ( mdr - Rp - Xi ) / Xi ;
      expl_grafts[i] = exp( -exp_arg * exp_arg ) ;
    }
  }
  else
    for ( i=0 ; i<ML ; i++ )
      expl_rhop[i] = expl_grafts[i] = 0.0 ;

  cout << "Rp, Xi = " << Rp << " " << Xi << endl;

  Vnp = real( integ_trapPBC( gammap  ) ) ;
  
  write_data_bin( "grafts", grafts ) ;
  write_data_bin( "gamma" , gammap ) ;

  if ( explicit_np ) {
    write_data_bin( "expl_rhop" , expl_rhop ) ; 
    write_data_bin( "expl_grafts" , expl_grafts ) ;

    fft_fwd_wrapper( expl_rhop , expl_rhop ) ;
  }

  fft_fwd_wrapper( gammap , gammap  ) ; 

  fft_fwd_wrapper( gamma_vir , gamma_vir  ) ; 

  complex<double> norm ;
  if ( sigma > 0.0 ) {
    norm = integ_trapPBC( grafts  ) + integ_trapPBC( expl_grafts ) ;
    fft_fwd_wrapper( grafts, grafts ) ;
  }

  for ( i=0 ; i<ML ; i++ ) {
    gammap[i] *= V ;
    gamma_vir[i] *= V ;
    if ( sigma > 0.0 ) {
      grafts[i] *= V / norm ;
      expl_grafts[i] *= 1.0 / norm ;
    }
  }

  cout << "Nanoparticle core volume: " << Vnp << endl;

  if ( sigma > 0.0 ) {
    if ( Dim == 2 )
      ng_per_np = sigma * 2.0 * PI * Rp ;
    else if ( Dim == 3 )
      ng_per_np = sigma * 4.0 * PI * Rp * Rp ;

    Vnp += double( Ng ) * ng_per_np ;

    cout << "Volume of core + grafts: " << Vnp << endl;
    cout << "grafts per NP: " << ng_per_np << endl << endl;

  }

  else
    ng_per_np = 0.0 ;

}

complex<double> gnp_dQdV( 
    complex<double> *WP,  // Bare particle potential
    complex<double> *SMWP,  // Potential felt by particle cores
    complex<double> *WG,  // Potential felt by grafts
    complex<double> *grafts,  // Graft distribution for one particle
    complex<double> *gammaNP,  // nanoparticle shape function
    complex<double> *gam_vir,  // nanoparticle shape function
    complex<double> **q, complex<double> **qdag, // Propagators for grafts
    complex<double> *rp , complex<double> *smrp, // Nanoparticle densities
    complex<double> *rhg, // Grafted chain densities
    double ngrafts_per_np, // Number of grafted chains 
    complex<double> npar,    // number of particles
    int Ng // Length of grafts
    )
{

  if ( sigma > 0.0 )
    die("dQdV not set up to deal with grafted particles!");

  int i, j ;
  complex<double> dQdV ;
  
  fft_fwd_wrapper( WP, tmp ) ;

  for ( i=0 ; i<ML ; i++ ) {
    tmp[i] *= gam_vir[i] ;
  }
  
  fft_bck_wrapper( tmp, tmp  ) ;

  for ( i=0 ; i<ML ; i++ )
    tmp[i] *= exp( -SMWP[i] ) ;

  dQdV = integ_trapPBC( tmp ) / V ;

  return dQdV ;

}
