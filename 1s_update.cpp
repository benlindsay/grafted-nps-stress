#include "globals.h"
void generate_1s_noise( complex<double>* , double ) ;


void update_1s( ) {

  int i ;
  complex<double> evar , F, A , numer, denom ;

  if (nD > 0.0) {
    fft_fwd_wrapper(rhoda, rhoda);
    fft_fwd_wrapper(rhodb, rhodb);
  }
   
  if (nAH > 0.0) 
    fft_fwd_wrapper(rhoha, rhoha); 

  if (sigma > 0.0 && nFP > 0.0)
    fft_fwd_wrapper(rhoga, rhoga);

  if (iter==0)
    write_data_bin("pre_pre_rhoga_exp", rhoga_exp);

  if (sigma > 0.0 && n_exp_nr > 0.0)
    fft_fwd_wrapper(rhoga_exp, rhoga_exp);

  if (do_fld_np)
    fft_fwd_wrapper(rho_fld_np, rho_fld_np);
 
  fft_fwd_wrapper(wpl, wpl);

  if (chiN > 0.0) {
    fft_fwd_wrapper(wabm, wabm);
    fft_fwd_wrapper(wabp, wabp);
  }

  if (do_CL) 
    generate_1s_noise(etap, lam_pl);

  // Update w+ field //
  for (i=0; i<ML; i++) {
    // The 1 at the k=0 mode takes care of the delta function arising from
    // taking a fourier transform of a constant
    if (i == 0 && myrank == 0)
      evar = 1.0 ;
    else
      evar = 0.0 ;
    
    F = (kappaN <= 0.0 ? 0.0 : C*wpl[i]/kappaN) 
      + I * C * (surfH[i] + exp_nrH[i] - evar)
      + I * rho_fld_np[i] / double(N)
      + I * hhat[i] / double(N) * (rhoha[i] + rhoda[i] +
                                   rhoga[i] + rhoga_exp[i] + rhodb[i]);
    A = (kappaN <= 0.0 ? 0.0 : C/kappaN)
      + nD * double(N) * hhat[i] * hhat[i] * (gaa[i] + 2.0 * gab[i] + gbb[i]) / V
      + ( nAH * double(Nah * Nah) / double(N)
          + ng_per_np * nP * double(Ng * Ng) / double(N) )
        * hhat[i] * hhat[i] * gd[i] / V;
    numer = wpl[i] - lam_pl * (F - A * wpl[i]);
    if (do_CL) 
      numer += etap[i];
    denom = 1.0 + lam_pl * A;
    wpl[i] = numer / denom; 
  }

  if (iter==0) {
    write_data_bin("pre_rhodb", rhodb);
    write_data_bin("pre_rhoda", rhoda);
    write_data_bin("pre_rhoha", rhoha);
    write_data_bin("pre_rhoga", rhoga);
    write_data_bin("pre_rhoga_exp", rhoga_exp);
    write_data_bin("pre_rho_fld_np", rho_fld_np);
    write_data_bin("pre_exp_nrH", exp_nrH);
  }
  
  // Update AB fields //
  if (chiN > 0.0) {
    if (do_CL) {
      generate_1s_noise(etap, lam_pl);
      generate_1s_noise(etam, lam_mi);
    }

    for ( i=0 ; i<ML ; i++ ) {
      // AB+ //
      F = 2.0 * C * wabp[i] / chiN
        + I * hhat[i] / double(N) * (rhoda[i] + rhodb[i] +
                                     rhoha[i] + rhoga[i] + rhoga_exp[i] );
      // If NP is A chemistry instead of neutral, add NP contributions
      if (np_chem == 1)
        F += I * rho_fld_np[i] / double(N) + I * C * exp_nrH[i];

      A = 2.0 * C / chiN 
        + nD*double(N)*hhat[i]*hhat[i] * (gaa[i] + 2.0*gab[i] + gbb[i]) / V 
        + ( nAH * double(Nah * Nah) / double(N)
            + ng_per_np * nP * double(Ng * Ng) / double(N) )
          * hhat[i] * hhat[i] * gd[i] / V;
      numer = wabp[i] - lam_pl * ( F - A * wabp[i] ) ;
      if ( do_CL )
        numer += etap[i] ;
      denom = 1.0 + lam_pl * A ;
      wabp[i] = numer / denom ;

      // AB- //
      F = 2.0 * C * wabm[i] / chiN
          + hhat[i] / double(N) * ( rhodb[i] - rhoda[i] - rhoha[i] -
                                    rhoga[i] - rhoga_exp[i]);
      // If NP is A chemistry instead of neutral, add NP contributions
      if (np_chem == 1)
        F += - I * rho_fld_np[i] / double(N) - I * C * exp_nrH[i];
      A = 2.0 * C / chiN ;
      numer = wabm[i] - lam_mi * ( F - A * wabm[i] ) ;
      if ( do_CL ) 
        numer += etam[i] ;
      denom = 1.0 + lam_mi * A ;
      wabm[i] = numer / denom ;
    }
  }
    
  else {
    for ( i=0 ; i<ML ; i++ )
      wabp[i] = wabm[i] = 0.0 ;
  }

  if (iter==0) {
    write_data_bin("post_rhodb", rhodb);
    write_data_bin("post_rhoda", rhoda);
    write_data_bin("post_rhoha", rhoha);
    write_data_bin("post_rhoga", rhoga);
    write_data_bin("post_rhoga_exp", rhoga_exp);
    write_data_bin("post_rho_fld_np", rho_fld_np);
    write_data_bin("post_exp_nrH", exp_nrH);
  }

  fft_bck_wrapper( wpl , wpl ) ;

  if ( chiN > 0.0 ) {
    fft_bck_wrapper( wabm , wabm ) ;
    fft_bck_wrapper( wabp , wabp ) ;
  }

  calc_poly_density() ;
}


// Generates Gaussian noise in k-space with appropriate statistics //
// for the 1s updating scheme. 
void generate_1s_noise( complex<double> *et, double lambda ) {

  int i ;
  double scale = sqrt( 2.0 * lambda ) ;
  for ( i=0 ; i<Dim ; i++ )
    scale /= sqrt( dx[i] ) ;

  for ( i=0 ; i<ML ; i++ ) 
    et[i] = scale * gasdev2() ;

  fft_fwd_wrapper( et , et ) ;

}
