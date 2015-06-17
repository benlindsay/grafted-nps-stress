#include "globals.h"

complex<double> solvent_density( complex<double>* , complex<double>*, double ) ;
complex<double> diblock_discrete( complex<double>* , complex<double>*,
    complex<double>** , complex<double>** , int , int );
void integrate_diblock_discrete( complex<double>** , complex<double>** ,
    complex<double>, complex<double> , complex<double>* , complex<double>* ,
    complex<double>* , complex<double>*, int , int );
complex<double> homopolymer_discrete( complex<double>*, complex<double>**,
    int );
void integrate_homopoly_discrete(  complex<double>** , complex<double>,
    complex<double> , complex<double>* , complex<double>*,  int );
complex<double> brush_homopoly_discrete(complex<double>*, complex<double>*, 
    complex<double>**, complex<double>**, complex<double>*, int);
void integrate_brush_homopoly_discrete(complex<double>**, complex<double>**, 
    complex<double>, complex<double>, complex<double>*, complex<double>*, int );
complex<double> grafted_nanoparticles( complex<double>*, complex<double>*,
    complex<double>*, complex<double>**, complex<double>**, complex<double>*,
    complex<double>*, complex<double>*, double , double , int ) ;

void calc_poly_density() {

  int i;

  // Iterate over all points
  for (i=0; i<ML; i++) {

    // Add all fields that make up wa and wb
    wa[i] = (I * (wpl[i] + wabp[i]) - wabm[i]) / double(N);
    wb[i] = (I * (wpl[i] + wabp[i]) + wabm[i]) / double(N);

    // If doing a film calculation, add surface density times appropriate
    // lambda to wa and wb
    if (do_film) {
      // Get global index from index on current processor
      int i_global = unstack_stack(i);
      int nn[Dim];
      // Fill nn with global x, y [, and z] indices
      unstack(i_global, nn);
      // Get x[Dim-1] (z in 3D) position
      double z = dx[Dim-1] * double(nn[Dim-1]);
      double lamA, lamB;
      // Set lamA and lamB to the values associated with whichever surface
      // is closer
      if (z > L[Dim-1]/2.0) {
        lamA = top_wall_lamA;
        lamB = top_wall_lamB;
      }
      else {
        lamA = bot_wall_lamA;
        lamB = bot_wall_lamB;
      }
      // Add surface density times appropriate lambda to each field
      wa[i] += lamA * rho_surf[i];
      wb[i] += lamB * rho_surf[i];
    }

    // If doing explicit nanorods, add nanorod interaction terms to fields
    if (n_exp_nr > 0) {
      wa[i] += exp_nr_chiAPN / double(N) * rho_exp_nr[i];
      wb[i] += exp_nr_chiBPN / double(N) * rho_exp_nr[i];
    }
  }

  // Convolve wa and wb with h to get smeared fields smwa and smwb
  fft_fwd_wrapper(wa, smwa);
  fft_fwd_wrapper(wb, smwb);
  for (i=0; i<ML; i++) {
    smwa[i] *= hhat[i];
    smwb[i] *= hhat[i];
  }
  fft_bck_wrapper(smwa, smwa);
  fft_bck_wrapper(smwb, smwb);

  ///////////////////////////////
  // Polymer matrix propagator //
  ///////////////////////////////
  if (nD > 0.0) {
    Qd = diblock_discrete(smwa, smwb, qd, qddag, N, Nda);
    integrate_diblock_discrete(qd, qddag, Qd, nD, rhoda, rhodb,
                               smwa, smwb, N, Nda);
  }
  else 
    Qd = 1.0;
  
  if (nAH > 0.0) {
    Qha = homopolymer_discrete(smwa, qha, Nah);
    integrate_homopoly_discrete(qha, Qha, nAH, rhoha, smwa, Nah);
  }
  else 
    Qha = 1.0;
} // End calc_poly_density

///////////////////////////////////////////////////////
// Calculate density of a discrete diblock copolymer //
///////////////////////////////////////////////////////
void integrate_diblock_discrete(complex<double>** q, complex<double> **qdag, 
    complex<double> Q, complex<double> nK, complex<double>* rda, 
    complex<double>* rdb, complex<double>* wA,
    complex<double>* wB, int Ns, int Na) {

  int i, j;

  complex<double> factor;

  factor = nK / Q / V;

  for (i=0; i<ML; i++) {
    rda[i] = rdb[i] = 0.0;
    for (j=0; j<Ns; j++) {
      if (Na > 0 && j < Na)
        rda[i] += q[j][i] * qdag[Ns-j-1][i];
      else
        rdb[i] += q[j][i] * qdag[Ns-j-1][i];
    }

    rda[i] *= exp(wA[i]) * factor;
    rdb[i] *= exp(wB[i]) * factor;

  }//for ( i=0 ; i<ML
} // End integrate_diblock_discrete

void integrate_homopoly_discrete(complex<double>** q,
    complex<double> Q, complex<double> nK, complex<double>* rh,
    complex<double>* W, int N) {

  int i, j;

  complex<double> factor;

  factor = nK / Q / V;

  for (i=0; i<ML; i++) {
    rh[i] = 0.0;

    for (j=0; j<N; j++) {
        rh[i] += q[j][i] * q[N-j-1][i];
    }

    rh[i] *= exp(W[i]) * factor;
  } // for ( i=0 ; i<ML
} // End integrate_homopoly_discrete
