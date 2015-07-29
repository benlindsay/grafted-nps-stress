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
void generate_smwp_iso(complex<double>*, complex<double>*,
    complex<double>*);
void generate_smwp_aniso(complex<double>*, complex<double>***,
    complex<double>***);
void integ_sphere_posits(complex<double>***, complex<double>*);
complex<double> np_density_sphere(void);
complex<double> np_density_rod(void);

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

  //////////////////////////////////////////////////////////////
  // Nanoparticle partition function and density calculations //
  //////////////////////////////////////////////////////////////

  // Generate smwp, which is the wa field convolved with Gamma (assuming the
  // nanoparticles are chemically identical to A).
  if (do_fld_np) {
    if (np_type == 1) {
      // For spherical particles, no orientation dependence
      generate_smwp_iso(wa, Gamma_iso, smwp_iso);
      Qp = np_density_sphere();
    }
    else if (np_type == 2) {
      // Orientation dependence for rods or other anisotropic particles
      generate_smwp_aniso(wa, Gamma_aniso, smwp_aniso);
      Qp = np_density_rod();
    }
  }
  else {
    Qp = 1.0;
    smwp_min = 0.0;
  }
  
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
    } // j

    rda[i] *= exp(wA[i]) * factor;
    rdb[i] *= exp(wB[i]) * factor;

  } // i
} // integrate_diblock_discrete

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

// Generate particle field convolved with Gamma for isotropic particles
//
// Given a w field (i.e. wa) in real space and smGamma, the Gamma function
// in k-space, this generates smwp, the particle field convolved with Gamma,
// as a function of position
void generate_smwp_iso(complex<double>* w, complex<double>* smGamma,
                         complex<double>* smwp) {

  int i, j, k;

  // Fourier transform w to prepare for convolution
  fft_fwd_wrapper(w, w);
  for(i=0; i<ML; i++){
    // Multiply w and Gamma in k-space. Extra V is necessary
    smwp[i] = w[k] * smGamma[i] * V;  
    // Inverse fourier transform smwp to conclude convolution
    fft_bck_wrapper(smwp, smwp);
  } // i

  // Shift smwp so the minimum is 0. This avoids numerical difficulties.
  // This is corrected for when Qp is calculated. First, find current
  // processor's minimum:
  double smwp_min_local = 1000.0;
  for (i=0; i<ML; i++) {
    if (real(smwp[i]) < smwp_min_local)
      smwp_min_local = real(smwp[i]) ;
  }

#ifdef PAR
  // Next, find the global minimum. Replace smwp_min in this processor with
  // that global minimum across all processors.
  MPI_Allreduce(&smwp_min_local, &smwp_min, 1,
                MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  // Finally, subtract off the global minimum from all values in smwp
  for (i=0; i<ML; i++) {
    smwp[i] -= smwp_min;
  }

} // generate_smwp_iso

// Generate particle field convolved with Gamma for anisotropic particles
//
// Given a w field (i.e. wa) in real space and smGamma, the Gamma function
// in k-space, this generates smwp, the particle field convolved with Gamma,
// as a function of orientation and position
void generate_smwp_aniso(complex<double>* w, complex<double>*** smGamma,
                         complex<double>*** smwp) {

  int i, j, k;

  // Fourier transform w to prepare for convolution
  fft_fwd_wrapper(w, w);
  for(i=0; i<Nu; i++){
    for(j=0; j<2*Nu; j++){
      for(k=0; k<ML; k++){
        // Multiply w and Gamma in k-space. Extra V is not necessary because
        // it was already multiplied into Gamma during initialization. It's
        // cheaper that way.
        smwp[i][j][k] = w[k] * smGamma[i][j][k];
      } // k
      // Inverse fourier transform smwp to conclude convolution
      fft_bck_wrapper(smwp[i][j], smwp[i][j]);
    } // j  
  } // i

  // Shift smwp so the minimum is 0. This avoids numerical difficulties.
  // This is corrected for when Qp is calculated. First, find current
  // processor's minimum:
  double smwp_min = 1000.0;
  for (i=0; i<Nu; i++) {
    for (j=0; j<2*Nu; j++) {
      for (k=0; k<ML; k++) {
        if (real(smwp[i][j][k]) < smwp_min)
          smwp_min = real(smwp[i][j][k]) ;
      }
    }
  }

#ifdef PAR
  // Next, find the global minimum. Replace smwp_min in this processor with
  // that global minimum across all processors.
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Allreduce(&smwp_min, &smwp_min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
#endif

  // Finally, subtract off the global minimum from all values in smwp
  // Store the exponential of the negative smeared field for use in calculating
  // the partition function and such
  for (i=0; i<Nu; i++) {
    for (j=0; j<2*Nu; j++) {
      for (k=0; k<ML; k++) {
        smwp[i][j][k] -= smwp_min;
        exp_neg_smwp[i][j][k] = exp(-smwp[i][j][k]);
      }
    }
  }

} // generate_smwp_aniso

// Calculate density of field-based nanospheres
complex<double> np_density_sphere() {
  return 0.0;
}

// Calculate density of field-based nanorods
complex<double> np_density_rod() {
  int i, j, k;
  complex<double> Qrod;

  // Store integral over u of exp(-smwp) in rho_fld_np_c. rho_fld_np_c will
  // still be missing a factor of 1 / (4*PI*V*Qp)
  integ_sphere_posits(exp_neg_smwp, rho_fld_np_c);

  // Store integral over u and r of exp(-smwp) in the variable Qrod.
  // Qrod will still need a factor of 1 / (4*PI*V) before it represents Qp,
  // so currently Qrod = 4*PI*V*Qp (see Koski JCP Paper Eq. 34)
  Qrod = integ_trapPBC(rho_fld_np_c);

  for (i=0; i<ML; i++) {
    // Equation: rho_fld_np_c = nP / (4 * PI * V * Qp) * exp(-smwp)
    // (see Eq 35 in Koski JCP paper)
    rho_fld_np_c[i] *= nFP / Qrod;

    // Set rho_fld_np to 0 to prepare for the += stuff coming up
    rho_fld_np[i] = 0.0;
  }

  for (i=0; i<Nu; i++) {
    for (j=0; j<2*Nu; j++) {
      fft_fwd_wrapper(exp_neg_smwp[i][j], tmp);
      for (k=0; k<ML; k++) {
        rho_fld_np[k] += tmp[k] * Gamma_aniso[i][j][k] * theta_weights[i]
                         * phi_weights[j];
      }
    }
  }
  
  for (i=0; i<ML; i++) {
    // Qrod still represents 4*PI*V*Qp here
    rho_fld_np[i] *= nFP / Qrod;
  }

  // Convert rho_fld_np back to real space
  fft_bck_wrapper(rho_fld_np, rho_fld_np);

  // Now multiply Qrod by 1/(4*PI*V) so it actually represents Qp
  Qrod = Qrod / (4.0 * PI * V);

  // Sanity check
  complex<double> np_check = integ_trapPBC(rho_fld_np_c);
  complex<double> npVp_check = integ_trapPBC(rho_fld_np);
  complex<double> C_check_easy = (nD + nP*V_1_fld_np/double(N)) / V;
  complex<double> C_check_hard = (nD + npVp_check / double(N)) / V;
  if (myrank==0) {
    printf("nP=%lf, np_check=%lf\n", nP, real(np_check));
    printf("nP*Vp = %lf, npVp_check=%lf\n",
        nP*V_1_fld_np, real(npVp_check) );
    printf("C = %lf, C_check_easy=%lf, C_check_hard=%lf\n",
                 C,  real(C_check_easy), real(C_check_hard) );
  }
  
  return Qrod;
}
