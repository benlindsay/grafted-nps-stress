#include "globals.h"

complex<double> wpl_part() ;
complex<double> wab_part() ;


complex<double> calc_H() {

  int i;

  Hcur = wpl_part() + wab_part()  ;

  // Add diblock part if applicable
  if (nD > 0.0)
    Hcur += - nD*log(Qd) ;

  // Add homopolymer part if applicable
  if (nAH > 0.0)
    Hcur += - nAH * log(Qha);

  // Add field-based nanoparticle part if applicable. Note that smwp_min,
  // which was originally the minimum value of the smeared wp field (smwp),
  // is subtracted off here because when we subtracted smwp_min from smwp,
  // that made the value stored in Qp actually Qp * exp(smwp_min).
  if (do_fld_np)
    Hcur += - nP * (log(Qp) - smwp_min);

  // Add contribution of chains grafted to explicit nanoparticles
  if (n_exp_nr > 0 && sigma > 0.0)
    Hcur += - double(n_exp_nr) * ng_per_np * real(Qga_exp);

  // Exit if H is NaN
  if ( Hcur != Hcur ) {
    if (myrank == 0) {
      printf("Hcur is NaN!\n");
      printf("real(Hcur) = %lf\n", real(Hcur));
      printf("real(wpl_part) = %lf\n", real(wpl_part()) );
      printf("real(wab_part)=%lf\n", real(wab_part()) );
      printf("real(-nD*log(Qd))=%lf\n", real(-nD*log(Qd)) );
      printf("real(-nAH*log(Qha))=%lf\n", real(-nAH*log(Qha)) );
      printf("real(-nP*[log(Qp)-smwp_min])=%lf\n",
              real(-nP*(log(Qp)-smwp_min)) );
      printf("real(-n_exp_grafts*Qga_exp) = %lf\n",
              - double(n_exp_nr) * ng_per_np * real(Qga_exp) );
      printf("real(Qga_exp) = %lf\n", real(Qga_exp));
    }
    exit(1);
  }

  return Hcur ;
}

complex<double> wpl_part() {

  int i;

  for (i=0; i<ML; i++) {
    tmp2[i] = ( kappaN<=0.0 ? 0.0 : wpl[i]*wpl[i]*C/kappaN/2.0 )
      - I * C * ( 1.0 - rho_surf[i] - rho_exp_nr[i] ) * wpl[i];
  }

  return integ_trapPBC(tmp2);

}

complex<double> wab_part() {

  if ( chiN == 0.0 ) {
     return 0.0 ;
  }
  else {
    int i ;
 
    for ( i=0 ; i<ML ; i++ ) 
      tmp[i] = ( wabp[i] * wabp[i] + wabm[i] * wabm[i] ) * C / chiN ;
 
    return integ_trapPBC( tmp ) ;
  }
}


