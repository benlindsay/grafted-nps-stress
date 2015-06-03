#include "globals.h"

complex<double> wpl_part() ;
complex<double> wab_part() ;


complex<double> calc_H() {

  int i;

  Hcur = wpl_part() + wab_part()  ;
  if ( nD > 0.0 )
    Hcur += - nD*log(Qd) ;
  if ( nAH > 0.0 )
    Hcur += - nAH * log(Qha);
  if ( Hcur != Hcur )
    cout << Hcur << " " << wpl_part() << " " << wab_part() << endl;

  return Hcur ;
}


complex<double> wpl_part() {

  int i ;

  for ( i=0 ; i<ML ; i++ ) 
    tmp2[i] = ( kappaN <= 0.0 ? 0.0 : wpl[i] * wpl[i] * C / kappaN / 2.0 )
      - I * C * ( 1.0 - rho_surf[i] - rho_exp_nr[i] ) * wpl[i] ;

  return integ_trapPBC( tmp2 ) ;

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


