#include "globals.h"
#include <cmath>
using namespace std ;

void calc_eigenvals( double**, double* ) ;
double calc_trace( double** ) ;
double calc_determ( double** ) ;


void inv_mat(double a[3][3], double ainv[3][3]) {

  double den;
  int i,j;
        
  den = a[0][0]*a[1][1]*a[2][2] - a[0][0]*a[1][2]*a[2][1] -
    a[1][0]*a[0][1]*a[2][2] + a[1][0]*a[0][2]*a[2][1] +
    a[2][0]*a[0][1]*a[1][2] - a[2][0]*a[0][2]*a[1][1];

  ainv[0][0] = a[1][1]*a[2][2]-a[1][2]*a[2][1];
  ainv[0][1] = a[0][2]*a[2][1]-a[0][1]*a[2][2];
  ainv[0][2] = a[0][1]*a[1][2]-a[0][2]*a[1][1];

  ainv[1][0] = a[1][2]*a[2][0]-a[1][0]*a[2][2];
  ainv[1][1] = a[0][0]*a[2][2]-a[0][2]*a[2][0];
  ainv[1][2] = a[0][2]*a[1][0]-a[0][0]*a[1][2];

  ainv[2][0] = a[1][0]*a[2][1]-a[1][1]*a[2][0];
  ainv[2][1] = a[0][1]*a[2][0]-a[0][0]*a[2][1];
  ainv[2][2] = a[0][0]*a[1][1]-a[0][1]*a[1][0];

  for (i=0; i<3; i++)
    for (j=0; j<3; j++)
      ainv[i][j] *= (1.0 / den);
                              
                              
}


double calc_determ( double **A ) {

  double D;

  D = A[0][0] * ( A[1][1] * A[2][2] - A[1][2] * A[2][1] ) +
      A[0][1] * ( A[1][2] * A[2][0] - A[2][2] * A[0][1] ) + 
      A[0][2] * ( A[0][1] * A[2][1] - A[1][1] * A[2][0] ) ;

  return D;

}

double calc_trace( double **A ) {

  return A[0][0] + A[1][1] + A[2][2] ;

}

void calc_eigenvals( double **A , double *lam ) {

  double p, q, **B, **qI, r, phi ;
  double d1, d2, d3 ;
  int i, j;
  
  B = new double*[3] ;
  qI = new double*[3] ;
  
  for ( i=0 ; i<3 ; i++ ) {
    B[i] = new double [3] ;
    qI[i] = new double[3] ;
  }


  p = A[0][1] * A[0][1] + A[0][2]*A[0][2] + A[1][2] * A[1][2] ;

  // A already diagonal //
  if ( p == 0.0 ) {
    lam[0] = A[0][0] ;
    lam[1] = A[1][1] ;
    lam[2] = A[2][2] ;
  }


  
  else {
    q = calc_trace(A) / 3.0 ;
    
    for ( i=0 ; i<3 ; i++ ) {
      for ( j=0 ; j<3 ; j++ ) 
        qI[i][j] = 0.0 ;
      qI[i][i] = q ;
    }

    d1 = A[0][0] - q ;
    d2 = A[1][1] - q ;
    d3 = A[2][2] - q ;

    p = d1*d1 + d2*d2 + d3*d3 + 2 * p ;

    p = sqrt( p / 6.0 ) ;

    for ( i=0 ; i<3 ; i++ ) 
      for ( j=0 ; j<3 ; j++ )
        B[i][j] = ( A[i][j] - qI[i][j] ) / p ;

    r = calc_determ( B ) / 2.0 ;

    if ( r <= -1.0 )
      phi = PI / 3.0 ;
    else if ( r >= 1.0 )
      phi = 0.0 ;
    else
      phi = acos( r ) / 3.0 ;
 
    // Set up so that lam[0] >= lam[1] >= lam[2]
    lam[0] = q + 2.0 * p * cos( phi ) ;
    lam[2] = q + 2.0 * p * cos(phi + PI * 2.0 / 3.0 ) ;
    lam[1] = 3.0 * q - lam[0] - lam[2] ;

  }

  for (i = 0; i < 3; i++) {
    delete[] qI[i];
    delete[] B[i];
  }
  delete[] qI;
  delete[] B;

}


void calc_symmetric_matrix_eigenvecs(double** A, double* lam, double** V) {

  // if (myrank == 0) {
  //   printf("in calc_symmetric_matric_eigenvecs\n");
  //   fflush(stdout);
  // }

  calc_eigenvals( A , lam ) ;

  double ai[3][3], aj[3][3] , I[3][3] , P[3][3] , nm ;
  int i, j, k ;
  
  for ( i=0 ; i<3 ; i++ )  {
    for ( j=0 ; j<3 ; j++ ) 
      I[i][j] = 0.0 ;
    I[i][i] = 1.0 ;
  }

  // First eigenvector //
  for ( i=0 ; i<3 ; i++ ) 
    for ( j=0 ; j<3 ; j++ ) {
      P[i][j] = 0.0 ;
      ai[i][j] = A[i][j] - lam[1] * I[i][j] ;
      aj[i][j] = A[i][j] - lam[2] * I[i][j] ;
    }

  for ( i=0 ; i<3 ; i++ ) 
    for ( j=0 ; j<3 ; j++ ) 
      for ( k=0 ; k<3 ; k++ ) 
        P[i][j] += ai[i][k] * aj[k][j] ;

  nm = 0.0 ;
  for ( i=0; i<3 ; i++ )
    nm += P[0][i]*P[0][i] ;

  nm = sqrt( nm ) ;
  
  if ( P[0][2] < 0.0 )
    nm *= -1.0 ;

  for ( i=0 ; i<3 ; i++ )
    V[0][i] = P[0][i] / nm ;


  // if (myrank == 0) {
  //   printf("after first eigenvector\n");
  //   fflush(stdout);
  // }

  // Second eigenvector //
  for ( i=0 ; i<3 ; i++ ) 
    for ( j=0 ; j<3 ; j++ ) {
      P[i][j] = 0.0 ;
      ai[i][j] = A[i][j] - lam[0] * I[i][j] ;
      aj[i][j] = A[i][j] - lam[2] * I[i][j] ;
    }

  for ( i=0 ; i<3 ; i++ ) 
    for ( j=0 ; j<3 ; j++ ) 
      for ( k=0 ; k<3 ; k++ ) 
        P[i][j] += ai[i][k] * aj[k][j] ;

  nm = 0.0 ;
  for ( i=0; i<3 ; i++ )
    nm += P[0][i]*P[0][i] ;

  nm = sqrt( nm ) ;
  
  if ( P[0][2] < 0.0 )
    nm *= -1.0 ;

  for ( i=0 ; i<3 ; i++ )
    V[1][i] = P[0][i] / nm ;


  // if (myrank == 0) {
  //   printf("after second eigenvector\n");
  //   fflush(stdout);
  // }

  // Third eigenvector //
  for ( i=0 ; i<3 ; i++ ) 
    for ( j=0 ; j<3 ; j++ ) {
      P[i][j] = 0.0 ;
      ai[i][j] = A[i][j] - lam[0] * I[i][j] ;
      aj[i][j] = A[i][j] - lam[1] * I[i][j] ;
    }

  for ( i=0 ; i<3 ; i++ ) 
    for ( j=0 ; j<3 ; j++ ) 
      for ( k=0 ; k<3 ; k++ ) 
        P[i][j] += ai[i][k] * aj[k][j] ;

  nm = 0.0 ;
  for ( i=0; i<3 ; i++ )
    nm += P[0][i]*P[0][i] ;

  nm = sqrt( nm ) ;
  
  if ( P[0][2] < 0.0 )
    nm *= -1.0 ;

  for ( i=0 ; i<3 ; i++ )
    V[2][i] = P[0][i] / nm ;

  // if (myrank == 0) {
  //   printf("after third eigenvector\n");
  //   fflush(stdout);
  // }


}

