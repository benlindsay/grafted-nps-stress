#include "globals.h"
void random_particles( void ) ;
void allocate_particles( void ) ;

void read_input() {

  FILE *inp;
  inp = fopen("bcp.input", "r");

  if (inp==NULL) { 
    printf("Failed to open bcp.input!");
    exit(1);
  }

  char tt[80];
  double dm1, dm2;
  int di1, di2,i, j;

  if (myrank == 0) printf("---Reading bcp.input---\n");

  // Main Parameters
  fgets( tt, 80, inp );
  fscanf(inp,"%d %d", &Nda , &Ndb );     fgets(tt,80,inp);
  fscanf( inp , "%lf %d" , &phiH , &Nah) ; fgets( tt , 80 , inp ) ;
  fscanf( inp  , "%lf" , &a_smear ) ; fgets( tt , 80 ,  inp ) ;
  fscanf(inp,"%lf", &C );      fgets(tt,80,inp);
  
  fgets( tt ,80 , inp ) ;
  fgets( tt ,80 , inp ) ;

  // Interaction Parameters
  fscanf( inp , "%lf" ,&chiN ) ; fgets( tt , 80 , inp ) ;
  fscanf( inp , "%lf" ,&kappaN ) ; fgets( tt , 80 , inp ) ;

  fgets( tt ,80 , inp ) ;
  fgets( tt ,80 , inp ) ;

  // Initial condition flags
  fscanf( inp , "%d %lf %lf %lf" , &ic_flag[0] ,
          &ic_pre[0] , &ic_dir[0] , &ic_period[0] ) ; 
  fgets(tt , 80 , inp ) ;
  fscanf( inp , "%d %lf %lf %lf" , &ic_flag[1] ,
          &ic_pre[1] , &ic_dir[1] , &ic_period[1] ) ; 
  fgets(tt , 80 , inp ) ;
  fscanf(inp, "%d", &keep_fields); fgets( tt ,80 , inp ) ;

  fgets( tt ,80 , inp ) ;
  fgets( tt ,80 , inp ) ;

  // Numerical parameters
  for ( i=0 ; i<Dim ; i++ ) 
    fscanf( inp, "%d" , &Nx[i] );
  fgets( tt, 80 , inp );
  for ( i=0 ; i<Dim ; i++ ) 
    fscanf( inp, "%lf" , &L[i] );
  fgets( tt, 80 , inp );
  if ( myrank == 0 ) 
    for (i=0; i<Dim ; i++)  
      printf("Nx%d: %d Lx%d: %lf\n" , i, Nx[i] , i, L[i] ) ;

  // Brent's method parameters
  fscanf(inp, "%d", &do_brent); fgets(tt, 80, inp);
  fscanf(inp, "%lf %lf %lf %lf", &L_low, &L_high, &L_step, &brent_tol);
  fgets(tt, 80, inp);

  // Lambda/time step parameters
  fscanf(inp,"%lf", &lam_pl ); fgets(tt,80,inp);
  fscanf(inp,"%lf", &lam_mi ); fgets(tt,80,inp);

  // Iteration number parameters
  fscanf(inp,"%d",&itermax);     fgets(tt,80,inp);
  fscanf(inp,"%d",&print_freq);  fgets(tt,80,inp);
  fscanf(inp,"%d %d",&sample_freq , &sample_wait );  fgets(tt,80,inp);
  fscanf(inp,"%d", &update_scheme); fgets(tt,80,inp);
  fscanf( inp , "%d" , &do_CL ) ; fgets( tt , 80 , inp ) ;
  
  fgets( tt ,80 , inp ) ;
  fgets( tt ,80 , inp ) ;

  // Film parameters
  fscanf( inp, "%d" , &do_film) ; fgets( tt, 80, inp) ;
  fscanf( inp, "%lf %lf", &wallT, &wallXi ) ; fgets(tt,80,inp);
  fscanf( inp, "%lf %lf", &wall_lamA, &wall_lamB ) ; fgets( tt, 80, inp);

  fgets( tt ,80 , inp ) ;
  fgets( tt ,80 , inp ) ;
  
  // Explicit nanorod parameters
  fscanf( inp, "%d" , &n_exp_nr) ; fgets( tt, 80, inp) ;
  fscanf(inp, "%lf %lf %lf", &L_nr, &R_nr, &xi_nr); fgets(tt, 80, inp);
  fscanf( inp, "%lf %lf", &exp_nr_chiAPN, &exp_nr_chiBPN ); fgets(tt, 80, inp);
  for (j=0; j<2; j++) {
    for (i=0; i<Dim; i++)
      fscanf(inp, "%lf", &exp_nr_c[j][i]);
    fgets( tt, 80, inp);

    dm1=0;
    for (i=0; i<Dim; i++) {
      fscanf(inp, "%lf", &exp_nr_u[j][i]);
      // Get vector magnitude to normalize as it's scanned in
      dm1 += exp_nr_u[j][i]*exp_nr_u[j][i];
    }
    fgets(tt, 80, inp);
    dm2 = sqrt(dm1);

    // Normalize nanorod orientation vector
    for (i=0; i<Dim; i++)
      exp_nr_u[j][i] /= dm2;
  }
  fclose(inp);
  printf("---Reading of bcp.input complete---\n");
}
