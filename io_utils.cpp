#include "globals.h"
void accumulate_all_averages(void);
void read_one_resume_file(FILE*, complex<double>*);
void write_avg_data(char*, complex<double>*);
void write_avg_kdata(char* , complex<double>*);
void write_kdata(char*,complex<double>*);
void write_avg_data_bin(char*, complex<double>*);
void write_data_bin(char*, complex<double>*);
void write_data(char*, complex<double>*);

void write_outputs() {
  int i , j;

  if (nD > 0.0) {
    fft_fwd_wrapper(rhoda, tmp);
    for (i=0; i<ML; i++) 
      tmp[i] *= hhat[i];
    fft_bck_wrapper(tmp, tmp);
    write_data_bin("rhoda", tmp);

    write_data_bin("rhoda_c", rhoda);
    
    fft_fwd_wrapper(rhodb, tmp);
    for (i=0; i<ML; i++) 
      tmp[i] *= hhat[i];
    fft_bck_wrapper(tmp, tmp);
    write_data_bin("rhodb", tmp);
    
    write_data_bin("rhodb_c", rhodb);
  }

  if (nAH > 0.0) {
    fft_fwd_wrapper(rhoha, tmp);
    for (i=0; i<ML; i++) 
      tmp[i] *= hhat[i];
    fft_bck_wrapper(tmp, tmp);
    write_data_bin("rhoha", tmp);
    write_data_bin("rhoha_c", rhoha);
  }

  if (n_exp_nr > 0 && iter == 1) {
    write_data_bin("rho_exp_nr", rho_exp_nr);
  }

  if (do_fld_np) {
    write_data_bin("rho_fld_np", rho_fld_np);
    write_data_bin("rho_fld_np_c", rho_fld_np_c);
  }

  if (do_CL && iter >= sample_freq) {
    int frame;
    char nm[20];
    frame = iter / print_freq;
    if (nD > 0.0) {
      write_avg_data_bin("avg_rhoda", avg_rhoda);
      write_avg_data_bin("avg_rhodb", avg_rhodb);
    }
        
    if (nAH > 0.0) {
      write_avg_data_bin("avg_rhoha", avg_rhoha);
    }
  }

  write_data_bin("wpl", wpl);
  write_data_bin("wabp", wabp);
  write_data_bin("wabm", wabm);
  
} // end write_outputs

// If running in parallel, this routine searches each file for
// the starting spatial position for this specific processor, 
// then it reads in ML data points.
void read_one_resume_file(FILE *inp, complex<double> *w ) {
  int i,j, nn[Dim];

  double dr, di, bgn[Dim], dm[Dim];

#ifdef PAR
  unstack(unstack_stack(0), nn);
  for (i=0; i<Dim; i++)
    bgn[i] = dx[i] * double(nn[i]);
#endif

  for (i=0; i<M; i++) {
    for (j=0; j<Dim; j++)
      fscanf(inp,"%lf ", &dm[j]);
    fscanf(inp, "%lf %lf\n", &dr, &di);

#ifndef PAR
    w[i] = dr + I*di;
#else
    double diff = 0.0;

    for (j=0; j<Dim; j++)
      diff += fabs(bgn[j] - dm[j]);

    if (diff < 0.001) {
      w[0] = dr + I*di;
      printf("Processor %d reading from line %d\n", myrank, i+1);
      break;
    }
#endif

  }//for (i=0; i<M[0];

#ifdef PAR
  for (i=1; i<ML; i++) {
    for (j=0; j<Dim; j++)
      fscanf(inp,"%lf ", &dm[j]);
    fscanf(inp, "%lf %lf\n", &dr, &di);
    w[i] = dr + I * di;
  }
#endif

} // end read_one_resume_file

void read_resume_files() {
  FILE *inp;
  char nm[20];
  int iter_file_flag = 0, i , j;

  inp = fopen("wpl.res","r");
  if (inp!=NULL) {
    read_one_resume_file(inp, wpl);
    fclose(inp);
    printf("Read wpl.res\n");
  }

  inp = fopen("wabp.res","r");
  if (inp!=NULL) {
    read_one_resume_file(inp, wabp);
    fclose(inp);
    printf("Read wabp.res\n");
  }

  inp = fopen("wabm.res","r");
  if (inp!=NULL) {
    read_one_resume_file(inp, wabm);
    fclose(inp);
    printf("Read wabm.res\n");
  }
} // end read_resume_files

void write_fft_cpx(char* nm, fftw_complex *dt) {

  int i,j, nn[Dim];
  FILE *otp;
  otp = fopen(nm, "w");

  for (i=0; i<ML; i++) {
    unstack( unstack_stack(i) , nn );
    fprintf(otp,"%lf ", double(nn[0])*dx[0] )  ;
    for (j=1; j<Dim; j++)
      fprintf(otp,"%lf ", double(nn[j])*dx[j]);
    fprintf(otp,"%1.16e %1.16e\n", dt[i][0], dt[i][1] );

    if (Dim==2 && nn[0]==Nx[0]-1)
      fprintf(otp,"\n");

  }

  fclose(otp);
} // end write_fft_cpx

void write_kdata(char* nmi, complex<double> *dt) {

  int i,j, nn[Dim] ;
  double k2, kv[Dim];
  FILE *otp;
  char nm[20];
#ifndef PAR
  sprintf( nm , "%s.dat", nmi);
#else
  sprintf( nm , "%s.p%d.dat" , nmi, myrank);
#endif
  otp = fopen(nm, "w");

  for (i=0; i<ML; i++) {
    k2 = get_k(i, kv);
    fprintf(otp, "%lf ", kv[0]);
    for (j=1; j<Dim; j++)
      fprintf(otp, "%lf ", kv[j]);
    fprintf(otp,"%1.16e %1.16e %1.16e %1.16e\n", abs(dt[i]), sqrt(k2),
            real(dt[i]), imag(dt[i]) );

    if (Dim==2 && nn[0]==Nx[0]-1)
      fprintf(otp, "\n");
  }

  fclose(otp);
} // end write_kdata


void write_avg_kdata(char* nmi, complex<double> *dt) {

  int i,j, nn[Dim] ;
  double k2, kv[Dim];
  FILE *otp;
  char nm[20];
#ifndef PAR
  sprintf( nm , "%s.dat", nmi);
#else
  sprintf( nm , "%s.p%d.dat" , nmi, myrank);
#endif
  otp = fopen(nm, "w");

  for (i=0; i<ML; i++) {
    k2 = get_k( i , kv ) ;
    fprintf(otp,"%lf ", kv[0] ) ;
    for (j=1; j<Dim; j++)
      fprintf(otp,"%lf ", kv[j] ) ;
    fprintf(otp,"%1.16e %1.16e %1.16e %1.16e\n", abs(dt[i])/n_samples, sqrt(k2),
        real(dt[i])/n_samples, imag(dt[i])/n_samples );

    if (Dim==2 && nn[0]==Nx[0]-1)
      fprintf(otp,"\n");

  }

  fclose(otp);
} // end write_avg_kdata


void write_avg_data(char* nmi, complex<double> *dt) {

  int i,j, nn[Dim];
  FILE *otp;
  char nm[20];
#ifndef PAR
  sprintf( nm , "%s.dat", nmi);
#else
  sprintf( nm , "%s.p%d.dat" , nmi, myrank);
#endif
  otp = fopen(nm, "w");

  for (i=0; i<ML; i++) {
    unstack( unstack_stack(i) , nn );
    fprintf(otp,"%lf ", double(nn[0])*dx[0]  ) ;
    for (j=1; j<Dim; j++)
      fprintf(otp,"%lf ", double(nn[j])*dx[j]);
    fprintf(otp,"%1.16e %1.16e\n", real(dt[i])/n_samples, imag(dt[i])/n_samples );

    if (Dim==2 && nn[0]==Nx[0]-1)
      fprintf(otp,"\n");

  }

  fclose(otp);
} // end write_avg_data

void write_avg_data_bin(char* nmi, complex<double> *dt) {

  int i,j, nn[Dim];
  FILE *otp;
  char nm[20];
  
  sprintf( nm , "%s.p%d.bin" , nmi, myrank);
  otp = fopen(nm, "wb");

  // The main cpu writes the number of processors, grid spacing, Nx, Ny, Nz
  if ( myrank == 0 ) {
    int dm = Dim ;
    fwrite( &dm , sizeof( int ) , 1 , otp ) ;
    fwrite( Nx , sizeof( int ) , dm , otp ) ;
    fwrite( L , sizeof( double ) , dm , otp ) ;
    fwrite( &nprocs , sizeof( int ) , 1 , otp ) ;
  }

  fwrite( &ML , sizeof( int ) , 1 , otp ) ;
  for ( i=0 ; i<ML ; i++ )
    tmp[i] = dt[i] / n_samples ;
  fwrite( tmp , sizeof( complex<double> ) , ML , otp ) ;

  fclose(otp);
} // end write_avg_data_bin

void write_data_bin(char* nmi, complex<double> *dt) {

  int i,j, nn[Dim];
  FILE *otp;
  char nm[20];
  
  sprintf(nm, "%s.p%d.bin", nmi, myrank);
  otp = fopen(nm, "wb");

  // The main cpu writes the number of processors, grid spacing, Nx, Ny, Nz
  if (myrank == 0) {
    int dm = Dim;
    fwrite(&dm, sizeof(int), 1, otp);
    fwrite(Nx, sizeof(int), dm, otp);
    fwrite(L, sizeof(double), dm, otp);
    fwrite(&nprocs, sizeof(int), 1, otp);
  }

  fwrite(&ML, sizeof(int), 1, otp);
  fwrite(dt, sizeof(complex<double>), ML, otp);

  fclose(otp);
} // end write_data_bin

void write_data(char* nmi, complex<double> *dt) {

  int i,j, nn[Dim];
  FILE *otp;
  char nm[20];
#ifndef PAR
  sprintf(nm, "%s.dat", nmi);
#else
  sprintf(nm, "%s.p%d.dat", nmi, myrank);
#endif
  otp = fopen(nm, "w");

  for (i=0; i<ML; i++) {
    unstack(unstack_stack(i), nn);
    fprintf(otp, "%lf ", double(nn[0])*dx[0] );
    for (j=1; j<Dim; j++)
      fprintf(otp, "%lf ", double(nn[j])*dx[j] );
    fprintf(otp, "%1.16e %1.16e\n", real(dt[i]), imag(dt[i]) );

    if (Dim==2 && nn[0]==Nx[0]-1)
      fprintf(otp,"\n");
  }

  fclose(otp);
} // end write_data


void append_data(char* nmi, complex<double> *dt) {

  int i,j, nn[Dim];
  FILE *otp;
  char nm[20];
#ifndef PAR
  sprintf(nm, "%s.dat", nmi);
#else
  sprintf(nm, "%s.p%d.dat", nmi, myrank);
#endif

  otp = fopen(nm, "a");

  for (i=0; i<ML; i++) {
    unstack(unstack_stack(i), nn);
    for (j=0; j<Dim; j++)
      fprintf(otp, "%lf ", double(nn[j])*dx[j] );
    fprintf(otp, "%1.12e %1.12e\n", real(dt[i]), imag(dt[i]) );

    if (Dim==2 && nn[0]==Nx[0]-1)
      fprintf(otp,"\n");
  }

  fclose(otp);
} // end append_data
