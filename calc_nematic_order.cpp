// calc_nematic_order.cpp
//
// Copyright (c) 2017 Ben Lindsay <benjlindsay@gmail.com>

#include "globals.h"

void field_gradient_2(complex<double>*, complex<double>*, int, int);
void field_gradient_2_5pt(complex<double>*, complex<double>*, int);
void calc_symmetric_matrix_eigenvecs(double**, double*, double**);

void calc_one_nematic_term(
      complex<double> *nematic_order_element,
      int dir1,
      int dir2,
      complex<double> Q_local,
      complex<double> **q_local,
      complex<double> **qdag_local,
      int N_local) {
  // No matter what component we're looking at right now, b2 (kuhn length
  // squared) depends on diblock chain length
  double b2 = 6.0 / (Nda + Ndb - 1.0);
  // complex<double> factor = V * b2 / (36.0 * N * Qd);
  complex<double> factor = b2 / (36.0 * Q_local * (N_local - 1.0));
  for (int i = 0; i < ML; i++) {
    nematic_order_element[i] = tmp[i] = tmp2[i] = 0.0;
  }
  for (int n = 0; n < N_local - 1; n++) {
    field_gradient_2(qdag_local[N_local-2-n], tmp, dir1, dir2);
    for (int i = 0; i < ML; i++) {
      nematic_order_element[i] += q_local[n][i] * tmp[i];
    }
    field_gradient_2(q_local[n], tmp, dir1, dir2);
    for (int i = 0; i < ML; i++) {
      nematic_order_element[i] += tmp[i] * qdag_local[N_local-2-n][i];
    }
    field_gradient(q_local[n], tmp, dir1);
    field_gradient(qdag_local[N_local-2-n], tmp2, dir2);
    for (int i = 0; i < ML; i++) {
      nematic_order_element[i] += - tmp[i] * tmp2[i];
    }
    field_gradient(q_local[n], tmp, dir2);
    field_gradient(qdag_local[N_local-2-n], tmp2, dir1);
    for (int i = 0; i < ML; i++) {
      nematic_order_element[i] += - tmp[i] * tmp2[i];
    }
    if (dir1 == dir2) {
      for (int d = 0; d < Dim; d++) {
        field_gradient_2(qdag_local[N_local-2-n], tmp, d, d);
        for (int i = 0; i < ML; i++) {
          nematic_order_element[i] += -1.0/3.0 * q_local[n][i] * tmp[i];
        }
        field_gradient_2(q_local[n], tmp, d, d);
        for (int i = 0; i < ML; i++) {
          nematic_order_element[i] += -1.0/3.0 * tmp[i]
                                      * qdag_local[N_local-2-n][i];
        }
        field_gradient(q_local[n], tmp, d);
        field_gradient(qdag_local[N_local-2-n], tmp2, d);
        for (int i = 0; i < ML; i++) {
          nematic_order_element[i] += 2.0/3.0 * tmp[i] * tmp2[i];
        }
      }
    }
  }
  for (int i = 0; i < ML; i++) {
    nematic_order_element[i] *= factor;
  }
}

void calc_eigenvecs(complex<double> **nematic_order) {
    double **A = new double*[Dim];
    double **V = new double*[Dim];
    double *lam = new double[Dim];
    int i_nematic;
    for (int i = 0; i < ML; i++) {
      for (int d = 0; d < Dim; d++) {
        A[d] = new double[Dim];
        V[d] = new double[Dim];
      }
      i_nematic = 0;
      for (int d1 = 0; d1 < Dim; d1++) {
        for (int d2 = d1; d2 < Dim; d2++) {
          A[d1][d2] = real(nematic_order[i_nematic][i]);
          A[d2][d1] = A[d1][d2];
          i_nematic++;
        }
      }
      calc_symmetric_matrix_eigenvecs(A, lam, V);
      for (int d = 0; d < Dim; d++) {
        eigenvec[d][i] = V[0][d];
      }
    }
    for (int d = 0; d < Dim; d++) {
      delete[] A[d];
      delete[] V[d];
    }
    delete[] A;
    delete[] V;
    delete[] lam;
}

// just bond stress for now...
void calc_nematic_order(
      complex<double> **nematic_order,
      complex<double> Q_local,
      complex<double> **q_local,
      complex<double> **qdag_local,
      int N_local) {
  int i_nematic = 0;
  for (int dir1 = 0; dir1 < Dim; dir1++) {
    for (int dir2 = dir1; dir2 < Dim; dir2++) {
      calc_one_nematic_term(nematic_order[i_nematic], dir1, dir2, Q_local,
                            q_local, qdag_local, N_local);
      i_nematic++;
    }
  }
  if (nematic_order_output_mode == 1 && Dim == 3) {
    if (myrank == 0) {
      printf("about to enter calc_eigenvecs()\n");
      fflush(stdout);
    }
    calc_eigenvecs(nematic_order);
    if (myrank == 0) {
      printf("done with calc_eigenvecs()\n");
      fflush(stdout);
    }
  }
}
