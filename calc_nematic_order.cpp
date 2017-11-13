// calc_nematic_order.cpp
//
// Copyright (c) 2017 Ben Lindsay <benjlindsay@gmail.com>

#include "globals.h"

void field_gradient_2(complex<double>*, complex<double>*, int, int);
void field_gradient_2_5pt(complex<double>*, complex<double>*, int);

void calc_one_nematic_term(complex<double> *nematic_order_element,
                           int dir1, int dir2) {
  double b2 = 6.0 / (N - 1.0);
  // complex<double> factor = V * b2 / (36.0 * N * Qd);
  complex<double> factor = b2 / 36.0;
  for (int i = 0; i < ML; i++) {
    nematic_order_element[i] = tmp[i] = tmp2[i] = 0.0;
  }
  for (int n = 0; n < N - 1; n++) {
    field_gradient_2(qddag[N-2-n], tmp, dir1, dir2);
    for (int i = 0; i < ML; i++) {
      nematic_order_element[i] += qd[n][i] * tmp[i];
    }
    field_gradient_2(qd[n], tmp, dir1, dir2);
    for (int i = 0; i < ML; i++) {
      nematic_order_element[i] += tmp[i] * qddag[N-2-n][i];
    }
    field_gradient(qd[n], tmp, dir1);
    field_gradient(qddag[N-2-n], tmp2, dir2);
    for (int i = 0; i < ML; i++) {
      nematic_order_element[i] += - tmp[i] * tmp2[i];
    }
    field_gradient(qd[n], tmp, dir2);
    field_gradient(qddag[N-2-n], tmp2, dir1);
    for (int i = 0; i < ML; i++) {
      nematic_order_element[i] += - tmp[i] * tmp2[i];
    }
    if (dir1 == dir2) {
      for (int d = 0; d < Dim; d++) {
        field_gradient_2(qddag[N-2-n], tmp, d, d);
        for (int i = 0; i < ML; i++) {
          nematic_order_element[i] += -1.0/3.0 * qd[n][i] * tmp[i];
        }
        field_gradient_2(qd[n], tmp, d, d);
        for (int i = 0; i < ML; i++) {
          nematic_order_element[i] += -1.0/3.0 * tmp[i] * qddag[N-2-n][i];
        }
        field_gradient(qd[n], tmp, d);
        field_gradient(qddag[N-2-n], tmp2, d);
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

// just bond stress for now...
void calc_nematic_order(complex<double> **nematic_order) {

  int i_nematic = 0;
  for (int dir1 = 0; dir1 < Dim - 1; dir1++) {
    for (int dir2 = dir1+1; dir2 < Dim; dir2++) {
      calc_one_nematic_term(nematic_order[i_nematic], dir1, dir2);
      i_nematic++;
    }
  }
}
