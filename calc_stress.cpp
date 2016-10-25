// calc_stress.cpp
//
// Copyright (c) 2016 Ben Lindsay <benjlindsay@gmail.com>

#include "globals.h"

void field_gradient_2(complex<double>*, complex<double>*, int);

void calc_dH_dL_diblock(complex<double> dH_dL_diblock[Dim]) {
  double b2 = 6.0 / (N - 1.0);
  for (int d = 0; d < Dim; d++) {
    dH_dL_diblock[d] = 0;
    complex<double> factor = nD * b2 / (3.0 * Qd * L[d] * V);
    for (int j = 0; j < N-1; j++) {
      field_gradient_2(qd[j], tmp, d);
      for (int i = 0; i < ML; i++) {
        tmp[i] *= qddag[N-j-1][i];
        if (j < Nda) {
          tmp[i] *= exp(smwa[i]);
        } else {
          tmp[i] *= exp(smwb[i]);
        }
      }
      dH_dL_diblock[d] += integ_trapPBC(tmp);
    }
    dH_dL_diblock[d] *= factor;
  }
}

void calc_dH_dL_grafts(complex<double> dH_dL_grafts[Dim]) {
  double b2 = 6.0 / (N - 1.0);
  double n_GA = ng_per_np * (n_exp_nr);
  for (int d = 0; d < Dim; d++) {
    dH_dL_grafts[d] = 0;
    complex<double> factor = n_GA * b2 / (3.0 * L[d]);
    for (int j = 0; j < N-1; j++) {
      field_gradient_2(qg[j], tmp, d);
      for (int i = 0; i < ML; i++) {
        tmp[i] *= qgdag_exp[N-j-1][i];
        tmp[i] *= exp(smwa[i]);
      }
      dH_dL_grafts[d] += integ_trapPBC(tmp);
    }
    dH_dL_grafts[d] *= factor;
  }
}

// just bond stress for now...
void calc_stress(complex<double> stress_diblock[Dim],
                 complex<double> stress_grafts[Dim]) {
  complex<double> dH_dL_diblock[Dim];
  calc_dH_dL_diblock(dH_dL_diblock);
  complex<double> dH_dL_grafts[Dim];
  if (sigma > 0) {
    calc_dH_dL_grafts(dH_dL_grafts);
  }
  for (int d = 0; d < Dim; d++) {
    stress_diblock[d] = L[d] / V * dH_dL_diblock[d];
    stress_grafts[d] = L[d] / V * dH_dL_grafts[d];
  }
}
