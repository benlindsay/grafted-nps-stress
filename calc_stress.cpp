// calc_stress.cpp
//
// Copyright (c) 2016 Ben Lindsay <benjlindsay@gmail.com>

#include "globals.h"

void field_gradient_2(complex<double>*, complex<double>*, int);
void write_data(char*, complex<double>*);

void calc_dH_dL_diblock(complex<double> dH_dL_diblock[Dim]) {
  double b2 = 6.0 / (N - 1.0);
  for (int d = 0; d < Dim; d++) {
    complex<double> bond_factor = nD * b2 / (3.0 * Qd * L[d] * V);
    complex<double> smear_factor = a_squared / L[d];
    for (int i = 0; i < ML; i++) {
      diblock_stress[d][i] = 0;
    }
    for (int j = 0; j < N-1; j++) {
      field_gradient_2(qd[j], tmp, d);
      for (int i = 0; i < ML; i++) {
        tmp[i] *= qddag[N-j-1][i];
        if (j < Nda) {
          tmp[i] *= exp(smwa[i]);
        } else {
          tmp[i] *= exp(smwb[i]);
        }
        diblock_stress[d][i] += tmp[i];
      }
    }
    for (int i = 0; i < ML; i++) {
      diblock_stress[d][i] *= bond_factor;
    }
    field_gradient_2(smwa, tmp, d);
    field_gradient_2(smwb, tmp2, d);
    for (int i = 0; i < ML; i++) {
      tmp[i] *= rhoda[i];
      tmp2[i] *= rhodb[i];
      tmp[i] += tmp2[i];
      tmp[i] *= smear_factor;
      diblock_stress[d][i] += tmp[i];
    }
    dH_dL_diblock[d] = integ_trapPBC(diblock_stress[d]);
  }
}

void calc_dH_dL_diblock_other_way(complex<double> dH_dL_diblock[Dim]) {
  double b2 = 6.0 / (N - 1.0);
  for (int d = 0; d < Dim; d++) {
    dH_dL_diblock[d] = 0;
    complex<double> factor = nD * b2 / (3.0 * Qd * L[d] * V);
    for (int j = 1; j < N; j++) {
      field_gradient_2(qddag[N-j-1], tmp, d);
      for (int i = 0; i < ML; i++) {
        tmp[i] *= qd[j][i];
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
  // calc_dH_dL_diblock_other_way(dH_dL_diblock);
  // printf("dHdL_x other  way: %lf\n", dH_dL_diblock[0]);
  calc_dH_dL_diblock(dH_dL_diblock);
  // printf("dHdL_x normal way: %lf\n", dH_dL_diblock[0]);
  complex<double> dH_dL_grafts[Dim];
  if (sigma > 0) {
    calc_dH_dL_grafts(dH_dL_grafts);
  }
  for (int d = 0; d < Dim; d++) {
    stress_diblock[d] = L[d] / V * dH_dL_diblock[d];
    stress_grafts[d] = L[d] / V * dH_dL_grafts[d];
  }
}
