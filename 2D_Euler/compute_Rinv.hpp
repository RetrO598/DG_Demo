#pragma once
#include "global.hpp"
#include <cmath>
#include <vector>
inline void compute_Rinv(std::vector<std::vector<double>> &Rmat,
                         std::vector<std::vector<double>> &Rinv,
                         const double &rho, const double &u, const double &v,
                         const double &E, const double &n1, const double &n2) {
  double norm = std::sqrt(n1 * n1 + n2 * n2);
  double nf1, nf2;
  nf1 = n1 / norm;
  nf2 = n2 / norm;

  double pr = (E - 0.5 * rho * (u * u + v * v)) * (gasGamma - 1.0);
  double c = std::sqrt(gasGamma * pr / rho);
  double eH = (E + pr) / rho;
  double unf = u * nf1 + v * nf2;
  double ek = 0.5 * (u * u + v * v);

  Rmat[0][0] = 1.0;
  Rmat[1][0] = u - c * nf1;
  Rmat[2][0] = v - c * nf2;
  Rmat[3][0] = eH - c * unf;

  Rmat[0][1] = 1.0;
  Rmat[1][1] = u;
  Rmat[2][1] = v;
  Rmat[3][1] = ek;

  Rmat[0][2] = 1.0;
  Rmat[1][2] = u + c * nf1;
  Rmat[2][2] = v + c * nf2;
  Rmat[3][2] = eH + c * unf;

  Rmat[0][3] = 0.0;
  Rmat[1][3] = nf2;
  Rmat[2][3] = -nf1;
  Rmat[3][3] = u * nf2 - v * nf1;

  Rinv[0][0] = ((gasGamma - 1.0) * ek + c * unf) * 0.5 / (c * c);
  Rinv[1][0] = (c * c - (gasGamma - 1.0) * ek) / (c * c);
  Rinv[2][0] = ((gasGamma - 1.0) * ek - c * unf) * 0.5 / (c * c);
  Rinv[3][0] = v * nf1 - u * nf2;

  Rinv[0][1] = ((1.0 - gasGamma) * u - c * nf1) * 0.5 / (c * c);
  Rinv[1][1] = (gasGamma - 1.0) * u / (c * c);
  Rinv[2][1] = ((1.0 - gasGamma) * u + c * nf1) * 0.5 / (c * c);
  Rinv[3][1] = nf2;

  Rinv[0][2] = ((1.0 - gasGamma) * v - c * nf2) * 0.5 / (c * c);
  Rinv[1][2] = (gasGamma - 1.0) * v / (c * c);
  Rinv[2][2] = ((1.0 - gasGamma) * v + c * nf2) * 0.5 / (c * c);
  Rinv[3][2] = -nf1;

  Rinv[0][3] = (gasGamma - 1.0) * 0.5 / (c * c);
  Rinv[1][3] = (1.0 - gasGamma) / (c * c);
  Rinv[2][3] = (gasGamma - 1.0) * 0.5 / (c * c);
  Rinv[3][3] = 0.0;
}