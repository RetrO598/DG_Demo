#pragma once
#include "global.hpp"
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <xtensor.hpp>
inline void eigenvalueMm(double &Amax, double &Amin, const double &rho,
                         const double &rhou, const double &rhov,
                         const double &E, const double &n1, const double &n2) {
  double u = rhou / rho;
  double v = rhov / rho;
  double un = u * n1 + v * n2;
  double p = (gasGamma - 1.0) * (E - 0.5 * rho * (u * u + v * v));
  double c = std::sqrt(std::abs(gasGamma * p / rho));

  Amax = un + c;
  Amin = un - c;
}

inline double calculate_dt(const xt::xarray<double> &uh) {
  double alphax = 0.0;
  double alphay = 0.0;
  double alpha1 = 0.0;
  double alpha2 = 0.0;

  for (size_t j = 0; j < NY + 2; ++j) {
    for (size_t i = 0; i < NX + 2; ++i) {
      eigenvalueMm(alpha1, alpha2, uh(j, i, 0, 0), uh(j, i, 0, 1),
                   uh(j, i, 0, 2), uh(j, i, 0, 3), 1.0, 0.0);
      if (std::abs(alpha1) > alphax || std::abs(alpha2) > alphax) {
        alphax = std::max(std::abs(alpha1), std::abs(alpha2));
      }

      eigenvalueMm(alpha1, alpha2, uh(j, i, 0, 0), uh(j, i, 0, 1),
                   uh(j, i, 0, 2), uh(j, i, 0, 3), 0.0, 1.0);
      if (std::abs(alpha1) > alphay || std::abs(alpha2) > alphay) {
        alphay = std::max(std::abs(alpha1), std::abs(alpha2));
      }
    }
  }

  return CFL / (alphax / hx + alphay / hy);
}