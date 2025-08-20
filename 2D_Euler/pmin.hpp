#pragma once
#include "basis.hpp"
#include "global.hpp"
#include "pressure.hpp"
#include <cstddef>
#include <xtensor.hpp>
#include <xtensor/core/xtensor_forward.hpp>
inline double calculate_pmin(const xt::xarray<double> &uh) {
  double pmin = 1000000.0;

  for (size_t j = 1; j < NY + 1; ++j) {
    for (size_t i = 1; i < NX + 1; ++i) {
      xt::xarray<double> uhGLL({NumGLP, NumGLP, NumEq, 3}, 0.0);
      for (size_t n = 0; n < NumEq; ++n) {
        for (size_t d = 0; d < dimPk; ++d) {
          for (size_t j1 = 0; j1 < NumGLP; ++j1) {
            for (size_t i1 = 0; i1 < NumGLP; ++i1) {
              uhGLL(j1, i1, n, 0) +=
                  uh(j, i, d, n) * Poly(lambdaL[i1], lambda[j1], d);
              uhGLL(j1, i1, n, 1) +=
                  uh(j, i, d, n) * Poly(lambda[i1], lambdaL[j1], d);
              uhGLL(j1, i1, n, 2) +=
                  uh(j, i, d, n) * Poly(lambda[i1], lambda[j1], d);
            }
          }
        }
      }

      for (size_t j1 = 0; j1 < NumGLP; ++j1) {
        for (size_t i1 = 0; i1 < NumGLP; ++i1) {
          for (size_t d = 0; d < 3; ++d) {
            double p1 =
                pressure(uhGLL(j1, i1, 0, d), uhGLL(j1, i1, 1, d),
                         uhGLL(j1, i1, 2, d), uhGLL(j1, i1, 3, d), gasGamma);
            if (p1 < pmin) {
              pmin = p1;
            }
          }
        }
      }
    }
  }
  return pmin;
}