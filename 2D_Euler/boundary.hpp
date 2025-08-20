#pragma once
#include "global.hpp"
#include <cstddef>
#include <xtensor.hpp>
inline void set_bc(xt::xarray<double> &uh) {
  for (size_t j = 0; j < NY + 2; ++j) {
    for (size_t k = 0; k < dimPk; ++k) {
      for (size_t n = 0; n < NumEq; ++n) {
        uh(j, 0, k, n) = uh(j, 1, k, n);
        uh(j, NX + 1, k, n) = uh(j, NX, k, n);
      }
    }
  }

  for (size_t i = 0; i < NX + 2; ++i) {
    for (size_t k = 0; k < dimPk; ++k) {
      for (size_t n = 0; n < NumEq; ++n) {
        uh(0, i, k, n) = uh(1, i, k, n);
        uh(NY + 1, i, k, n) = uh(NY, i, k, n);
      }
    }
  }
}