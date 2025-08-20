#pragma once

#include "global.hpp"
#include <cstddef>
#include <cstdlib>
#include <xtensor.hpp>
inline double calculate_umax(const xt::xarray<double> &uh) {
  double umax = 0.0;
  for (size_t j = 1; j < NY + 1; ++j) {
    for (size_t i = 1; i < NX + 1; ++i) {
      if (std::abs(uh(j, i, 0, 0)) > umax) {
        umax = std::abs(uh(j, i, 0, 0));
      }
    }
  }
  return umax;
}