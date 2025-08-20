#pragma once
#include "global.hpp"
#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include <vector>

inline void minmod(const int &direction, const std::vector<double> &deltaU,
                   const std::vector<double> &deltaUL,
                   const std::vector<double> &deltaUR,
                   std::vector<double> &deltaUmod) {

  double hd = (direction == 1) ? hx : hy;
  for (size_t i = 0; i < NumEq; ++i) {
    if (std::abs(deltaU[i]) <= M * hd * hd) {
      deltaUmod[i] = deltaU[i];
    } else {
      int a = (deltaU[i] >= 0.0) ? 1 : -1;
      int b = (deltaUR[i] >= 0.0) ? 1 : -1;
      int c = (deltaUL[i] >= 0.0) ? 1 : -1;
      int s = (a + b + c) / 3;
      if (std::abs(s) == 1) {
        deltaUmod[i] = s * std::min({std::abs(deltaU[i]), std::abs(deltaUL[i]),
                                     std::abs(deltaUR[i])});
      } else {
        deltaUmod[i] = 0.0;
      }
    }
  }
}
