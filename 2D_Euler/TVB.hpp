#pragma once

#include "compute_Rinv.hpp"
#include "global.hpp"
#include "matmul.hpp"
#include "minmod.hpp"
#include <array>
#include <cstddef>
#include <cstdlib>
#include <fstream>
#include <vector>
#include <xtensor.hpp>
inline void TVB_limiter(xt::xarray<double> &uh) {

  auto uhmod = uh;
  std::vector<std::vector<double>> Rmat(4, {0.0, 0.0, 0.0, 0.0});
  std::vector<std::vector<double>> Rinv(4, {0.0, 0.0, 0.0, 0.0});

  std::vector<double> deltaUR(4, 0.0);
  std::vector<double> deltaUL(4, 0.0);
  std::vector<double> deltaUR1(4, 0.0);
  std::vector<double> deltaUL1(4, 0.0);
  std::vector<double> deltaUU1(4, 0.0);
  std::vector<double> deltaUD1(4, 0.0);
  std::vector<double> deltaUmod(4, 0.0);
  std::vector<double> deltaU(4, 0.0);
  std::vector<double> deltaUR1mod(4, 0.0);
  std::vector<double> deltaUL1mod(4, 0.0);
  std::vector<double> deltaUU1mod(4, 0.0);
  std::vector<double> deltaUD1mod(4, 0.0);

  for (size_t j = 1; j < NY + 1; ++j) {
    for (size_t i = 1; i < NX + 1; ++i) {
      std::array<int, 4> change{0, 0, 0, 0};
      double rho = uh(j, i, 0, 0);
      double u = uh(j, i, 0, 1) / rho;
      double v = uh(j, i, 0, 2) / rho;
      double E = uh(j, i, 0, 3);

      compute_Rinv(Rmat, Rinv, rho, u, v, E, 1.0, 0.0);

      for (size_t n = 0; n < NumEq; ++n) {
        deltaUR[n] = uh(j, i + 1, 0, n) - uh(j, i, 0, n);
        deltaUL[n] = uh(j, i, 0, n) - uh(j, i - 1, 0, n);
        deltaUR1[n] = uh(j, i, 1, n) + 2.0 / 3.0 * uh(j, i, 3, n);
        deltaUL1[n] = uh(j, i, 1, n) - 2.0 / 3.0 * uh(j, i, 3, n);
      }

      deltaUR = matMul(Rinv, deltaUR, 4);
      deltaUL = matMul(Rinv, deltaUL, 4);
      deltaU = matMul(Rinv, deltaUR1, 4);

      minmod(1, deltaU, deltaUL, deltaUR, deltaUmod);

      deltaUR1mod = matMul(Rmat, deltaUmod, 4);

      for (size_t n = 0; n < NumEq; ++n) {
        if (std::abs(deltaUR1mod[n] - deltaUR1[n]) > 1e-6) {
          change[n] = 1;
        }
      }

      deltaU = matMul(Rinv, deltaUL1, 4);

      minmod(1, deltaU, deltaUL, deltaUR, deltaUmod);

      deltaUL1mod = matMul(Rmat, deltaUmod, 4);

      for (size_t n = 0; n < NumEq; ++n) {
        if (std::abs(deltaUL1mod[n] - deltaUL1[n]) > 1e-6) {
          change[n] = 1;
        }
      }

      compute_Rinv(Rmat, Rinv, rho, u, v, E, 0.0, 1.0);
      for (size_t n = 0; n < NumEq; ++n) {
        deltaUR[n] = uh(j + 1, i, 0, n) - uh(j, i, 0, n);
        deltaUL[n] = uh(j, i, 0, n) - uh(j - 1, i, 0, n);
        deltaUU1[n] = uh(j, i, 2, n) + 2.0 / 3.0 * uh(j, i, 5, n);
        deltaUD1[n] = uh(j, i, 2, n) - 2.0 / 3.0 * uh(j, i, 5, n);
      }

      deltaUR = matMul(Rinv, deltaUR, 4);
      deltaUL = matMul(Rinv, deltaUL, 4);
      deltaU = matMul(Rinv, deltaUU1, 4);

      minmod(2, deltaU, deltaUL, deltaUR, deltaUmod);

      deltaUU1mod = matMul(Rmat, deltaUmod, 4);

      for (size_t n = 0; n < NumEq; ++n) {
        if (std::abs(deltaUU1mod[n] - deltaUU1[n]) > 1e-6) {
          change[n] = 1;
        }
      }

      deltaU = matMul(Rinv, deltaUD1, 4);

      minmod(2, deltaU, deltaUL, deltaUR, deltaUmod);
      deltaUD1mod = matMul(Rmat, deltaUmod, 4);

      for (size_t n = 0; n < NumEq; ++n) {
        if (std::abs(deltaUD1mod[n] - deltaUD1[n]) > 1e-6) {
          change[n] = 1;
        }
      }

      for (size_t n = 0; n < NumEq; ++n) {
        if (change[n] == 1) {
          uhmod(j, i, 4, n) = 0.0;
          uhmod(j, i, 1, n) = 0.5 * (deltaUR1mod[n] + deltaUL1mod[n]);
          uhmod(j, i, 2, n) = 0.5 * (deltaUU1mod[n] + deltaUD1mod[n]);
          uhmod(j, i, 3, n) = 3.0 / 4.0 * (deltaUR1mod[n] - deltaUL1mod[n]);
          uhmod(j, i, 5, n) = 3.0 / 4.0 * (deltaUU1mod[n] - deltaUD1mod[n]);
        }
      }
    }
  }

  uh = uhmod;
}