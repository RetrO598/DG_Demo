#pragma once

#include "HLL.hpp"
#include "basis.hpp"
#include "calculateDt.hpp"
#include "global.hpp"
#include "output.hpp"
#include <cstddef>
#include <xtensor.hpp>
#include <xtensor/core/xtensor_forward.hpp>

inline xt::xarray<double> dgOperator(const xt::xarray<double> &uh) {
  xt::xarray<double> du({NY + 2, NX + 2, dimPk, NumEq}, 0.0);
  xt::xarray<double> rhoM({NumGLP, NumGLP}, 0.0);
  xt::xarray<double> uM({NumGLP, NumGLP}, 0.0);
  xt::xarray<double> vM({NumGLP, NumGLP}, 0.0);
  xt::xarray<double> EM({NumGLP, NumGLP}, 0.0);
  xt::xarray<double> pM({NumGLP, NumGLP}, 0.0);
  xt::xarray<double> Fx({NY + 2, NX + 2, dimPk, NumEq}, 0.0);
  xt::xarray<double> Fy({NY + 2, NX + 2, dimPk, NumEq}, 0.0);
  for (size_t j = 1; j < NY + 1; ++j) {
    for (size_t i = 1; i < NX + 1; ++i) {
      xt::xarray<double> uGint({NumGLP, NumGLP, NumEq}, 0.0);
      for (size_t j1 = 0; j1 < NumGLP; ++j1) {
        for (size_t i1 = 0; i1 < NumGLP; ++i1) {
          for (size_t n = 0; n < NumEq; ++n) {
            for (size_t k = 0; k < dimPk; ++k) {
              uGint(j1, i1, n) +=
                  uh(j, i, k, n) * Poly(lambda[i1], lambda[j1], k);
            }
          }
        }
      }

      for (size_t j1 = 0; j1 < NumGLP; ++j1) {
        for (size_t i1 = 0; i1 < NumGLP; ++i1) {
          rhoM(j1, i1) = uGint(j1, i1, 0);
          uM(j1, i1) = uGint(j1, i1, 1) / rhoM(j1, i1);
          vM(j1, i1) = uGint(j1, i1, 2) / rhoM(j1, i1);
          EM(j1, i1) = uGint(j1, i1, 3);

          pM(j1, i1) =
              (gasGamma - 1.0) * (EM(j1, i1) - 0.5 * rhoM(j1, i1) *
                                                   (uM(j1, i1) * uM(j1, i1) +
                                                    vM(j1, i1) * vM(j1, i1)));

          Fx(j1, i1, 0) = uGint(j1, i1, 1);
          Fx(j1, i1, 1) = rhoM(j1, i1) * uM(j1, i1) * uM(j1, i1) + pM(j1, i1);
          Fx(j1, i1, 2) = rhoM(j1, i1) * uM(j1, i1) * vM(j1, i1);
          Fx(j1, i1, 3) = uM(j1, i1) * (EM(j1, i1) + pM(j1, i1));

          Fy(j1, i1, 0) = uGint(j1, i1, 2);
          Fy(j1, i1, 1) = rhoM(j1, i1) * uM(j1, i1) * vM(j1, i1);
          Fy(j1, i1, 2) = rhoM(j1, i1) * vM(j1, i1) * vM(j1, i1) + pM(j1, i1);
          Fy(j1, i1, 3) = vM(j1, i1) * (EM(j1, i1) + pM(j1, i1));
        }
      }

      for (size_t j1 = 0; j1 < NumGLP; ++j1) {
        for (size_t i1 = 0; i1 < NumGLP; ++i1) {
          for (size_t k = 1; k < dimPk; ++k) {
            for (size_t n = 0; n < NumEq; ++n) {
              du(j, i, k, n) +=
                  0.25 * weight[i1] * weight[j1] *
                  (Fx(j1, i1, n) * PolyX(lambda[i1], lambda[j1], k) +
                   Fy(j1, i1, n) * PolyY(lambda[i1], lambda[j1], k));
            }
          }
        }
      }
    }
  }

  xt::xarray<double> UL({NY, NX + 1, NumGLP, NumEq}, 0.0);
  xt::xarray<double> UR({NY, NX + 1, NumGLP, NumEq}, 0.0);
  xt::xarray<double> UU({NY + 1, NX, NumGLP, NumEq}, 0.0);
  xt::xarray<double> UD({NY + 1, NX, NumGLP, NumEq}, 0.0);
  xt::xarray<double> FL({NY, NX + 1, NumGLP, NumEq}, 0.0);
  xt::xarray<double> FR({NY, NX + 1, NumGLP, NumEq}, 0.0);
  xt::xarray<double> FU({NY + 1, NX, NumGLP, NumEq}, 0.0);
  xt::xarray<double> FD({NY + 1, NX, NumGLP, NumEq}, 0.0);

  for (size_t j = 0; j < NY; ++j) {
    for (size_t i = 0; i < NX + 1; ++i) {
      for (size_t k = 0; k < NumGLP; ++k) {
        for (size_t d = 0; d < dimPk; ++d) {
          UL(j, i, k, 0) += uh(j + 1, i, d, 0) * Poly(1.0, lambda[k], d);
          UL(j, i, k, 1) += uh(j + 1, i, d, 1) * Poly(1.0, lambda[k], d);
          UL(j, i, k, 2) += uh(j + 1, i, d, 2) * Poly(1.0, lambda[k], d);
          UL(j, i, k, 3) += uh(j + 1, i, d, 3) * Poly(1.0, lambda[k], d);
          UR(j, i, k, 0) += uh(j + 1, i + 1, d, 0) * Poly(-1.0, lambda[k], d);
          UR(j, i, k, 1) += uh(j + 1, i + 1, d, 1) * Poly(-1.0, lambda[k], d);
          UR(j, i, k, 2) += uh(j + 1, i + 1, d, 2) * Poly(-1.0, lambda[k], d);
          UR(j, i, k, 3) += uh(j + 1, i + 1, d, 3) * Poly(-1.0, lambda[k], d);
        }
      }
    }
  }

  for (size_t j = 0; j < NY; ++j) {
    for (size_t i = 0; i < NX + 1; ++i) {
      for (size_t k = 0; k < NumGLP; ++k) {
        double rhoL = UL(j, i, k, 0);
        double uL = UL(j, i, k, 1) / rhoL;
        double vL = UL(j, i, k, 2) / rhoL;
        double EL = UL(j, i, k, 3);
        double pL = (gasGamma - 1.0) * (EL - 0.5 * rhoL * (uL * uL + vL * vL));
        FL(j, i, k, 0) = UL(j, i, k, 1);
        FL(j, i, k, 1) = rhoL * uL * uL + pL;
        FL(j, i, k, 2) = rhoL * uL * vL;
        FL(j, i, k, 3) = uL * (EL + pL);
        double rhoR = UR(j, i, k, 0);
        double uR = UR(j, i, k, 1) / rhoR;
        double vR = UR(j, i, k, 2) / rhoR;
        double ER = UR(j, i, k, 3);
        double pR = (gasGamma - 1.0) * (ER - 0.5 * rhoR * (uR * uR + vR * vR));
        FR(j, i, k, 0) = UR(j, i, k, 1);
        FR(j, i, k, 1) = rhoR * uR * uR + pR;
        FR(j, i, k, 2) = rhoR * uR * vR;
        FR(j, i, k, 3) = uR * (ER + pR);
      }
    }
  }

  for (size_t j = 0; j < NY + 1; ++j) {
    for (size_t i = 0; i < NX; ++i) {
      for (size_t k = 0; k < NumGLP; ++k) {
        for (size_t d = 0; d < dimPk; ++d) {
          UU(j, i, k, 0) += uh(j + 1, i + 1, d, 0) * Poly(lambda[k], -1.0, d);
          UU(j, i, k, 1) += uh(j + 1, i + 1, d, 1) * Poly(lambda[k], -1.0, d);
          UU(j, i, k, 2) += uh(j + 1, i + 1, d, 2) * Poly(lambda[k], -1.0, d);
          UU(j, i, k, 3) += uh(j + 1, i + 1, d, 3) * Poly(lambda[k], -1.0, d);
          UD(j, i, k, 0) += uh(j, i + 1, d, 0) * Poly(lambda[k], 1.0, d);
          UD(j, i, k, 1) += uh(j, i + 1, d, 1) * Poly(lambda[k], 1.0, d);
          UD(j, i, k, 2) += uh(j, i + 1, d, 2) * Poly(lambda[k], 1.0, d);
          UD(j, i, k, 3) += uh(j, i + 1, d, 3) * Poly(lambda[k], 1.0, d);
        }
      }
    }
  }

  for (size_t j = 0; j < NY + 1; ++j) {
    for (size_t i = 0; i < NX; ++i) {
      for (size_t k = 0; k < NumGLP; ++k) {
        double rhoU = UU(j, i, k, 0);
        double uU = UU(j, i, k, 1) / rhoU;
        double vU = UU(j, i, k, 2) / rhoU;
        double EU = UU(j, i, k, 3);
        double pU = (gasGamma - 1.0) * (EU - 0.5 * rhoU * (uU * uU + vU * vU));
        FU(j, i, k, 0) = UU(j, i, k, 2);
        FU(j, i, k, 1) = rhoU * uU * vU;
        FU(j, i, k, 2) = rhoU * vU * vU + pU;
        FU(j, i, k, 3) = vU * (EU + pU);
        double rhoD = UD(j, i, k, 0);
        double uD = UD(j, i, k, 1) / rhoD;
        double vD = UD(j, i, k, 2) / rhoD;
        double ED = UD(j, i, k, 3);
        double pD = (gasGamma - 1.0) * (ED - 0.5 * rhoD * (uD * uD + vD * vD));
        FD(j, i, k, 0) = UD(j, i, k, 2);
        FD(j, i, k, 1) = rhoD * uD * vD;
        FD(j, i, k, 2) = rhoD * vD * vD + pD;
        FD(j, i, k, 3) = vD * (ED + pD);
      }
    }
  }

  xt::xarray<double> Fxhat({NY, NX + 1, NumGLP, NumEq}, 0.0);
  xt::xarray<double> Fyhat({NY + 1, NX, NumGLP, NumEq}, 0.0);

  for (size_t j = 0; j < NY; ++j) {
    for (size_t i = 0; i < NX + 1; ++i) {
      for (size_t k = 0; k < NumGLP; ++k) {
        double SRmax = 0.0;
        double SRmin = 0.0;
        double SLmax = 0.0;
        double SLmin = 0.0;
        eigenvalueMm(SRmax, SRmin, UR(j, i, k, 0), UR(j, i, k, 1),
                     UR(j, i, k, 2), UR(j, i, k, 3), 1.0, 0.0);
        eigenvalueMm(SLmax, SLmin, UL(j, i, k, 0), UL(j, i, k, 1),
                     UL(j, i, k, 2), UL(j, i, k, 3), 1.0, 0.0);
        double SR = std::max(SRmax, SLmax);
        double SL = std::min(SRmin, SLmin);
        for (size_t n = 0; n < NumEq; ++n) {
          double FR1 = FR(j, i, k, n);
          double FL1 = FL(j, i, k, n);
          double UR1 = UR(j, i, k, n);
          double UL1 = UL(j, i, k, n);
          Fxhat(j, i, k, n) = HLL_FLUX(SL, SR, FL1, FR1, UL1, UR1);
        }
      }
    }
  }

  for (size_t j = 0; j < NY + 1; ++j) {
    for (size_t i = 0; i < NX; ++i) {
      for (size_t k = 0; k < NumGLP; ++k) {
        double SRmax = 0.0;
        double SRmin = 0.0;
        double SLmax = 0.0;
        double SLmin = 0.0;
        eigenvalueMm(SRmax, SRmin, UU(j, i, k, 0), UU(j, i, k, 1),
                     UU(j, i, k, 2), UU(j, i, k, 3), 0.0, 1.0);
        eigenvalueMm(SLmax, SLmin, UD(j, i, k, 0), UD(j, i, k, 1),
                     UD(j, i, k, 2), UD(j, i, k, 3), 0.0, 1.0);
        double SR = std::max(SRmax, SLmax);
        double SL = std::min(SRmin, SLmin);
        for (size_t n = 0; n < NumEq; ++n) {
          double FR1 = FU(j, i, k, n);
          double FL1 = FD(j, i, k, n);
          double UR1 = UU(j, i, k, n);
          double UL1 = UD(j, i, k, n);
          Fyhat(j, i, k, n) = HLL_FLUX(SL, SR, FL1, FR1, UL1, UR1);
        }
      }
    }
  }

  for (size_t j = 1; j < NY + 1; ++j) {
    for (size_t i = 1; i < NX + 1; ++i) {
      for (size_t k = 0; k < dimPk; ++k) {
        for (size_t n = 0; n < NumEq; ++n) {
          for (size_t j1 = 0; j1 < NumGLP; ++j1) {
            du(j, i, k, n) =
                du(j, i, k, n) -
                0.5 / hx * weight[j1] *
                    (Fxhat(j - 1, i, j1, n) * Poly(1.0, lambda[j1], k) -
                     Fxhat(j - 1, i - 1, j1, n) * Poly(-1.0, lambda[j1], k));
          }
        }
      }
    }
  }

  for (size_t j = 1; j < NY + 1; ++j) {
    for (size_t i = 1; i < NX + 1; ++i) {
      for (size_t k = 0; k < dimPk; ++k) {
        for (size_t n = 0; n < NumEq; ++n) {
          for (size_t i1 = 0; i1 < NumGLP; ++i1) {
            du(j, i, k, n) =
                du(j, i, k, n) -
                0.5 / hy * weight[i1] *
                    (Fyhat(j, i - 1, i1, n) * Poly(lambda[i1], 1.0, k) -
                     Fyhat(j - 1, i - 1, i1, n) * Poly(lambda[i1], -1.0, k));
          }
        }
      }
    }
  }

  for (size_t j = 0; j < NY + 2; ++j) {
    for (size_t i = 0; i < NX + 2; ++i) {
      for (size_t k = 0; k < dimPk; ++k) {
        for (size_t n = 0; n < NumEq; ++n) {
          du(j, i, k, n) = du(j, i, k, n) / mass[k];
        }
      }
    }
  }
  return du;
}