#include "Lh.hpp"
#include "TVB.hpp"
#include "basis.hpp"
#include "boundary.hpp"
#include "calculateDt.hpp"
#include "global.hpp"
#include "init.hpp"
#include "output.hpp"
#include "pmin.hpp"
#include "umax.hpp"
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <ios>
#include <iostream>
#include <vector>
#include <xtensor.hpp>
#include <xtensor/core/xtensor_forward.hpp>

int main() {
  std::vector<double> Xc(NX, 0.0);
  std::vector<double> Yc(NY, 0.0);

  for (size_t i = 0; i < NX; ++i) {
    Xc[i] = xa + (i + 0.5) * hx;
  }
  for (size_t i = 0; i < NY; ++i) {
    Yc[i] = ya + (i + 0.5) * hy;
  }

  xt::xarray<double> uh({NY + 2, NX + 2, dimPk, NumEq}, 0.0);
  xt::xarray<double> du({NY + 2, NX + 2, dimPk, NumEq}, 0.0);
  xt::xarray<double> u1({NY + 2, NX + 2, dimPk, NumEq}, 0.0);
  xt::xarray<double> u2({NY + 2, NX + 2, dimPk, NumEq}, 0.0);
  xt::xarray<double> uh0({NY + 2, NX + 2, dimPk, NumEq}, 0.0);

  for (size_t j = 1; j < NY + 1; ++j) {
    for (size_t i = 1; i < NX + 1; ++i) {
      for (size_t k = 0; k < dimPk; ++k) {
        for (size_t j1 = 0; j1 < NumGLP; ++j1) {
          for (size_t i1 = 0; i1 < NumGLP; ++i1) {
            double x = Xc[i - 1] + hx1 * lambda[i1];
            double y = Yc[j - 1] + hy1 * lambda[j1];
            uh(j, i, k, 0) += 0.25 * weight[i1] * weight[j1] *
                              Poly(lambda[i1], lambda[j1], k) * U1(x, y) /
                              mass[k];
            uh(j, i, k, 1) += 0.25 * weight[i1] * weight[j1] *
                              Poly(lambda[i1], lambda[j1], k) * U2(x, y) /
                              mass[k];
            uh(j, i, k, 2) += 0.25 * weight[i1] * weight[j1] *
                              Poly(lambda[i1], lambda[j1], k) * U3(x, y) /
                              mass[k];
            uh(j, i, k, 3) += 0.25 * weight[i1] * weight[j1] *
                              Poly(lambda[i1], lambda[j1], k) * U4(x, y) /
                              mass[k];
          }
        }
      }
    }
  }

  set_bc(uh);

  TVB_limiter(uh);

  // output(uh);

  auto umax = calculate_umax(uh);
  auto pmin = calculate_pmin(uh);

  double t = 0.0;
  std::cout << std::fixed << std::setprecision(16) << "t = " << t
            << ", umax = " << umax << ", pmin = " << pmin << std::endl;
  while (t < tend) {
    auto dt = calculate_dt(uh);
    if (t + dt >= tend) {
      dt = tend - t;
      t = tend;
    } else {
      t += dt;
    }

    set_bc(uh);

    du = dgOperator(uh);
    u1 = uh + dt * du;

    uh0 = uh;
    uh = u1;

    TVB_limiter(uh);

    set_bc(uh);

    du = dgOperator(uh);

    uh = u1 + dt * du;

    u2 = 3.0 / 4.0 * uh0 + 1.0 / 4.0 * uh;
    uh = u2;
    TVB_limiter(uh);

    set_bc(uh);
    du = dgOperator(uh);
    uh = u2 + dt * du;
    uh = 1.0 / 3.0 * uh0 + 2.0 / 3.0 * uh;
    TVB_limiter(uh);

    umax = calculate_umax(uh);
    pmin = calculate_pmin(uh);

    std::cout << std::fixed << std::setprecision(16) << "t = " << t
              << ", umax = " << umax << ", pmin = " << pmin << std::endl;
  }
  output(uh);
}
