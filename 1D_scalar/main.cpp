#include "matplotlibcpp.h"
#include <algorithm>
#include <array>
#include <blaze/Math.h>
#include <cmath>
#include <cstdlib>
#include <vector>

namespace plt = matplotlibcpp;

#define sgn(x) ((x) > 0 ? 1 : -1)

constexpr double xleft = -1.0;
constexpr double xright = 1.0;
constexpr int NX = 160;
constexpr int k = 2;
constexpr int NumGLP = 5;
constexpr int dimPK = k + 1;
constexpr double CFL = 0.2f;
constexpr double tend = 0.4;
constexpr double speed = 1.0;
constexpr double deltaX = (xright - xleft) / NX;
constexpr double halfDeltax = deltaX * 0.5;

constexpr double max = 1.0;
constexpr double min = 0.0;

constexpr std::array<double, 5> lambda = {
    -0.9061798459386639927976269, -0.5384693101056830910363144, 0.0,
    0.5384693101056830910363144, 0.9061798459386639927976269};

constexpr std::array<double, 5> weight = {
    0.2369268850561890875142640, 0.4786286704993664680412915,
    0.5688888888888888888888889, 0.4786286704993664680412915,
    0.2369268850561890875142640};

constexpr std::array<double, 5> lambdaL = {-1.0, -0.6546536707079771437983, 0.0,
                                           0.654653670707977143798, 1.0};

constexpr std::array<double, 3> mass = {1.0, 1.0 / 3.0, 4.0 / 45.0};

std::vector<double> vecAdd2(const std::vector<double> &vec1,
                            const std::vector<double> &vec2,
                            const double &alpha, const double &beta) {
  std::vector<double> result(vec1.size(), 0.0);
  for (int i = 0; i < vec1.size(); ++i) {
    result[i] = vec1[i] * alpha + vec2[i] * beta;
  }
  return result;
}

std::vector<double> vecAdd3(const std::vector<double> &vec1,
                            const std::vector<double> &vec2,
                            const std::vector<double> &vec3,
                            const double &alpha, const double &beta,
                            const double &gamma) {
  std::vector<double> result(vec1.size(), 0.0);
  for (int i = 0; i < vec1.size(); ++i) {
    result[i] = vec1[i] * alpha + vec2[i] * beta + vec3[i] * gamma;
  }
  return result;
}

template <typename Func>
void init(std::vector<double> &vec, const double &xleft, const double &xright,
          const int &NX, const int &numGLP, Func func) {
  double deltaX = (xright - xleft) / NX;
  double halfDeltaX = 0.5 * deltaX;

  std::vector<double> xCenter(NX);

  for (int i = 0; i < xCenter.size(); ++i) {
    xCenter[i] = (2 * i + 1) * halfDeltaX + xleft;
  }

  for (int i = 0; i < NX; ++i) {
    for (int j = 0; j < numGLP; ++j) {
      vec[j + i * numGLP] = func(xCenter[i] + halfDeltaX * lambda[j]);
    }
  }
}

void plot2D(const std::vector<double> &x, const std::vector<double> &y,
            const std::string &name) {
  plt::figure_size(1200, 780);
  plt::plot(x, y);
  plt::save(name, 500);
}

void plotAve(const std::vector<double> &x, const std::vector<double> &uh,
             const int &nx, const int &dimpk, const std::string &name) {
  std::vector<double> y(nx, 0.0);
  for (int i = 0; i < nx; ++i) {
    y[i] = uh[i * dimpk];
  }
  plot2D(x, y, name);
}

double polyPK(const double &x, const int &order) {
  if (order == 0) {
    return 1.0;
  } else if (order == 1) {
    return x;
  } else {
    return x * x - 1.0 / 3.0;
  }
}

double polyxPK(const double &x, const double &halfDeltax, const int &order) {
  if (order == 0) {
    return 0.0;
  } else if (order == 1) {
    return 1.0 / halfDeltax;
  } else {
    return 2.0 * x / halfDeltax;
  }
}

void projection(const std::vector<double> &input, std::vector<double> &output,
                const int &nx, const int &numGLP, const int &dimpk) {
  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < numGLP; ++j) {
      output[i * numGLP + j] = 0.0;
      for (int k = 0; k < dimpk; ++k) {
        output[i * numGLP + j] += input[i * dimpk + k] * polyPK(lambda[j], k);
      }
    }
  }
}
template <typename Func>
std::vector<double> dgOperator(const std::vector<double> &uh, const int &nx,
                               const int &numGLP, const int &dimpk,
                               const double &alpha, Func func) {
  std::vector<double> uhG(nx * numGLP, 0.0);
  std::vector<double> du(nx * dimpk, 0.0);
  projection(uh, uhG, nx, numGLP, dimpk);
  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < dimpk; ++j) {
      for (int k = 0; k < numGLP; ++k) {
        du[i * dimpk + j] += 0.5 * weight[k] * func(uhG[i * numGLP + k]) *
                             polyxPK(lambda[k], halfDeltax, j);
      }
    }
  }

  std::vector<double> fluxLeft((nx + 1), 0.0);
  std::vector<double> fluxRight((nx + 1), 0.0);
  std::vector<double> fhat((nx + 1), 0.0);

  for (int i = 0; i < nx + 1; ++i) {
    for (int j = 0; j < dimpk; ++j) {
      if (i == 0) {
        fluxLeft[i] += uh[(nx - 1) * dimpk + j] * polyPK(1.0, j);
        fluxRight[i] += uh[i * dimpk + j] * polyPK(-1.0, j);
      } else if (i == nx) {
        fluxLeft[i] += uh[(i - 1) * dimpk + j] * polyPK(1.0, j);
        fluxRight[i] += uh[j] * polyPK(-1.0, j);
      } else {
        fluxLeft[i] += uh[(i - 1) * dimpk + j] * polyPK(1.0, j);
        fluxRight[i] += uh[i * dimpk + j] * polyPK(-1.0, j);
      }
    }
  }

  for (int i = 0; i < nx + 1; ++i) {
    fhat[i] = 0.5 * (func(fluxLeft[i]) + func(fluxRight[i]) -
                     alpha * (fluxRight[i] - fluxLeft[i]));
  }

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < dimpk; ++j) {
      du[i * dimpk + j] = (du[i * dimpk + j] - (polyPK(1.0, j) * fhat[i + 1] -
                                                polyPK(-1.0, j) * fhat[i]) /
                                                   deltaX) /
                          mass[j];
    }
  }
  return du;
}

double minmod(const double &var1, const double &var2, const double &var3) {
  if (std::abs(var1) < (deltaX * deltaX)) {
    return var1;
  } else {
    if ((sgn(var1) == sgn(var2)) && (sgn(var1) == sgn(var3))) {
      auto a1 = std::abs(var1);
      auto a2 = std::abs(var2);
      auto a3 = std::abs(var3);
      auto min = std::min(a1, a2);
      return sgn(var1) * std::min(min, a3);
    } else {
      return 0.0;
    }
  }
}

void TVDlimiter(std::vector<double> &uh, const int &nx, const int &numGLP,
                const int &dimpk) {
  for (int i = 0; i < nx; ++i) {
    double deltaUR = uh[i * dimpk + 1] + 2.0 / 3.0 * uh[i * dimpk + 2];
    double deltaUL = uh[i * dimpk + 1] - 2.0 / 3.0 * uh[i * dimpk + 2];
    double deltaURM = 0.0;
    double deltaULM = 0.0;
    if (i == 0) {
      deltaURM = uh[(i + 1) * dimpk] - uh[i * dimpk];
      deltaULM = uh[i * dimpk] - uh[(nx - 1) * dimpk];
    } else if (i == (nx - 1)) {
      deltaURM = uh[0] - uh[i * dimpk];
      deltaULM = uh[i * dimpk] - uh[(i - 1) * dimpk];
    } else {
      deltaURM = uh[(i + 1) * dimpk] - uh[i * dimpk];
      deltaULM = uh[i * dimpk] - uh[(i - 1) * dimpk];
    }
    double deltaURM1 = minmod(deltaUR, deltaURM, deltaULM);
    double deltaULM1 = minmod(deltaUL, deltaURM, deltaULM);

    uh[i * dimpk + 1] = (deltaURM1 + deltaULM1) * 0.5;
    uh[i * dimpk + 2] = 0.75 * (deltaURM1 - deltaULM1);
  }
}

void maxLimiter(std::vector<double> uh, const int &nx, const int &numGLP,
                const int &dimpk, const double &max, const double &min) {
  std::vector<double> uGLP(nx * numGLP, 0.0);
  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < numGLP; ++j) {
      for (int k = 0; k < dimpk; ++k) {
        uGLP[i * numGLP + j] += uh[i * dimpk + k] * polyPK(lambdaL[j], k);
      }
    }
  }

  for (int i = 0; i < nx; ++i) {
    double theta = 1.0;
    double thetaq = 0.0;
    for (int j = 0; j < numGLP; ++j) {
      thetaq =
          std::min(std::abs((max - uh[i * dimpk]) / (uGLP[i * numGLP + j]) -
                            uh[i * dimpk]),
                   std::abs((min - uh[i * dimpk]) /
                            (uGLP[i * numGLP + j] - uh[i * dimpk])));
      thetaq = std::min(thetaq, 1.0);
      if (thetaq < theta) {
        theta = thetaq;
      }
    }
    uh[i * dimpk + 1] *= theta;
    uh[i * dimpk + 2] *= theta;
  }
}

int main() {
  std::vector<double> x(NX * NumGLP);
  std::vector<double> xplot(NX);
  std::vector<double> uInit(NX * NumGLP);
  for (int i = 0; i < x.size(); ++i) {
    x[i] = (xright - xleft) / (NX * NumGLP) * i + xleft;
  }
  for (int i = 0; i < xplot.size(); ++i) {
    xplot[i] = (xright - xleft) / NX * i + xleft;
  }
  init(uInit, xleft, xright, NX, NumGLP, [](double x) {
    if (x >= -0.5 && x <= 0) {
      return 1.0;
    } else {
      return 0.0;
    }
  });

  std::vector<double> uh(NX * dimPK, 0.0);
  for (int i = 0; i < NX; ++i) {
    for (int j = 0; j < dimPK; ++j) {
      for (int k = 0; k < NumGLP; ++k) {
        uh[i * dimPK + j] += 1.0 / (2.0 * mass[j]) * weight[k] *
                             uInit[i * NumGLP + k] * polyPK(lambda[k], j);
      }
    }
  }

  std::vector<double> uCell(NX * NumGLP, 0.0);
  projection(uh, uCell, NX, NumGLP, dimPK);
  plot2D(x, uCell, "init.png");
  double dt = CFL * (deltaX) / speed;

  int NSTEPS = ceil(tend / dt);

  for (int t = 0; t < NSTEPS; ++t) {
    std::cout << "t = " << t << "\n";
    auto du = dgOperator(uh, NX, NumGLP, dimPK, 2.35, [](double x) {
      return (4.0 * x * x) / (4.0 * x * x + (1.0 - x) * (1.0 - x));
    });
    auto u1 = vecAdd2(uh, du, 1.0, dt);
    TVDlimiter(u1, NX, NumGLP, dimPK);
    // maxLimiter(u1, NX, NumGLP, dimPK, max, min);

    du = dgOperator(u1, NX, NumGLP, dimPK, 2.35, [](double x) {
      return (4.0 * x * x) / (4.0 * x * x + (1.0 - x) * (1.0 - x));
    });
    auto u2 = vecAdd3(uh, u1, du, 0.75, 0.25, 0.25 * dt);
    TVDlimiter(u2, NX, NumGLP, dimPK);
    // maxLimiter(u2, NX, NumGLP, dimPK, max, min);

    du = dgOperator(u2, NX, NumGLP, dimPK, 2.35, [](double x) {
      return (4.0 * x * x) / (4.0 * x * x + (1.0 - x) * (1.0 - x));
    });
    uh = vecAdd3(uh, u2, du, 1.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0 * dt);
    TVDlimiter(uh, NX, NumGLP, dimPK);
    // maxLimiter(uh, NX, NumGLP, dimPK, max, min);
  }
  projection(uh, uCell, NX, NumGLP, dimPK);
  plotAve(xplot, uh, NX, dimPK, "result.png");
}