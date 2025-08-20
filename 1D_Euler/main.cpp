#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <tuple>
#include <utility>
#include <vector>

#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

#define sgn(x) ((x) > 0 ? 1 : -1)

constexpr double xleft = 0.0;
constexpr double xright = 1.0;
constexpr int NX = 601;
constexpr int k = 2;
constexpr int NumGLP = 5;
constexpr int dimPK = k + 1;
constexpr int NumEq = 3;
constexpr double CFL = 0.15;
constexpr double tend = 0.012;
constexpr double speed = 1.0;
constexpr double deltaX = (xright - xleft) / NX;
constexpr double halfDeltax = deltaX * 0.5;

constexpr double gasgamma = 1.4;

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

double initPressure(const double &x) {
  if (x < 0.5) {
    return 1000.0;
  } else {
    return 0.01;
  }
}

double initRho(const double &x) {
  if (x < 0.5) {
    return 1.0;
  } else {
    return 1.0;
  }
}

double initVel(const double &x) {
  if (x < 0.5) {
    return 10.0;
    ;
  } else {
    return 0.0;
  }
}

double flux1(const double &u1, const double &u2, const double &u3) {
  return u2;
}

double flux2(const double &u1, const double &u2, const double &u3) {
  double p = (gasgamma - 1.0) * (u3 - 0.5 * u2 * u2 / u1);
  return u2 * u2 / u1 + p;
}

double flux3(const double &u1, const double &u2, const double &u3) {
  double p = (gasgamma - 1.0) * (u3 - 0.5 * u2 * u2 / u1);
  return u2 / u1 * (u3 + p);
}

std::pair<double, double> waveSpeed(const double &u1, const double &u2,
                                    const double &u3) {
  double p = (gasgamma - 1.0) * (u3 - 0.5 * u2 * u2 / u1);
  double c = std::sqrt(gasgamma * p / u1);
  double SR = u2 / u1 + c;
  double SL = u2 / u1 - c;
  return std::make_pair(SR, SL);
}

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

void init(std::vector<double> &vec, const double &xleft, const double &xright,
          const int &NX, const int &numGLP, const int &numEq) {
  double deltaX = (xright - xleft) / NX;
  double halfDeltaX = 0.5 * deltaX;

  std::vector<double> xCenter(NX);

  for (int i = 0; i < xCenter.size(); ++i) {
    xCenter[i] = (2 * i + 1) * halfDeltaX + xleft;
  }

  for (int i = 0; i < NX; ++i) {
    for (int j = 0; j < numGLP; ++j) {
      vec[i * numGLP * numEq + j * numEq] =
          initRho(xCenter[i] + halfDeltaX * lambda[j]);
      vec[i * numGLP * numEq + j * numEq + 1] =
          initRho(xCenter[i] + halfDeltaX * lambda[j]) *
          initVel(xCenter[i] + halfDeltaX * lambda[j]);
      vec[i * numGLP * numEq + j * numEq + 2] =
          initPressure(xCenter[i] + halfDeltaX * lambda[j]) / (gasgamma - 1.0) +
          0.5 * initRho(xCenter[i] + halfDeltaX * lambda[j]) *
              initVel(xCenter[i] + halfDeltaX * lambda[j]) *
              initVel(xCenter[i] + halfDeltaX * lambda[j]);
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
             const int &nx, const int &dimpk, const int &numEq,
             const std::string &name) {
  std::vector<double> rho(nx, 0.0);
  std::vector<double> vel(nx, 0.0);
  std::vector<double> p(nx, 0.0);
  for (int i = 0; i < nx; ++i) {
    rho[i] = uh[i * dimpk * numEq];
    vel[i] = uh[i * dimpk * numEq + 1] / rho[i];
    p[i] =
        (uh[i * dimpk * numEq + 2] -
         0.5 * uh[i * dimpk * numEq + 1] * uh[i * dimpk * numEq + 1] / rho[i]) *
        (gasgamma - 1.0);
  }
  plot2D(x, rho, "rho_" + name);
  plot2D(x, vel, "u_" + name);
  plot2D(x, p, "p_" + name);
}

void plotProject(const std::vector<double> &x, const std::vector<double> &uCell,
                 const int &nx, const int &numGLP, const int &numEq,
                 const std::string &name) {
  std::vector<double> rho(nx * numGLP, 0.0);
  std::vector<double> vel(nx * numGLP, 0.0);
  std::vector<double> p(nx * numGLP, 0.0);
  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < numGLP; ++j) {
      rho[i * numGLP + j] = uCell[i * numGLP * numEq + j * numEq];
      vel[i * numGLP + j] = uCell[i * numGLP * numEq + j * numEq + 1] / rho[i];
      p[i * numGLP + j] =
          (uCell[i * numGLP * numEq + j * numEq + 2] -
           0.5 * uCell[i * numGLP * numEq + j * numEq + 1] *
               uCell[i * numGLP * numEq + j * numEq + 1] / rho[i]) *
          (gasgamma - 1.0);
    }
  }
  plot2D(x, rho, "rho_" + name);
  plot2D(x, vel, "u_" + name);
  plot2D(x, p, "p_" + name);
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
                const int &nx, const int &numGLP, const int &dimpk,
                const int &numEq) {
  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < numGLP; ++j) {
      for (int n = 0; n < numEq; ++n) {
        output[i * numGLP * numEq + j * numEq + n] = 0.0;
        for (int k = 0; k < dimpk; ++k) {
          output[i * numGLP * numEq + j * numEq + n] +=
              input[i * dimpk * numEq + k * numEq + n] * polyPK(lambda[j], k);
        }
      }
    }
  }
}

auto fluxHLLC(const std::tuple<double, double, double> &UL,
              const std::tuple<double, double, double> &UR,
              const std::tuple<double, double, double> &FL,
              const std::tuple<double, double, double> &FR, const double &SL,
              const double &SR) {
  if (SL >= 0) {
    return FL;
  } else if (SR <= 0) {
    return FR;
  }
  auto u1Left = std::get<0>(UL);
  auto u2Left = std::get<1>(UL);
  auto u3Left = std::get<2>(UL);
  auto u1Right = std::get<0>(UR);
  auto u2Right = std::get<1>(UR);
  auto u3Right = std::get<2>(UR);
  auto flux1Left = std::get<0>(FL);
  auto flux2Left = std::get<1>(FL);
  auto flux3Left = std::get<2>(FL);
  auto flux1Right = std::get<0>(FR);
  auto flux2Right = std::get<1>(FR);
  auto flux3Right = std::get<2>(FR);
  auto velLeft = u2Left / u1Left;
  auto velRight = u2Right / u1Right;
  auto pLeft = (gasgamma - 1.0) * (u3Left - 0.5 * u1Left * velLeft * velLeft);
  auto pRight =
      (gasgamma - 1.0) * (u3Right - 0.5 * u1Right * velRight * velRight);
  auto Sstar = (pRight - pLeft + u1Left * velLeft * (SL - velLeft) -
                u1Right * velRight * (SR - velRight)) /
               (u1Left * (SL - velLeft) - u1Right * (SR - velRight));
  auto Fhat1 = 0.0;
  auto Fhat2 = 0.0;
  auto Fhat3 = 0.0;
  if (SL < 0 && Sstar >= 0) {
    Fhat1 = (Sstar * (SL * u1Left - flux1Left) +
             SL * (pLeft + u1Left * (SL - velLeft) * (Sstar - velLeft)) * 0.0) /
            (SL - Sstar);
    Fhat2 = (Sstar * (SL * u2Left - flux2Left) +
             SL * (pLeft + u1Left * (SL - velLeft) * (Sstar - velLeft)) * 1.0) /
            (SL - Sstar);
    Fhat3 =
        (Sstar * (SL * u3Left - flux3Left) +
         SL * (pLeft + u1Left * (SL - velLeft) * (Sstar - velLeft)) * Sstar) /
        (SL - Sstar);
  } else {
    Fhat1 =
        (Sstar * (SR * u1Right - flux1Right) +
         SR * (pRight + u1Right * (SR - velRight) * (Sstar - velRight)) * 0.0) /
        (SR - Sstar);
    Fhat2 =
        (Sstar * (SR * u2Right - flux2Right) +
         SR * (pRight + u1Right * (SR - velRight) * (Sstar - velRight)) * 1.0) /
        (SR - Sstar);
    Fhat3 = (Sstar * (SR * u3Right - flux3Right) +
             SR * (pRight + u1Right * (SR - velRight) * (Sstar - velRight)) *
                 Sstar) /
            (SR - Sstar);
  }
  return std::make_tuple(Fhat1, Fhat2, Fhat3);
}

std::vector<double> dgOperator(const std::vector<double> &uh, const int &nx,
                               const int &numGLP, const int &dimpk,
                               const int &numEq) {
  std::vector<double> uhG(nx * numGLP * numEq, 0.0);
  std::vector<double> du(nx * dimpk * numEq, 0.0);
  projection(uh, uhG, nx, numGLP, dimpk, numEq);
  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < dimpk; ++j) {
      for (int k = 0; k < numGLP; ++k) {
        du[i * dimpk * numEq + j * numEq] +=
            0.5 * weight[k] *
            flux1(uhG[i * numGLP * numEq + k * numEq],
                  uhG[i * numGLP * numEq + k * numEq + 1],
                  uhG[i * numGLP * numEq + k * numEq + 2]) *
            polyxPK(lambda[k], halfDeltax, j);
        du[i * dimpk * numEq + j * numEq + 1] +=
            0.5 * weight[k] *
            flux2(uhG[i * numGLP * numEq + k * numEq],
                  uhG[i * numGLP * numEq + k * numEq + 1],
                  uhG[i * numGLP * numEq + k * numEq + 2]) *
            polyxPK(lambda[k], halfDeltax, j);
        du[i * dimpk * numEq + j * numEq + 2] +=
            0.5 * weight[k] *
            flux3(uhG[i * numGLP * numEq + k * numEq],
                  uhG[i * numGLP * numEq + k * numEq + 1],
                  uhG[i * numGLP * numEq + k * numEq + 2]) *
            polyxPK(lambda[k], halfDeltax, j);
      }
    }
  }

  std::vector<double> uLeft((nx + 1) * numEq, 0.0);
  std::vector<double> uRight((nx + 1) * numEq, 0.0);
  std::vector<double> fhat((nx + 1) * numEq, 0.0);

  for (int i = 0; i < nx + 1; ++i) {
    for (int j = 0; j < dimpk; ++j) {
      for (int n = 0; n < numEq; ++n) {
        if (i == 0) {
          uLeft[i * numEq + n] +=
              uh[i * dimpk * numEq + j * numEq + n] * polyPK(1.0, j);
          uRight[i * numEq + n] +=
              uh[i * dimpk * numEq + j * numEq + n] * polyPK(-1.0, j);
        } else if (i == nx) {
          uLeft[i * numEq + n] +=
              uh[(i - 1) * dimpk * numEq + j * numEq + n] * polyPK(1.0, j);
          uRight[i * numEq + n] +=
              uh[(i - 1) * dimpk * numEq + j * numEq + n] * polyPK(-1.0, j);
        } else {
          uLeft[i * numEq + n] +=
              uh[(i - 1) * dimpk * numEq + j * numEq + n] * polyPK(1.0, j);
          uRight[i * numEq + n] +=
              uh[i * dimpk * numEq + j * numEq + n] * polyPK(-1.0, j);
        }
      }
    }
  }
  for (int i = 0; i < nx + 1; ++i) {
    auto UL = std::make_tuple(uLeft[i * numEq], uLeft[i * numEq + 1],
                              uLeft[i * numEq + 2]);
    auto UR = std::make_tuple(uRight[i * numEq], uRight[i * numEq + 1],
                              uRight[i * numEq + 2]);
    auto FL = std::make_tuple(
        flux1(std::get<0>(UL), std::get<1>(UL), std::get<2>(UL)),
        flux2(std::get<0>(UL), std::get<1>(UL), std::get<2>(UL)),
        flux3(std::get<0>(UL), std::get<1>(UL), std::get<2>(UL)));
    auto FR = std::make_tuple(
        flux1(std::get<0>(UR), std::get<1>(UR), std::get<2>(UR)),
        flux2(std::get<0>(UR), std::get<1>(UR), std::get<2>(UR)),
        flux3(std::get<0>(UR), std::get<1>(UR), std::get<2>(UR)));
    auto [SLmax, SLmin] =
        waveSpeed(std::get<0>(UL), std::get<1>(UL), std::get<2>(UL));
    auto [SRmax, SRmin] =
        waveSpeed(std::get<0>(UR), std::get<1>(UR), std::get<2>(UR));
    auto SR = std::max(SRmax, SLmax);
    auto SL = std::min(SLmin, SRmin);

    auto Fhat = fluxHLLC(UL, UR, FL, FR, SL, SR);
    fhat[i * numEq] = std::get<0>(Fhat);
    fhat[i * numEq + 1] = std::get<1>(Fhat);
    fhat[i * numEq + 2] = std::get<2>(Fhat);
  }

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < dimpk; ++j) {
      for (int n = 0; n < numEq; ++n) {
        du[i * dimpk * numEq + j * numEq + n] =
            (du[i * dimpk * numEq + j * numEq + n] -
             (polyPK(1.0, j) * fhat[(i + 1) * numEq + n] -
              polyPK(-1.0, j) * fhat[i * numEq + n]) /
                 deltaX) /
            mass[j];
      }
    }
  }
  return du;
}

auto minmod(const std::vector<double> &var1, const std::vector<double> &var2,
            const std::vector<double> &var3, const int &numEq) {
  std::vector<double> result(numEq, 0.0);
  for (int n = 0; n < numEq; ++n) {
    if (std::abs(var1[n]) < (deltaX * deltaX)) {
      result[n] = var1[n];
    } else {
      if ((sgn(var1[n]) == sgn(var2[n])) && (sgn(var1[n]) == sgn(var3[n]))) {
        auto a1 = std::abs(var1[n]);
        auto a2 = std::abs(var2[n]);
        auto a3 = std::abs(var3[n]);
        auto min = std::min(a1, a2);
        result[n] = sgn(var1[n]) * std::min(min, a3);
      } else {
        result[n] = 0.0;
      }
    }
  }
  return result;
}

bool invertMatrix(const std::vector<std::vector<double>> &A,
                  std::vector<std::vector<double>> &inverse) {
  int n = A.size();
  inverse = std::vector<std::vector<double>>(n, std::vector<double>(n, 0.0));

  // 构造增广矩阵 [A | I]
  std::vector<std::vector<double>> aug(n, std::vector<double>(2 * n));
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      aug[i][j] = A[i][j];
    }
    aug[i][i + n] = 1.0; // 单位矩阵部分
  }

  // 高斯-约旦消元法
  for (int i = 0; i < n; ++i) {
    // 寻找主元素
    double pivot = aug[i][i];
    if (abs(pivot) < 1e-10) {
      // 尝试行交换
      bool found = false;
      for (int k = i + 1; k < n; ++k) {
        if (abs(aug[k][i]) > 1e-10) {
          swap(aug[i], aug[k]);
          pivot = aug[i][i];
          found = true;
          break;
        }
      }
      if (!found)
        return false; // 矩阵不可逆
    }

    // 归一化主元素所在行
    for (int j = 0; j < 2 * n; ++j) {
      aug[i][j] /= pivot;
    }

    // 消去其他行
    for (int k = 0; k < n; ++k) {
      if (k != i) {
        double factor = aug[k][i];
        for (int j = 0; j < 2 * n; ++j) {
          aug[k][j] -= factor * aug[i][j];
        }
      }
    }
  }

  // 提取逆矩阵
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
      inverse[i][j] = aug[i][j + n];

  return true;
}

auto matMul(const std::vector<std::vector<double>> &mat,
            const std::vector<double> &vec, const int &n) {
  std::vector<double> result(n, 0.0);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      result[i] += mat[i][j] * vec[j];
    }
  }
  return result;
}

void TVDlimiter(std::vector<double> &uh, const int &nx, const int &numGLP,
                const int &dimpk, const int &numEq) {
  for (int i = 0; i < nx; ++i) {
    std::vector<double> deltaUR(numEq, 0.0);
    std::vector<double> deltaUL(numEq, 0.0);
    std::vector<double> deltaURM(numEq, 0.0);
    std::vector<double> deltaULM(numEq, 0.0);
    for (int n = 0; n < numEq; ++n) {
      deltaUR[n] = uh[i * dimpk * numEq + 1 * numEq + n] +
                   2.0 / 3.0 * uh[i * dimpk * numEq + 2 * numEq + n];
      deltaUL[n] = uh[i * dimpk * numEq + 1 * numEq + n] -
                   2.0 / 3.0 * uh[i * dimpk * numEq + 2 * numEq + n];
      if (i == 0) {
        deltaURM[n] =
            uh[(i + 1) * dimpk * numEq + n] - uh[i * dimpk * numEq + n];
        deltaULM[n] = uh[i * dimpk * numEq + n] - uh[i * dimpk * numEq + n];
      } else if (i == (nx - 1)) {
        deltaURM[n] = uh[i * dimpk * numEq + n] - uh[i * dimpk * numEq + n];
        deltaULM[n] =
            uh[i * dimpk * numEq + n] - uh[(i - 1) * dimpk * numEq + n];
      } else {
        deltaURM[n] =
            uh[(i + 1) * dimpk * numEq + n] - uh[i * dimpk * numEq + n];
        deltaULM[n] =
            uh[i * dimpk * numEq + n] - uh[(i - 1) * dimpk * numEq + n];
      }
    }
    double v = uh[i * dimpk * numEq + 1] / uh[i * dimpk * numEq];
    double p = (gasgamma - 1.0) * (uh[i * dimpk * numEq + 2] -
                                   0.5 * uh[i * dimpk * numEq] * v * v);
    double c = std::sqrt(gasgamma * p / uh[i * dimpk * numEq]);
    double H = (uh[i * dimpk * numEq + 2] + p) / uh[i * dimpk * numEq];
    std::vector<std::vector<double>> R{
        std::vector<double>{1.0, 1.0, 1.0},
        std::vector<double>{v - c, v, v + c},
        std::vector<double>{H - v * c, 0.5 * v * v, H + v * c}};
    std::vector<std::vector<double>> L(3, std::vector<double>{0.0, 0.0, 0.0});
    invertMatrix(R, L);

    auto LdeltaUR = matMul(L, deltaUR, numEq);
    auto LdeltaUL = matMul(L, deltaUL, numEq);
    auto LdeltaURM = matMul(L, deltaURM, numEq);
    auto LdeltaULM = matMul(L, deltaULM, numEq);
    auto deltaURM1 = minmod(LdeltaUR, LdeltaURM, LdeltaULM, numEq);
    auto deltaULM1 = minmod(LdeltaUL, LdeltaURM, LdeltaULM, numEq);
    deltaURM1 = matMul(R, deltaURM1, numEq);
    deltaULM1 = matMul(R, deltaULM1, numEq);

    for (int n = 0; n < numEq; ++n) {
      uh[i * dimpk * numEq + 1 * dimpk + n] =
          (deltaURM1[n] + deltaULM1[n]) * 0.5;
      uh[i * dimpk * numEq + 2 * dimpk + n] =
          0.75 * (deltaURM1[n] - deltaULM1[n]);
    }
  }
}

double pressure(const double &u1, const double &u2, const double &u3) {
  return (gasgamma - 1.0) * (u3 - 0.5 * u2 * u2 / u1);
}

void ppLimiter(std::vector<double> uh, const int &nx, const int &numGLP,
               const int &dimpk, const int &numEq) {
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
  std::vector<double> uInit(NX * NumGLP * NumEq);
  for (int i = 0; i < x.size(); ++i) {
    x[i] = (xright - xleft) / (NX * NumGLP) * i + xleft;
  }
  for (int i = 0; i < xplot.size(); ++i) {
    xplot[i] = (xright - xleft) / NX * i + xleft;
  }
  init(uInit, xleft, xright, NX, NumGLP, NumEq);

  std::vector<double> uh(NX * dimPK * NumEq, 0.0);
  for (int i = 0; i < NX; ++i) {
    for (int j = 0; j < dimPK; ++j) {
      for (int k = 0; k < NumGLP; ++k) {
        for (int n = 0; n < NumEq; ++n) {
          uh[i * dimPK * NumEq + j * NumEq + n] +=
              1.0 / (2.0 * mass[j]) * weight[k] *
              uInit[i * NumGLP * NumEq + k * NumEq + n] * polyPK(lambda[k], j);
        }
      }
    }
  }

  std::vector<double> uCell(NX * NumGLP * NumEq, 0.0);
  projection(uh, uCell, NX, NumGLP, dimPK, NumEq);
  double alpha = 1.0;
  for (int i = 0; i < NX; ++i) {
    double alpha1 = waveSpeed(uh[i * dimPK * NumEq], uh[i * dimPK * NumEq + 1],
                              uh[i * dimPK * NumEq + 2])
                        .first;
    if (alpha1 > alpha) {
      alpha = alpha1;
    }
  }

  double dt = CFL * (deltaX) / alpha;
  int NSTEPS = ceil(tend / dt);
  std::cout << "dt: " << dt << "\n";
  std::cout << "NSTEPS: " << NSTEPS << "\n";
  for (int t = 0; t < NSTEPS; ++t) {
    if (t % 500 == 0) {
      std::cout << "t = " << t << "\n";
    }
    auto du = dgOperator(uh, NX, NumGLP, dimPK, NumEq);
    auto u1 = vecAdd2(uh, du, 1.0, dt);
    TVDlimiter(u1, NX, NumGLP, dimPK, NumEq);

    du = dgOperator(u1, NX, NumGLP, dimPK, NumEq);
    auto u2 = vecAdd3(uh, u1, du, 0.75, 0.25, 0.25 * dt);
    TVDlimiter(u2, NX, NumGLP, dimPK, NumEq);

    du = dgOperator(u2, NX, NumGLP, dimPK, NumEq);
    uh = vecAdd3(uh, u2, du, 1.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0 * dt);
    TVDlimiter(uh, NX, NumGLP, dimPK, NumEq);
  }

  plotAve(xplot, uh, NX, dimPK, NumEq, "result.png");
}