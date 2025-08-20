#pragma once
#include "global.hpp"
#include <fstream>
#include <xtensor.hpp>
#include <xtensor/core/xtensor_forward.hpp>
inline void output(const xt::xarray<double> &uh) {
  std::ofstream fout1("u1.txt");
  fout1 << std::scientific << std::setprecision(15);
  for (size_t k = 0; k < dimPk; ++k) {
    for (size_t j = 1; j < NY + 1; ++j) {
      for (size_t i = 1; i < NX + 1; ++i) {
        fout1 << uh(j, i, k, 0) << "\n";
      }
    }
  }
  fout1.close();

  std::ofstream fout2("u2.txt");
  fout2 << std::scientific << std::setprecision(15);
  for (size_t k = 0; k < dimPk; ++k) {
    for (size_t j = 1; j < NY + 1; ++j) {
      for (size_t i = 1; i < NX + 1; ++i) {
        fout2 << uh(j, i, k, 1) << "\n";
      }
    }
  }
  fout2.close();

  std::ofstream fout3("u3.txt");
  fout3 << std::scientific << std::setprecision(15);
  for (size_t k = 0; k < dimPk; ++k) {
    for (size_t j = 1; j < NY + 1; ++j) {
      for (size_t i = 1; i < NX + 1; ++i) {
        fout3 << uh(j, i, k, 2) << "\n";
      }
    }
  }
  fout3.close();

  std::ofstream fout4("u4.txt");
  fout4 << std::scientific << std::setprecision(15);
  for (size_t k = 0; k < dimPk; ++k) {
    for (size_t j = 1; j < NY + 1; ++j) {
      for (size_t i = 1; i < NX + 1; ++i) {
        fout4 << uh(j, i, k, 3) << "\n";
      }
    }
  }
  fout4.close();
}