#pragma once
#include <vector>
inline auto matMul(const std::vector<std::vector<double>> &mat,
                   const std::vector<double> &vec, const int &n) {
  std::vector<double> result(n, 0.0);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      result[i] += mat[i][j] * vec[j];
    }
  }
  return result;
}