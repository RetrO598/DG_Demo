#pragma once

inline double pressure(const double &rho, const double &rhou,
                       const double &rhov, const double &E,
                       const double &gamma) {
  return (gamma - 1.0) * (E - 0.5 * (rhou * rhou + rhov * rhov) / rho);
}