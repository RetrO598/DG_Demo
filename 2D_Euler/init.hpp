#pragma once
#include "global.hpp"
inline double rho1(const double &x, const double &y) {
  double rho1 = 0.0;
  if (x > 0.8 && y > 0.8) {
    rho1 = 1.5;
  } else if (x <= 0.8 && y > 0.8) {
    rho1 = 0.5323;
  } else if (x <= 0.8 && y <= 0.8) {
    rho1 = 0.138;
  } else {
    rho1 = 0.5323;
  }

  return rho1;
}

inline double ux1(const double &x, const double &y) {
  double ux1 = 0.0;
  if (x > 0.8 && y > 0.8) {
    ux1 = 0.0;
  } else if (x <= 0.8 && y > 0.8) {
    ux1 = 1.206;
  } else if (x <= 0.8 && y <= 0.8) {
    ux1 = 1.206;
  } else {
    ux1 = 0.0;
  }
  return ux1;
}

inline double uy1(const double &x, const double &y) {
  double uy1 = 0.0;
  if (x > 0.8 && y > 0.8) {
    uy1 = 0.0;
  } else if (x <= 0.8 && y > 0.8) {
    uy1 = 0.0;
  } else if (x <= 0.8 && y <= 0.8) {
    uy1 = 1.206;
  } else {
    uy1 = 1.206;
  }
  return uy1;
}

inline double pp1(const double &x, const double &y) {
  double pp1 = 0.0;
  if (x > 0.8 && y > 0.8) {
    pp1 = 1.5;
  } else if (x <= 0.8 && y > 0.8) {
    pp1 = 0.3;
  } else if (x <= 0.8 && y <= 0.8) {
    pp1 = 0.029;
  } else {
    pp1 = 0.3;
  }
  return pp1;
}

inline double U1(const double &x, const double &y) { return rho1(x, y); }

inline double U2(const double &x, const double &y) {
  return rho1(x, y) * ux1(x, y);
}

inline double U3(const double &x, const double &y) {
  return rho1(x, y) * uy1(x, y);
}

inline double U4(const double &x, const double &y) {
  return pp1(x, y) / (gasGamma - 1.0) +
         0.5 * rho1(x, y) * (ux1(x, y) * ux1(x, y) + uy1(x, y) * uy1(x, y));
}