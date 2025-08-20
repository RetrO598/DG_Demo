#pragma once

inline double HLL_FLUX(const double &SL, const double &SR, const double &FL1,
                       const double &FR1, const double &UL1,
                       const double &UR1) {
  double Fhat1 = 0.0;
  if (SR < 0) {
    Fhat1 = FR1;
  } else if (SL > 0) {
    Fhat1 = FL1;
  } else {
    Fhat1 = (SR * FL1 - SL * FR1 + SL * SR * (UR1 - UL1)) / (SR - SL);
  }
  return Fhat1;
}