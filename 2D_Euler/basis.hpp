#pragma once
#include "global.hpp"
#include <sys/stat.h>

template <int T> struct poly {
  // static constexpr float value(const float &x, const float &y) { return 0; }
};

template <> struct poly<0> {
  static constexpr double value(const double &x, const double &y) {
    return 1.0;
  }
};

template <> struct poly<1> {
  static constexpr double value(const double &x, const double &y) { return x; }
};

template <> struct poly<2> {
  static constexpr double value(const double &x, const double &y) { return y; }
};

template <> struct poly<3> {
  static constexpr double value(const double &x, const double &y) {
    return x * x - 1.0 / 3.0;
  }
};

template <> struct poly<4> {
  static constexpr double value(const double &x, const double &y) {
    return x * y;
  }
};

template <> struct poly<5> {
  static constexpr double value(const double &x, const double &y) {
    return y * y - 1.0 / 3.0;
  }
};

template <int T> struct polyx {};

template <> struct polyx<0> {
  static constexpr double value(const double &x, const double &y) {
    return 0.0;
  }
};

template <> struct polyx<1> {
  static constexpr double value(const double &x, const double &y) {
    return 1.0 / hx1;
  }
};

template <> struct polyx<2> {
  static constexpr double value(const double &x, const double &y) {
    return 0.0;
  }
};

template <> struct polyx<3> {
  static constexpr double value(const double &x, const double &y) {
    return 2.0 * x / hx1;
  }
};

template <> struct polyx<4> {
  static constexpr double value(const double &x, const double &y) {
    return y / hx1;
  }
};

template <> struct polyx<5> {
  static constexpr double value(const double &x, const double &y) {
    return 0.0;
  }
};

template <int T> struct polyy {};

template <> struct polyy<0> {
  static constexpr double value(const double &x, const double &y) {
    return 0.0;
  }
};

template <> struct polyy<1> {
  static constexpr double value(const double &x, const double &y) {
    return 0.0;
  }
};

template <> struct polyy<2> {
  static constexpr double value(const double &x, const double &y) {
    return 1.0 / hy1;
  }
};

template <> struct polyy<3> {
  static constexpr double value(const double &x, const double &y) {
    return 0.0;
  }
};

template <> struct polyy<4> {
  static constexpr double value(const double &x, const double &y) {
    return x / hy1;
  }
};

template <> struct polyy<5> {
  static constexpr double value(const double &x, const double &y) {
    return 2.0 * y / hy1;
  }
};

inline double Poly(const double &x, const double &y, const int &k) {
  if (k == 0) {
    return 1.0;
  } else if (k == 1) {
    return x;
  } else if (k == 2) {
    return y;
  } else if (k == 3) {
    return x * x - 1.0 / 3.0;
  } else if (k == 4) {
    return x * y;
  } else {
    return y * y - 1.0 / 3.0;
  }
}

inline double PolyX(const double &x, const double &y, const int &k) {
  if (k == 0) {
    return 0.0;
  } else if (k == 1) {
    return 1.0 / hx1;
  } else if (k == 2) {
    return 0.0;
  } else if (k == 3) {
    return 2.0 * x / hx1;
  } else if (k == 4) {
    return y / hx1;
  } else {
    return 0.0;
  }
}

inline double PolyY(const double &x, const double &y, const int &k) {
  if (k == 0) {
    return 0.0;
  } else if (k == 1) {
    return 0.0;
  } else if (k == 2) {
    return 1.0 / hy1;
  } else if (k == 3) {
    return 0.0;
  } else if (k == 4) {
    return x / hy1;
  } else {
    return 2.0 * y / hy1;
  }
}