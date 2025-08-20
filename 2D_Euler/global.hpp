#pragma once

#include <array>
constexpr int NX = 480;
constexpr int NY = 480;
constexpr int NumEq = 4;
constexpr int NumGLP = 5;
constexpr int k = 2;
constexpr int dimPk = (k + 1) * (k + 2) / 2;

constexpr double gasGamma = 1.4;

constexpr double xb = 1.0;
constexpr double xa = 0.0;
constexpr double ya = 0.0;
constexpr double yb = 1.0;
constexpr double hx = (xb - xa) / NX;
constexpr double hy = (yb - ya) / NY;
constexpr double hx1 = 0.5 * hx;
constexpr double hy1 = 0.5 * hy;
constexpr double M = 100.0;
constexpr double CFL = 0.2;

constexpr double tend = 0.8;

constexpr std::array<double, NumGLP> lambda = {
    -0.9061798459386639927976269, -0.5384693101056830910363144, 0.0,
    0.5384693101056830910363144, 0.9061798459386639927976269};

constexpr std::array<double, NumGLP> weight = {
    0.2369268850561890875142640, 0.4786286704993664680412915,
    0.5688888888888888888888889, 0.4786286704993664680412915,
    0.2369268850561890875142640};

constexpr std::array<double, NumGLP> lambdaL = {
    -1.0, -0.6546536707079771437983, 0.0, 0.654653670707977143798, 1.0};

constexpr std::array<double, dimPk> mass = {1.0,        1.0 / 3.0, 1.0 / 3.0,
                                            4.0 / 45.0, 1.0 / 9.0, 4.0 / 45.0};