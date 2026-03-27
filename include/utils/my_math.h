//
// Created by Mert Biyikli on 30.09.25.
//

#ifndef OBLIVIOUSROUTING_MY_MATH_H
#define OBLIVIOUSROUTING_MY_MATH_H
#include <cmath>
#include <stdexcept>

constexpr static double EPS = 1e-16;
constexpr static double SOFT_EPS = 1e-6;
constexpr static double VERY_SOFT_EPS = 1e-1;

inline double NegativeExponent(double base, int exp) {
    if (base == 0.0) {
        throw std::invalid_argument("Base cannot be zero for negative exponent.");
    }
    if (exp < 0) {
        throw std::invalid_argument("Exponent must be non-negative.");
    }
    return 1.0 / std::pow(base, exp);
}

#endif //OBLIVIOUSROUTING_MY_MATH_H