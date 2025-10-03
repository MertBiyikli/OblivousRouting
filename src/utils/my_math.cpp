//
// Created by Mert Biyikli on 30.09.25.
//

#include "my_math.h"
#include <cmath>
#include <stdexcept>

 double NegativeExponent(double base, int exp) {
    if (base == 0.0) {
        throw std::invalid_argument("Base cannot be zero for negative exponent.");
    }
    if (exp < 0) {
        throw std::invalid_argument("Exponent must be non-negative.");
    }
    return 1.0 / std::pow(base, exp);
}