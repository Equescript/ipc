#pragma once
#include "Array.h"
#include "parameters.h"

class GravityEnergy {
public:
    static double val(const VecX &x);
    static VecX grad(const VecX &x);
};
