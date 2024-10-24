#pragma once
#include "Array.h"
#include "parameters.h"

class InertiaEnergy {
public:
    static double val(const VecX &x, const VecX &x_tilde);
    static VecX grad(const VecX &x, const VecX &x_tilde);
    static void hess(const VecX &x, const VecX &x_tilde, Array<int> &I, Array<int> &J, Array<double> &V);
};
