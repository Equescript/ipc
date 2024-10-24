#pragma once
#include "Array.h"
#include "parameters.h"

class BarrierEnergy {
public:
    static double val(const VecX &x);
    static VecX grad(const VecX &x);
    static void hess(const VecX &x, Array<int> &I, Array<int> &J, Array<double> &V);
    static double init_step_size(const VecX &x, const VecX &p);
};
