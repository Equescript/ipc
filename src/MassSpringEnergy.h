#pragma once
#include "Array.h"
#include "parameters.h"

class MassSpringEnergy {
public:
    static double val(const VecX &x, const Array<Vec2i> &e);
    static VecX grad(const VecX &x, const Array<Vec2i> &e);
    static void hess(const VecX &x, const Array<Vec2i> &e, Array<int> &I, Array<int> &J, Array<double> &V);
};
