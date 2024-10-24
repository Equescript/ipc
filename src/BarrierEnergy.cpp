#include <math.h>
#include "Array.h"
#include "parameters.h"
#include "BarrierEnergy.h"

double BarrierEnergy::val(const VecX &x) {
    double sum = 0.0;
    for (int i=0; i<x.size()/2; i++) {
        double d = x[i*2+1] - _y_ground;
        if (d < _dhat) {
            double s = d / _dhat;
            sum += _dhat * _kappa / 2.0 * (s - 1.0) * log(s);
        }
        /* d = (Vec2(x[i*2], x[i*2+1]) - o).norm() - radius;
        if (d < _dhat) {
            double s = d / _dhat;
            sum += _dhat * _kappa / 2.0 * (s - 1.0) * log(s);
        } */
    }
    return sum;
}
VecX BarrierEnergy::grad(const VecX &x) {
    VecX g(x.size());
    for (int i=0; i<x.size()/2; i++) {
        double d = x[i*2+1] - _y_ground;
        if (d < _dhat) {
            double s = d / _dhat;
            g[i*2] = 0.0;
            g[i*2+1] = _dhat * (_kappa / 2.0 * (log(s) / _dhat + (s - 1.0) / d));
        }
        /* d = (Vec2(x[i*2], x[i*2+1]) - o).norm() - radius;
        if (d < _dhat) {
            double s = d / _dhat;
            g[i*2] += 0.0;
            g[i*2+1] += _dhat * (_kappa / 2.0 * (log(s) / _dhat + (s - 1.0) / d));
        } */
    }
    return g;
}
void BarrierEnergy::hess(const VecX &x, Array<int> &I, Array<int> &J, Array<double> &V) {
    for (int i=0; i<x.size()/2; i++) {
        I[i] = i * 2 + 1;
        J[i] = i * 2 + 1;
        double d = x[i*2+1] - _y_ground;
        if (d < _dhat) {
            double s = d / _dhat;
            V[i] = _dhat * _kappa / (2.0 * d * d * _dhat) * (d + _dhat);
        } else {
            V[i] = 0.0;
        }
        /* d = (Vec2(x[i*2], x[i*2+1]) - o).norm() - radius;
        if (d < _dhat) {
            double s = d / _dhat;
            V[i] += _dhat * _kappa / (2.0 * d * d * _dhat) * (d + _dhat);
        } else {
            V[i] += 0.0;
        } */
    }
}
double BarrierEnergy::init_step_size(const VecX &x, const VecX &p) {
    double alpha = 1.0;
    for (int i=0; i<x.size()/2; i++) {
        if (p[i*2+1] < 0) {
            // double d = (x[i*2+1] - _y_ground) + ((Vec2(x[i*2], x[i*2+1]) - o).norm() - radius);
            double d = x[i*2+1] - _y_ground;
            alpha = std::min(alpha, 0.9 * -d / p[i*2+1]);
        }
    }
    return alpha;
}

