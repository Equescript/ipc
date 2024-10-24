#include "Array.h"
#include "parameters.h"
#include "InertiaEnergy.h"

double InertiaEnergy::val(const VecX &x, const VecX &x_tilde) {
    double sum = 0.0;
    for (int i=0; i<x.size()/2; i++) {
        Vec2 diff(x[i*2] - x_tilde[i*2], x[i*2+1] - x_tilde[i*2+1]);
        sum += 0.5 * _m[i] * diff.dot(diff);
    }
    return sum;
}
VecX InertiaEnergy::grad(const VecX &x, const VecX &x_tilde) {
    VecX g(x.size());
    for (int i=0; i<x.size()/2; i++) {
        g[i*2] = _m[i] * (x[i*2] - x_tilde[i*2]);
        g[i*2+1] = _m[i] * (x[i*2+1] - x_tilde[i*2+1]);
    }
    return g;
}
void InertiaEnergy::hess(const VecX &x, const VecX &x_tilde, Array<int> &I, Array<int> &J, Array<double> &V) {
    for (int i=0; i<x.size()/2; i++) {
        for (int d=0; d<2; d++) {
            I[i * 2 + d] = i * 2 + d;
            J[i * 2 + d] = i * 2 + d;
            V[i * 2 + d] = _m[i];
        }
    }
}

