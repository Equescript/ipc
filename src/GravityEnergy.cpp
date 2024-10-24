#include "Array.h"
#include "parameters.h"
#include "GravityEnergy.h"

double GravityEnergy::val(const VecX &x) {
    double sum = 0.0;
    for (int i=0; i<x.size()/2; i++) {
        sum += -_m[i] * Vec2(x[i*2], x[i*2+1]).dot(_gravity);
    }
    return sum;
}
VecX GravityEnergy::grad(const VecX &x) {
    VecX g(x.size());
    for (int i=0; i<x.size()/2; i++) {
        g[i*2] = -_m[i] * _gravity[0];
        g[i*2+1] = -_m[i] * _gravity[1];
    }
    return g;
}
