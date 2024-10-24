#ifdef _WIN32
#define EXPORT extern "C" __declspec( dllexport )
#else
#define EXPORT extern "C" __attribute__((visibility ("default")))
#endif

#include <iostream>
#include "Array.h"
#include "parameters.h"
#include "BarrierEnergy.h"
#include "GravityEnergy.h"
#include "InertiaEnergy.h"
#include "MassSpringEnergy.h"

double h = 0.01;        // time step size in s

double _dhat = 0;
double _kappa = 0;
double _y_ground = 0;

Vec2 _gravity(0.0, -9.81);

Array<double> _m;
Array<double> _l2; // rest length squared
Array<double> _k; // spring stiffness

Vec2 o(0.0, 0.0);
double radius = 0.3;

class Scene {
public:
    Array<Vec2> x;
    Array<Vec2i> e;
    Scene() {}
    Scene(Array<Vec2> x, Array<Vec2i> e): x(x), e(e) {}
};

Scene SCENE;
VecX X;
VecX V;

// Array<Vec2> InertiaEnergy_g;
// Array<Vec2> MassSpringEnergy_g;
// Array<Vec2> GravityEnergy_g;
// Array<Vec2> BarrierEnergy_g;
// Array<Vec2> IP_g;

Array<int> I_In;
Array<int> J_In;
Array<double> V_In;
Array<int> I_MS;
Array<int> J_MS;
Array<double> V_MS;
Array<int> I_B;
Array<int> J_B;
Array<double> V_B;
Array<int> I_H;
Array<int> J_H;
Array<double> V_H;

EXPORT void init(Vec2* x_data, int x_length, Vec2i* e_data, int e_length, double* m_data, double* l2_data, double* k_data, double dhat, double kappa) {
    SCENE = Scene(
        Array<Vec2>(x_data, x_length),
        Array<Vec2i>(e_data, e_length)
    );
    _m = Array<double>(m_data, SCENE.x.length);
    _l2 = Array<double>(l2_data, SCENE.e.length);
    _k = Array<double>(k_data, SCENE.e.length);
    _dhat = dhat;
    _kappa = kappa;
    /* _contact_area = Array<double>(new double[SCENE.x.length], SCENE.x.length);
    double step = _side_len / _n_seg;
    for (int i=0; i<SCENE.x.length; i++) {
        _contact_area[i] = step;
    } */
    // InertiaEnergy_g = Array<Vec2>(new Vec2[SCENE.x.length], SCENE.x.length);
    // MassSpringEnergy_g = Array<Vec2>(new Vec2[SCENE.x.length], SCENE.x.length);
    // GravityEnergy_g = Array<Vec2>(new Vec2[SCENE.x.length], SCENE.x.length);
    // BarrierEnergy_g = Array<Vec2>(new Vec2[SCENE.x.length], SCENE.x.length);
    // IP_g = Array<Vec2>(new Vec2[SCENE.x.length], SCENE.x.length);
    int H_length = SCENE.x.length * 2 + SCENE.e.length * 16 + SCENE.x.length;
    I_In = Array<int>(new int[H_length], SCENE.x.length * 2);
    J_In = Array<int>(new int[H_length], SCENE.x.length * 2);
    V_In = Array<double>(new double[H_length], SCENE.x.length * 2);
    I_MS = Array<int>(I_In.data + SCENE.x.length * 2, SCENE.e.length * 16);
    J_MS = Array<int>(J_In.data + SCENE.x.length * 2, SCENE.e.length * 16);
    V_MS = Array<double>(V_In.data + SCENE.x.length * 2, SCENE.e.length * 16);
    I_B = Array<int>(I_In.data + SCENE.x.length * 2 + SCENE.e.length * 16, SCENE.x.length);
    J_B = Array<int>(J_In.data + SCENE.x.length * 2 + SCENE.e.length * 16, SCENE.x.length);
    V_B = Array<double>(V_In.data + SCENE.x.length * 2 + SCENE.e.length * 16, SCENE.x.length);
    I_H = Array<int>(I_In.data, H_length);
    J_H = Array<int>(J_In.data, H_length);
    V_H = Array<double>(V_In.data, H_length);
    X = VecX(SCENE.x.length*2);
    V = VecX(SCENE.x.length*2);
    memcpy(&X(0), SCENE.x.data, sizeof(double)*SCENE.x.length*2);
    memset(&V(0), 0, sizeof(double)*SCENE.x.length*2);
}

double IP_val(const VecX &x, const VecX &x_tilde, const Array<Vec2i> &e) {
    return InertiaEnergy::val(x, x_tilde) + h * h * (MassSpringEnergy::val(x, e) + GravityEnergy::val(x) + BarrierEnergy::val(x));
}

VecX IP_grad(const VecX &x, const VecX &x_tilde, const Array<Vec2i> &e) {
    return InertiaEnergy::grad(x, x_tilde) + h * h * (MassSpringEnergy::grad(x, e) + GravityEnergy::grad(x) + BarrierEnergy::grad(x));
}

typedef Eigen::SparseMatrix<double, Eigen::ColMajor> SparseMatrix;

SparseMatrix IP_hess(const VecX &x, const VecX &x_tilde, const Array<Vec2i> &e) {
    InertiaEnergy::hess(x, x_tilde, I_In, J_In, V_In);
    MassSpringEnergy::hess(x, e, I_MS, J_MS, V_MS);
    BarrierEnergy::hess(x, I_B, J_B, V_B);
    for (int i=0; i<V_MS.length; i++) {
        V_MS[i] *= h * h;
    }
    for (int i=0; i<V_B.length; i++) {
        V_B[i] *= h * h;
    }
    SparseMatrix H(x.size(), x.size());
    for (int i=0; i<I_H.length; i++) {
        H.insert(I_H[i], J_H[i]) = V_H[i];
    }
    H.makeCompressed();
    return H;
}

VecX search_dir(const VecX &x, const VecX &x_tilde, const Array<Vec2i> &e) {
    SparseMatrix projected_hess = IP_hess(x, x_tilde, e);
    VecX reshaped_grad = IP_grad(x, x_tilde, e);
    // eliminate DOF by modifying gradient and Hessian for DBC:
    Eigen::ConjugateGradient<SparseMatrix> solver;
    solver.compute(projected_hess);
    return solver.solve(-reshaped_grad);
}

void step_forward(VecX &x, const Array<Vec2i> &e, VecX &v) {
    VecX x_tilde = x + h * v;
    VecX x_n = x;
    // Newton loop
    int iter = 0;
    double E_last = IP_val(x, x_tilde, e);
    VecX p = search_dir(x, x_tilde, e);
    double tol = 1e-2;
    while (p.lpNorm<Eigen::Infinity>() / h > tol) {
        std::cout << "Iteration" << iter << ":" << std::endl;
        std::cout << "residual =" << p.lpNorm<Eigen::Infinity>() / h << std::endl;
        // filter line search
        double alpha = BarrierEnergy::init_step_size(x, p);  // avoid interpenetration and tunneling
        // double alpha = 1.0;
        while (IP_val(x+alpha*p, x_tilde, e) > E_last) {
            alpha /= 2.0;
        }
        std::cout << "step size =" << alpha << std::endl;

        x += alpha * p;
        E_last = IP_val(x, x_tilde, e);
        p = search_dir(x, x_tilde, e);
        iter++;
    }
    v = (x - x_n) / h;   // implicit Euler velocity update
}

EXPORT void step() {
    step_forward(X, SCENE.e, V);
    memcpy(SCENE.x.data, &X(0), sizeof(double)*SCENE.x.length*2);
}
