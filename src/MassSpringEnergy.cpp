#include "Array.h"
#include "parameters.h"
#include "MassSpringEnergy.h"

Mat2 outer_product(Vec2 &v, Vec2 &w) {
    Mat2 m;
    m << (v[0] * w[0]), (v[0] * w[1]),
         (v[1] * w[0]), (v[1] * w[1]);
    return m;
}

/* template<const int N>
Eigen::Matrix<double, N, N> make_PSD(Eigen::Matrix<double, N, N> &m) {
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, N, N>> solver(m);
    Eigen::Vector<double, N> lam = solver.eigenvalues();
    Eigen::Matrix<double, N, N> V = solver.eigenvectors();
    lam = (lam.array() > 0.0).select(lam, 0.0);
    return V * lam.asDiagonal() * V.transpose();
} */

Eigen::Matrix4d make_PSD(Eigen::Matrix4d &m) {
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix4d> solver(m);
    Eigen::Vector4d lam = solver.eigenvalues();
    Eigen::Matrix4d V = solver.eigenvectors();
    lam = (lam.array() > 0.0).select(lam, 0.0);
    return V * lam.asDiagonal() * V.transpose();
}


double MassSpringEnergy::val(const VecX &x, const Array<Vec2i> &e) {
    double sum = 0.0;
    for (int i=0; i<e.length; i++) {
        Vec2 diff(x[e[i][0]*2] - x[e[i][1]*2], x[e[i][0]*2+1] - x[e[i][1]*2+1]);
        double ratio = diff.dot(diff) / _l2[i] - 1.0;
        sum += _l2[i] * 0.5 * _k[i] * ratio * ratio;
    }
    return sum;
}
VecX MassSpringEnergy::grad(const VecX &x, const Array<Vec2i> &e) {
    VecX g(x.size());
    for (int i=0; i<e.length; i++) {
        Vec2 diff(x[e[i][0]*2] - x[e[i][1]*2], x[e[i][0]*2+1] - x[e[i][1]*2+1]);
        Vec2 g_diff = 2.0 * _k[i] * (diff.dot(diff) / _l2[i] - 1.0) * diff;
        g[e[i][0]*2] += g_diff[0];
        g[e[i][0]*2+1] += g_diff[1];
        g[e[i][1]*2] -= g_diff[0];
        g[e[i][1]*2+1] -= g_diff[1];
    }
    return g;
}
void MassSpringEnergy::hess(const VecX &x, const Array<Vec2i> &e, Array<int> &I, Array<int> &J, Array<double> &V) {
    assert(I.length == e.length * 16);
    for (int i=0; i<e.length; i++) {
        Vec2 diff(x[e[i][0]*2] - x[e[i][1]*2], x[e[i][0]*2+1] - x[e[i][1]*2+1]);
        Mat2 H_diff = 2.0 * _k[i] / _l2[i] * (2.0 * outer_product(diff, diff) + (diff.dot(diff) - _l2[i]) * Mat2::Identity());
        Eigen::Matrix4d symmetric_matrix;
        symmetric_matrix.topLeftCorner(2, 2) = H_diff;
        symmetric_matrix.topRightCorner(2, 2) = -H_diff;
        symmetric_matrix.bottomLeftCorner(2, 2) = -H_diff;
        symmetric_matrix.bottomRightCorner(2, 2) = H_diff;
        Eigen::Matrix4d H_local = make_PSD(symmetric_matrix);
        for (int nI=0; nI<2; nI++) {
            for (int nJ=0; nJ<2; nJ++) {
                int indStart = i * 16 + (nI * 2 + nJ) * 4;
                for (int r=0; r<2; r++) {
                    for (int c=0; c<2; c++) {
                        I[indStart + r * 2 + c] = e[i][nI] * 2 + r;
                        J[indStart + r * 2 + c] = e[i][nJ] * 2 + c;
                        V[indStart + r * 2 + c] = H_local(nI * 2 + r, nJ * 2 + c);
                    }
                }
            }
        }
    }
}
