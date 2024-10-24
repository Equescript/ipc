#pragma once
#include <Eigen/Eigen>
#include "Array.h"

typedef Eigen::Vector2d Vec2;
typedef Eigen::Vector2i Vec2i;
typedef Eigen::VectorXd VecX;
typedef Eigen::Matrix2d Mat2;

extern double h;

extern double _dhat;
extern double _kappa;
extern double _y_ground;

extern Vec2 _gravity;

extern Array<double> _m;
extern Array<double> _l2; // rest length squared
extern Array<double> _k; // spring stiffness
