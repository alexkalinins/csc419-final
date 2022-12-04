#include "cell_def_gradient.h"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <vector>
#include <iostream>

Eigen::Matrix3d edge_basis(const Eigen::MatrixX3d &V, const Eigen::Vector4i &e)
{
  Eigen::Matrix3d basis;
  basis.col(0) = (V.row(e(1)) - V.row(e(0))).transpose();
  basis.col(1) = (V.row(e(2)) - V.row(e(0))).transpose();
  basis.col(2) = (V.row(e(3)) - V.row(e(0))).transpose();
  return basis;
}

void cell_def_gradient(
    const Eigen::MatrixX3d &V_undef,
    const Eigen::MatrixX3d &V_def,
    const Eigen::MatrixX4i &F,
    std::vector<Eigen::Matrix3d> &grad)
{
  grad.resize(F.rows());

  // tet face:
  Eigen::Vector4i e;

  // deformed, undeformed edge basis matrices
  Eigen::Matrix3d mat_def, mat_undef;

  for (int i = 0; i < F.rows(); i++)
  {
    e = F.row(i);

    mat_undef = edge_basis(V_undef, e);
    mat_def = edge_basis(V_def, e);

    Eigen::Matrix3d F = mat_def * mat_undef.inverse();
    grad[i] = F;
  }
}