#include "cell_def_gradient.h"
#include <Eigen/Core>
#include <vector>

Eigen::Matrix3d edge_basis(const Eigen::MatrixXd &V, const Eigen::RowVector3d e)
{
  return V.row(e(1)) - V.row(e(0)), V.row(e(2)) - V.row(e(0)), V.row(e(3)) - V.row(e(0)); // todo T?
}

void cell_def_gradient(
    const Eigen::MatrixXd &V_undef,
    const Eigen::MatrixXd &V_def,
    const Eigen::MatrixXi &F,
    std::vector<Eigen::Matrix3d> &grad)
{
  grad.resize(F.rows());

  // tet face:
  Eigen::RowVector3d e;

  // deformed, undeformed edge basis matrices
  Eigen::Matrix3d mat_def, mat_undef;

  for (int i = 0; i < F.rows(); i++)
  {
    e = F.row(i);
    mat_undef = edge_basis(V_undef, e);
    mat_def = edge_basis(V_def, e);

    Eigen::Matrix3d F = mat_def * mat_undef.inverse();
    grad.push_back(F);
  }
}