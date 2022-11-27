#include "phong_deformation_mesh.h"
#include "vertex_def_gradient.h"
#include "cell_def_gradient.h"
#include <vector>
#include <igl/barycentric_coordinates.h>

Eigen::Matrix3d edge_basis(const Eigen::MatrixXd &V, const Eigen::RowVector3d e)
{
  return V.row(e(1)) - V.row(e(0)), V.row(e(2)) - V.row(e(0)), V.row(e(3)) - V.row(e(0)); // todo T?
}

void phong_deformation_point(
    const Eigen::MatrixXd &V_full,
    const Eigen::MatrixXd &tet_def,
    const Eigen::MatrixXd &tet_undef,
    const Eigen::MatrixXi &tet_F,
    const Eigen::VectorXi &E,
    Eigen::Vector3d &V_def)
{

  V_def.resize(V_full.size());

  Eigen::MatrixXd c_def(tet_F.rows(), 3), c_undef(tet_F.rows(), 3);

  // compute deformed and undeformed tet-mesh centroids
  for (int i = 0; i < tet_F.rows(); i++)
  {
    c_def.row(i) = (tet_def.row(tet_F(i, 0)) + tet_def.row(tet_F(i, 1)) + tet_def.row(tet_F(i, 2)) + tet_def.row(tet_F(i, 3))) / 4;
    c_undef.row(i) = (tet_undef.row(tet_F(i, 0)) + tet_undef.row(tet_F(i, 1)) + tet_undef.row(tet_F(i, 2)) + tet_undef.row(tet_F(i, 3))) / 4;
  }

  // compute cell gradients for each tet mesh face
  std::vector<Eigen::Matrix3d> F_cell;
  cell_def_gradient(tet_undef, tet_def, tet_F, F_cell);

  // compute vertex gradients for each vertex in tet mesh
  std::vector<Eigen::Matrix3d> F_vert;
  vertex_def_gradient(tet_undef, tet_F, c_def, F_cell, F_vert);

  Eigen::MatrixXd B(V_full.rows(), 4); // barycentric coordinates

  // compute barycentric coordinates
  for (int i = 0; i < V_full.rows(); i++)
  {
    Eigen::Vector4d bary;
    igl::barycentric_coordinates(V_full.row(i), tet_undef.row(E(i, 0)), tet_undef.row(E(i, 1)), tet_undef.row(E(i, 2)), tet_undef.row(E(i, 3)), bary);
    B.row(i) = bary;
  }

  Eigen::Vector3d f_c, f_v;

  // compute deformed vertices
  for (int i = 0; i < V_full.rows(); i++)
  {
    // cell deformation
    f_c = c_def.row(E(i)) + F_cell[E(i)] * (V_full.row(i) - c_undef.row(E(i)));

    // vertex deformation
    f_v.setZero();

    for (int j = 0; j < 4; j++)
    {
      f_v += B(i, j) * F_vert[E(i, j)] * (V_full.row(i) - tet_undef.row(E(i, j)));
    }

    // phong deformation
    V_def.row(i) = (f_c + f_v) / 2;
  }
}