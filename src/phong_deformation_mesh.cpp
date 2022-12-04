#include "phong_deformation_mesh.h"
#include "vertex_def_gradient.h"
#include "cell_def_gradient.h"
#include <vector>
#include <igl/barycentric_coordinates.h>
#include <Eigen/Dense>

void phong_deformation_mesh(
    const Eigen::MatrixX3d &V_full,
    const Eigen::MatrixX3d &tet_def,
    const Eigen::MatrixX3d &tet_undef,
    const Eigen::MatrixX4i &tet_F,
    const Eigen::VectorXi &E,
    Eigen::MatrixX3d &V_def)
{

  V_def.resize(V_full.rows(), 3);
  V_def.setZero();

  Eigen::MatrixX3d c_def(tet_F.rows(), 3), c_undef(tet_F.rows(), 3);

  c_def.setZero();
  c_undef.setZero();

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
  vertex_def_gradient(tet_def, tet_F, c_undef, F_cell, F_vert);

  Eigen::MatrixX4d B(V_full.rows(), 4); // barycentric coordinates
  B.setZero();

  // compute barycentric coordinates
  for (int i = 0; i < V_full.rows(); i++)
  {
    Eigen::RowVector4d bary;
    int face = E(i);

    igl::barycentric_coordinates(V_full.row(i), tet_undef.row(tet_F(face, 0)), tet_undef.row(tet_F(face, 1)), tet_undef.row(tet_F(face, 2)), tet_undef.row(tet_F(face, 3)), bary);
    B.row(i) = bary;
  }

  Eigen::Vector3d f_c, f_v;

  // compute deformed vertices
  for (int i = 0; i < V_full.rows(); i++)
  {
    
    f_v.setZero();
    f_c.setZero();

    int tet = E(i);

    Eigen::Vector3d diff_cent = V_full.row(i) - c_undef.row(tet);

    // cell deformation
    f_c = (F_cell[tet] * diff_cent) + c_def.row(E(i)).transpose();
    
    // vertex deformation
    Eigen::Vector3d diff_vert;

    for (int j = 0; j < 4; j++)
    {
      int v = tet_F(tet, j);
      diff_vert = V_full.row(i) - tet_undef.row(v);

      f_v += B(i, j) * (F_vert[v] * diff_vert + tet_def.row(v).transpose());
    }

    // phong deformation
    V_def.row(i) = (f_c + f_v) / 2;
  }
}