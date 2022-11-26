#include "phong_deformation.h"

Eigen::Matrix3d edge_basis(const Eigen::MatrixXd &V, const Eigen::RowVector3d e)
{
  return V.row(e(1)) - V.row(e(0)), V.row(e(2)) - V.row(e(0)), V.row(e(3)) - V.row(e(0)); // todo T?
}

void phong_deformation(
    const Eigen::MatrixXd &V_full, // full undeformed geom
    const Eigen::RowVectorXi &e,   // tet verteces into V_full
    const Eigen::MatrixXd &D,      // deformed tet geom
    const Eigen::RowVectorXi &eD,  // tet vertices into D
    const Eigen::Vector3d &X,      // undeformed query point
    Eigen::Vector3d &x)
{
  x = X;

  // undeformed edge basis matrix
  Eigen::Matrix3d V_undef = edge_basis(V_full, e);

  Eigen::Vector3d bary = V_undef.inverse() * (X - V_full.row(e(0))); // todo T?
  double B_0 = 1 - bary(0) - bary(1) - bary(2);

  Eigen::Vector4d B; // barycentric coordinates of tet
  B << B_0, bary(0), bary(1), bary(2);

  Eigen::Vector3d def_tet_centroid = (D.row(eD(1)) + D.row(eD(2)) + D.row(eD(3)) + D.row(eD(0))) / 4.0;             // centroid of deformed tet
  Eigen::Vector3d tet_centroid = (V_full.row(e(1)) + V_full.row(e(2)) + V_full.row(e(3)) + V_full.row(e(0))) / 4.0; // centroid of undeformed tet

  // deformed edge basis matrix;
  Eigen::Matrix3d v_def = edge_basis(D, eD);

  // deformation gradient (centroid)
  Eigen::Matrix3d F_c = v_def * V_undef.inverse();

  // linear interpolation of deformation
  Eigen::Vector3d f_c = def_tet_centroid + F_c * (X - tet_centroid);
}