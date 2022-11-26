#include "phong_deformation.h"

void phong_deformation(
    const Eigen::Vector3d &V,
    const Eigen::RowVectorXi &e,
    const Eigen::MatrixXd &D,
    const Eigen::RowVectorXi &eD,
    const Eigen::MatrixXd &X,
    Eigen::Vector3d &x)
{
  x = V;
}