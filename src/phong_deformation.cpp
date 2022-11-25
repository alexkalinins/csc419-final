#include "phong_deformation.h"

void phong_deformation(
    const Eigen::MatrixXd &X, // undeformed complex mesh geometry
    const Eigen::MatrixXd &V, // deformed tetrahedral mesh geometry
    const Eigen::MatrixXi &F, // deformed tetrahedral mesh topology
    Eigen::MatrixXd &x        // deformed X mesh geometry
)
{
  x = X;
}