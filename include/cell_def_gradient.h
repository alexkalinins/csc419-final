#ifndef CELL_DEF_GRAD_H
#define CELL_DEF_GRAD_H

#include <Eigen/Core>
#include <vector>

/**
 * Computes the cell deformation gradients for a deformed tetrahedral mesh.
 *
 * @param V_undef #V by 3 undeformed mesh geometry
 * @param V_def #V by 3 deformed mesh geometry
 * @param F #F by 4 tetrahedral mesh topology
 * @param grad OUTPUT: #F long list of 3 by 3 cell deformation gradients; one for each face in F
 */
void cell_def_gradient(
    const Eigen::MatrixXd &V_undef,
    const Eigen::MatrixXd &V_def,
    const Eigen::MatrixXi &F,
    std::vector<Eigen::Matrix3d> &grad);

#endif
