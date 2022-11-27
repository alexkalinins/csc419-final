#ifndef VERT_DEF_GRAD_H
#define VERT_DEF_GRAD_H

#include <Eigen/Core>
#include <vector>

/**
 * Computes the cell deformation gradients for a deformed tetrahedral mesh.
 *
 * @param V_undef #V by 3 deformed mesh geometry
 * @param F #F by 4 tetrahedral mesh topology
 * @param C #F by 3 list of cell centroids (deformed)
 * @param F_cell #F list of cell gradients
 * @param grad OUTPUT: #V long list of 3 by 3 vector deformation gradients;
 * one for each vector in V
 */
void vertex_def_gradient(
    const Eigen::MatrixX3d &V_undef,
    const Eigen::MatrixX4i &F,
    const Eigen::MatrixX3d &C,
    const std::vector<Eigen::Matrix3d> &F_cell,
    std::vector<Eigen::Matrix3d> &grad);

#endif
