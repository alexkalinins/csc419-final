#ifndef PHONG_DEF_H
#define PHONG_DEF_H

#include <Eigen/Core>

/**
 * Consider a tetrahedron with undeformed (material space) vertex positions, given an arbitrary interior point,
 * we compute its displaced position using Phong Deformation Model.
 *
 * @param V #V by 3 list of undeformed vertex positions (material space) Each row is a single undefored vertex position
 * @param e 1 by 4 vertex indices for this tetrahedron
 * @param D #D by 3 list of deformed vertex positions from the simpler tetrahadral mesh
 * @param eD 1 by 4 vertex indices for this tetrahedron in the deformed mesh
 * @param X 3 by 1 undeformed query point position inside this tetrahedron
 *
 * @param x OUTPUT: 3 by 1 deformed query point position
 */
void phong_deformation(
    const Eigen::Vector3d &V,
    const Eigen::RowVectorXi &e,
    const Eigen::MatrixXd &D,    
    const Eigen::RowVectorXi &eD,
    const Eigen::MatrixXd &X,
    Eigen::Vector3d &x);

#endif
