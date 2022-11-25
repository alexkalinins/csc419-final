#ifndef PHONG_DEF_H
#define PHONG_DEF_H

#include <Eigen/Core>

/**
 * Phong deformation is used to cheaply compute deformations on complex meshes. Physics simulations
 * and other deformations can be performed on a simpler tetrahedral mesh. Given such simpler deformed
 * mesh and the complex undeformed object (in material space; not necessarily tetrahedral), this
 * function applies the deformation from the simpler mesh to the complex mesh.
 *
 * This allows physics simulations to be less expensive, as they are performed on a simpler
 * tetrahedral mesh.
 *
 * @param X Undeformed complex mesh geometry (in material space)
 * @param V Deformed tetrahedral mesh geometry (in material space)
 * @param F Tetrahedral mesh faces (indeces into V)
 * @param x OUTPUT: Deformed complex mesh geometry (in material space)
 */
void phong_deformation(
    const Eigen::MatrixXd &X, // undeformed complex mesh geometry
    const Eigen::MatrixXd &V, // deformed tetrahedral mesh geometry
    const Eigen::MatrixXi &F, // deformed tetrahedral mesh topology
    Eigen::MatrixXd &x  // deformed X mesh geometry
);

#endif
