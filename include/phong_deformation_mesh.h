#ifndef PHONG_DEF_MESH_H
#define PHONG_DEF_MESH_H

#include <Eigen/Core>

/**
 * Performs phong deformation on V_full
 *
 * @param V_full #V by 3 matrix of full, undeformed mesh geometry
 * @param tet_def #D by 3 matrix of deformed tetrahedral mesh geometry
 * @param tet_undef #D by 3 matrix of undeformed tetrahedral mesh geometry
 * @param tet_F #F by 4 matrix of tetrahedral mesh topology
 * @param E #V by 1 vector mapping vertices in V to a face in the tet mesh.
 * @param V_def OUTPUT: #V by 3 matrix of deformed (full) mesh geometry
 */
void phong_deformation_point(
    const Eigen::MatrixX3d &V_full,
    const Eigen::MatrixX3d &tet_def,
    const Eigen::MatrixX3d &tet_undef,
    const Eigen::MatrixX4i &tet_F,
    const Eigen::VectorXi &E,
    Eigen::MatrixX3d &V_def);

#endif
