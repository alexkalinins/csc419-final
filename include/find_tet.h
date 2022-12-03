#ifndef FIND_TET_H
#define FIND_TET_H

#include <Eigen/Core>

/**
 * For each point in V_full, this function determines in which tetrahedron
 * it is contained
 *
 * @param V_full triangle mesh
 * @param tet_undef undeformed tetrahedron mesh
 * @param tet_F tetrahedron mesh topology
 * @param E OUTPUT #V by 1 vector mapping vertices in V to a face in the tet mesh.
 */
void find_tet(
    const Eigen::MatrixX3d &V_full,
    const Eigen::MatrixX3d &tet_undef,
    const Eigen::MatrixX4i &tet_F,
    Eigen::VectorXi &E);

#endif
