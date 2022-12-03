#include "find_tet.h"
#include <igl/barycentric_coordinates.h>
#include <limits>
#include <iostream>

void find_tet(
    const Eigen::MatrixX3d &V_full,
    const Eigen::MatrixX3d &tet_undef,
    const Eigen::MatrixX4i &tet_F,
    Eigen::VectorXi &E)
{
  E.resize(V_full.rows());

  int bad_points = 0; // how many points outside tets

  // loop over all vertices in V_full
  for (int i = 0; i < V_full.rows(); i++)
  {
    Eigen::MatrixXd p, a, b, c, d;
    Eigen::MatrixXd bary;
    Eigen::Vector4i face;

    p = V_full.row(i);

    int t = -1;

    // loop over all tetrahedra in tet_F
    for (int j = 0; j < tet_F.rows(); j++)
    {
      face = tet_F.row(j);
      a = tet_undef.row(face[0]);
      b = tet_undef.row(face[1]);
      c = tet_undef.row(face[2]);
      d = tet_undef.row(face[3]);

      igl::barycentric_coordinates(p, a, b, c, d, bary);

      if (bary.minCoeff() >= 0)
      {
        // point is inside tet when all barycentric coordinates are positive
        t = j;
        break;
      }
    }

    if (t == -1)
    {
      bad_points++;
      E(i) = 0;
    }
    else
    {
      E(i) = t;
    }
  }

  assert(bad_points == 0 && "Mesh contains vertices outside tetrahedra");
}