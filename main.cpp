#include <igl/min_quad_with_fixed.h>
#include <igl/read_triangle_mesh.h>
#include <igl/readOBJ.h>
#include <igl/opengl/glfw/Viewer.h>
#include <Eigen/Core>
#include <iostream>
#include "find_tet.h"
#include "phong_deformation_mesh.h"
#include <cmath>

/**
 * Rotate around z index based on vertex z-position.
 * @param V
 */
Eigen::MatrixX3d deform_mesh(const Eigen::MatrixX3d V)
{
  Eigen::MatrixX3d V_def(V.rows(), 3);

  double rot;

  for (int i = 0; i < V.rows(); i++)
  {
    rot = V(i, 1) * 0.2;

    Eigen::Matrix3d r_mat;

    r_mat << cos(rot), 0, sin(rot),
        0, 1, 0,
        -sin(rot), 0, cos(rot);

    Eigen::Vector3d result = r_mat * V.row(i).transpose();

    V_def.row(i) = result;
  }

  return V_def;
}

int main(int argc, char *argv[])
{
  // args order: full.obj; tet.obj; tet-def.obj

  enum View
  {
    FULL,
    TET,
    TET_DEF,
    FULL_DEF,
  } view = FULL;

  Eigen::MatrixXi ftc, fn; // unused
  Eigen::MatrixXd tc, n;

  // read full (triangle) mesh
  Eigen::MatrixXd V_full;
  Eigen::MatrixXi F_full;
  igl::read_triangle_mesh(argv[1], V_full, F_full);

  // read undeformed tet mesh
  Eigen::MatrixXd V_tet_undef;
  Eigen::MatrixXi F_tet;
  igl::readOBJ(argv[2], V_tet_undef, tc, n, F_tet, ftc, fn);

  V_full.resize(V_full.rows(), 3);
  V_tet_undef.resize(V_tet_undef.rows(), 3);
  F_tet.resize(F_tet.rows(), 4);

  // read deformed tet mesh
  Eigen::MatrixXd V_tet_def;
  V_tet_def.resize(V_tet_undef.rows(), 3);
  // V_tet_def = deform_mesh(V_tet_undef);
  V_tet_def = V_tet_undef;

  Eigen::VectorXi E;
  // find tetrahedron for each point in V_full

  std::cout << std::endl;
  std::cout << "Precomputing phong deformation" << std::endl;
  find_tet(V_full, V_tet_undef, F_tet, E);
  std::cout << "Done" << std::endl;

  V_full.resize(V_full.rows(), 3);
  V_tet_undef.resize(V_tet_undef.rows(), 3);
  F_tet.resize(F_tet.rows(), 4);

  igl::opengl::glfw::Viewer viewer;
  std::cout << R"(
f        Show full mesh (original)
m        Show tet mesh (after deform)
x        Apply deformation to tet mesh
r        Reset deformations
[space]  Show deformed full mesh (via phong deformation)
)";

  auto phong_deform = [&]()
  {
    Eigen::MatrixX3d V_full_def;

    // deform V_full
    phong_deformation_mesh(V_full, V_tet_def, V_tet_undef, F_tet, E, V_full_def);

    viewer.data().clear();
    viewer.data().set_mesh(V_full_def, F_full);
    viewer.data().compute_normals();
  };

  auto update = [&](){
switch (view)
    {
    case TET_DEF:
      viewer.data().set_vertices(V_tet_def);
      viewer.data().compute_normals();
      break;
    case FULL_DEF:
      phong_deform();
      break;
    default:
      break;
    };
  };  

  viewer.callback_key_pressed =
      [&](igl::opengl::glfw::Viewer &, unsigned int key, int mod)
  {
    switch (key)
    {
    case 'f':
      view = FULL;
      viewer.data().clear();
      viewer.data().set_mesh(V_full, F_full);
      break;
    case 'm':
      view = TET_DEF;
      viewer.data().clear();
      viewer.data().set_mesh(V_tet_def, F_tet);
      viewer.data().compute_normals();
      break;
    case ' ':
      view = FULL_DEF;
      viewer.data().clear();

      phong_deform();
      break;

    case 'x':
      V_tet_def = deform_mesh(V_tet_def);
      update();
      break;

    case 'r':
      V_tet_def = V_tet_undef;
      update();
      break;

    default:
      return false;
    }

    return true;
  };

  viewer.data().set_mesh(V_full, F_full);
  viewer.data().show_lines = true;
  viewer.core().is_animating = false;
  viewer.data().face_based = true;

  viewer.launch();

  return EXIT_SUCCESS;
}