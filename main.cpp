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
    rot = V(i, 1);
    

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
  // igl::readOBJ(argv[1], V_full, tc, n, F_full, ftc, fn);

  std::cout << "V_full size " << V_full.rows() << " " << V_full.cols() << std::endl;
  std::cout << "F_full size " << F_full.rows() << " " << F_full.cols() << std::endl;

  // read undeformed tet mesh
  Eigen::MatrixXd V_tet;
  Eigen::MatrixXi F_tet;
  igl::readOBJ(argv[2], V_tet, tc, n, F_tet, ftc, fn);

  std::cout << "V_tet size " << V_tet.rows() << " " << V_tet.cols() << std::endl;
  std::cout << "F_tet size " << F_tet.rows() << " " << F_tet.cols() << std::endl;

  V_full.resize(V_full.rows(), 3); // todo fix suzanne mesh
  V_tet.resize(V_tet.rows(), 3);
  F_tet.resize(F_tet.rows(), 4);

  // read deformed tet mesh
  Eigen::MatrixXd V_tet_def;
  // Eigen::MatrixXi F_tet_def;
  // igl::readOBJ(argv[3], V_tet_def, tc, n, F_tet_def, ftc, fn);
  V_tet_def.resize(V_tet.rows(), 3);
  V_tet_def = deform_mesh(V_tet);

  std::cout << "random row def " << V_tet_def.row(24) << std::endl;
  std::cout << "random row undef" << V_tet.row(24) << std::endl;

  igl::opengl::glfw::Viewer viewer;
  std::cout << R"(
f        Show full mesh
u        Show undeformed tet mesh
m        Show deformed tet mesh
[space]  Apply deformation to full mesh
)";

  auto phong_deform = [&]()
  {
    Eigen::VectorXi E;
    Eigen::MatrixX3d V_full_def;

    V_full.resize(V_full.rows(), 3);
    V_tet.resize(V_tet.rows(), 3);
    F_tet.resize(F_tet.rows(), 4);

    std::cout << "Computing tetrahedral belonging for each point" << std::endl;

    // find tetrahedron for each point in V_full
    find_tet(V_full, V_tet, F_tet, E);

    std::cout << "Found all tetrahedra" << std::endl;

    // deform V_full
    phong_deformation_mesh(V_full, V_tet_def, V_tet, F_tet, E, V_full_def);

    viewer.data().clear();
    viewer.data().set_mesh(V_full_def, F_full);
    viewer.data().compute_normals();
    std::cout << "Done" << std::endl;
  };

  viewer.callback_key_pressed =
      [&](igl::opengl::glfw::Viewer &, unsigned int key, int mod)
  {
    switch (key)
    {
    case 'f':
      viewer.data().clear();
      viewer.data().set_mesh(V_full, F_full);
      break;
    case 'u':
      viewer.data().clear();
      viewer.data().set_mesh(V_tet, F_tet);
      viewer.data().compute_normals();
      break;
    case 'm':
      viewer.data().clear();
      viewer.data().set_mesh(V_tet_def, F_tet);
      viewer.data().compute_normals();
      break;
    case ' ':
      viewer.data().clear();

      phong_deform();
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