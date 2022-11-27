#include "vertex_def_gradient.h"
#include <Eigen/Core>
#include <vector>
#include <Eigen/Cholesky>

void vertex_def_gradient(
    const Eigen::MatrixX3d &V_undef,
    const Eigen::MatrixX4i &F,
    const Eigen::MatrixX3d &C,
    const std::vector<Eigen::Matrix3d> &F_cell,
    std::vector<Eigen::Matrix3d> &grad)
{
  // faces adjacent to each vertex
  std::vector<std::vector<int>> face_adj(V_undef.rows());

  for (int i = 0; i < F.rows(); i++)
  {
    // for each vertex in the face
    for (int j = 0; j < F.cols(); j++)
    {
      face_adj[F(i, j)].push_back(i);
    }
  }

  Eigen::MatrixXd A;
  Eigen::Vector3d r_k, r_k_hat, b, lambda;
  int k_count;

  for (int i = 0; i < V_undef.rows(); i++)
  {
    k_count = face_adj[i].size();

    A.setZero();
    b.setZero();
    lambda.setZero();

    A.resize(3, k_count);

    assert(k_count == 3); // should always be 3

    Eigen::LDLT<Eigen::MatrixXd> ldlt;
    Eigen::Vector3d w, r_norms;

    Eigen::Matrix3d R_hats;
    Eigen::Vector3d W; // all weights at this vertex

    // vertex-to-centroid vectors
    for (int k = 0; k < k_count; k++)
    {
      r_k = C.row(face_adj[i][k]) - V_undef.row(i);
      r_k_hat = r_k.normalized();
      R_hats.row(k) = r_k_hat;
      r_norms[k] = r_k.norm();

      A += r_k_hat * r_k_hat.transpose(); 
      b -= r_k_hat;                             // sumation with negative in front
    }

    // regularization
    A += Eigen::MatrixXd::Identity(3, 3); // they use epsilon=1

    // solve for lambda
    ldlt.compute(A);
    lambda = ldlt.solve(b);

    W = R_hats * lambda + Eigen::Vector3d::Ones(); // w_k = 1 + r_k_hat * lambda

    // Sum over j (w_j/r_j)
    double w_over_j_sum = (W.array() / r_norms.array()).sum();

    // centroid-to-vertex weights
    Eigen::Vector3d w_prime = (W.array() / r_norms.array()) / w_over_j_sum;

    // compute gradient
    Eigen::Matrix3d F_i;
    F_i.setZero();

    // sum over each adjacent cell
    for (int k = 0; k < k_count; k++)
    {
      F_i += w_prime[k] * F_cell[face_adj[i][k]];
    }

    grad.push_back(F_i);
  }
}