#include "vertex_def_gradient.h"
#include <Eigen/Core>
#include <vector>
#include <Eigen/Cholesky>
#include <iostream>

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

  grad.resize(V_undef.rows());

  // constant size
  Eigen::Matrix3d A;
  Eigen::Vector3d r_k, r_k_hat, b;
  Eigen::Vector3d lambda;
  int k_count;

  for (int i = 0; i < V_undef.rows(); i++)
  {
    k_count = face_adj[i].size();

    A.setZero();
    b.setZero();

    Eigen::LDLT<Eigen::MatrixXd> ldlt;
    Eigen::VectorXd w(k_count), r_norms(k_count);

    Eigen::MatrixXd R_hats(k_count, 3);
    Eigen::VectorXd W(k_count); // all weights at this vertex

    // vertex-to-centroid vectors
    for (int k = 0; k < k_count; k++)
    {
      r_k = C.row(face_adj[i][k]) - V_undef.row(i);
      r_k_hat = r_k.normalized();
      R_hats.row(k) = r_k_hat;
      r_norms[k] = r_k.norm();

      A += r_k_hat * r_k_hat.transpose();
      b -= r_k_hat; // sumation with negative in front
    }

    // regularization
    A += Eigen::MatrixXd::Identity(3, 3); // they use epsilon=1

    // solve for lambda
    ldlt.compute(A);
    lambda = ldlt.solve(b);

    W = R_hats * lambda; // w_k = 1 + r_k_hat * lambda
    W += Eigen::VectorXd::Ones(k_count);

    // Sum over j (w_j/r_j)
    double w_over_j_sum = (W.array() / r_norms.array()).sum();

    // centroid-to-vertex weights
    Eigen::VectorXd w_prime = W.array() / r_norms.array();
    w_prime /= w_over_j_sum;

    // compute gradient
    Eigen::Matrix3d F_i;
    F_i.setZero();

    // sum over each adjacent cell
    for (int k = 0; k < k_count; k++)
    {
      F_i += w_prime[k] * F_cell[face_adj[i][k]];
    }

    grad[i] = F_i;
  }
}