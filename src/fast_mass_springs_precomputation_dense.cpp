#include "fast_mass_springs_precomputation_dense.h"
#include "signed_incidence_matrix_dense.h"
#include <Eigen/Dense>

bool fast_mass_springs_precomputation_dense(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & E,
  const double k,
  const Eigen::VectorXd & m,
  const Eigen::VectorXi & b,
  const double delta_t,
  Eigen::VectorXd & r,
  Eigen::MatrixXd & M,
  Eigen::MatrixXd & A,
  Eigen::MatrixXd & C,
  Eigen::LLT<Eigen::MatrixXd> & prefactorization)
{
  /////////////////////////////////////////////////////////////////////////////
  // Replace with your code
  Eigen::MatrixXd Q = Eigen::MatrixXd::Identity(V.rows(),V.rows());

  //resize
  // M = Eigen::MatrixXd::Identity(m.rows(),m.rows());
  // r.resize(E.rows());
  // C = Eigen::MatrixXd::Identity(b.rows(),V.rows());
  M = Eigen::MatrixXd::Zero(V.rows(), V.rows());
  r.resize(E.rows());
  C = Eigen::MatrixXd::Zero(b.rows(), V.rows());
  


  // find A from sined_incidence_matrix_debse
  signed_incidence_matrix_dense(V.rows(), E, A);


  // Mij = Mi when i = j
  for (int i = 0; i < m.size(); i++){
    for (int j = 0; j < m.size(); j++){
      if (i == j){
        M(i,j ) = m[i];
      }


  }}

  // r
  

  //since E contains 2 list of edge indices, want to find its distance between them in r
  for (int i = 0; i < E.rows(); i++) {

      r[i] = (V.row(E(i, 1)) - V.row(E(i, 0))).norm();
  }


  // selection matrix C
  // put fixed vertices of b into V
  for (int i = 0; i < b.rows(); i++) {
    C(i, b(i)) = 1;
  }

  // pinned vertex penalty
  // equation from handout
  double w=1e10;
  Eigen::MatrixXd term = Eigen::MatrixXd::Identity(V.rows(),V.rows());
  
  // Eigen::MatrixXd P = Eigen::MatrixXd::Identity(b.rows(),V.rows());
  // P = C * V;
  // term = w * P.transpose() * P;
  term = w /2 * (C.transpose() * C); //bcz it is precomputation, so don't need V
  // implement equation from the slide
  Q = (k * A.transpose() * A) + (M / pow(delta_t, 2)) + term;





  /////////////////////////////////////////////////////////////////////////////
  prefactorization.compute(Q);
  return prefactorization.info() != Eigen::NumericalIssue;
}
