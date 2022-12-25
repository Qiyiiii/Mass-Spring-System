#include "fast_mass_springs_precomputation_sparse.h"
#include "signed_incidence_matrix_sparse.h"
#include <vector>

bool fast_mass_springs_precomputation_sparse(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & E,
  const double k,
  const Eigen::VectorXd & m,
  const Eigen::VectorXi & b,
  const double delta_t,
  Eigen::VectorXd & r,
  Eigen::SparseMatrix<double>  & M,
  Eigen::SparseMatrix<double>  & A,
  Eigen::SparseMatrix<double>  & C,
  Eigen::SimplicialLLT<Eigen::SparseMatrix<double> > & prefactorization)
{
  /////////////////////////////////////////////////////////////////////////////
  // Replace with your code
  // Replace with your code
  std::vector<Eigen::Triplet<double> > ijv;
  int n = V.rows();
  Eigen::SparseMatrix<double> Q(n,n);





  //resize
  // M = Eigen::MatrixXd::Identity(m.rows(),m.rows());
  // r.resize(E.rows());
  // C = Eigen::MatrixXd::Identity(b.rows(),V.rows());
  M.resize(V.rows(), V.rows());
  r.resize(E.rows());
  C.resize(b.rows(), V.rows());
  


  // find A from sined_incidence_matrix_debse
  signed_incidence_matrix_sparse(V.rows(), E, A);


  // Mij = Mi when i = j
  for (int i = 0; i < m.size(); i++){
    for (int j = 0; j < m.size(); j++){
      if (i == j){
        ijv.emplace_back(Eigen::Triplet<double>(i, j, m(i)));
      }


  }}
  M.setFromTriplets(ijv.begin(),ijv.end());

  // r
  

  //since E contains 2 list of edge indices, want to find its distance between them in r
  for (int i = 0; i < E.rows(); i++) {

      r[i] = (V.row(E(i, 1)) - V.row(E(i, 0))).norm();
  }

  std::vector<Eigen::Triplet<double> > ijv1;
  // selection matrix C
  // put fixed vertices of b into V
  for (int i = 0; i < b.rows(); i++) {
    ijv1.emplace_back(Eigen::Triplet<double>(i, b(i), 1));
  }
  C.setFromTriplets(ijv1.begin(),ijv1.end());

  // pinned vertex penalty
  // equation from handout
  double w=1e10;
  Eigen::SparseMatrix<double> term(n,n);
  
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
