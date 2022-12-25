#include "fast_mass_springs_step_sparse.h"
#include <igl/matlab_format.h>

void fast_mass_springs_step_sparse(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & E,
  const double k,
  const Eigen::VectorXi & b,
  const double delta_t,
  const Eigen::MatrixXd & fext,
  const Eigen::VectorXd & r,
  const Eigen::SparseMatrix<double>  & M,
  const Eigen::SparseMatrix<double>  & A,
  const Eigen::SparseMatrix<double>  & C,
  const Eigen::SimplicialLLT<Eigen::SparseMatrix<double> > & prefactorization,
  const Eigen::MatrixXd & Uprev,
  const Eigen::MatrixXd & Ucur,
  Eigen::MatrixXd & Unext)
{
  //////////////////////////////////////////////////////////////////////////////
  // // Replace with your code
  // for(int iter = 0;iter < 50;iter++)
  // {
  //   const Eigen::MatrixXd l = Ucur;
  //   Unext = prefactorization.solve(l);
  //   double w = 1e10;
  //   Eigen::SparseMatrix<double> term;
  //   term  = w /2 *  (C.transpose() * C) * V; //linear term with V

  //   // equation from the tut slide
  //   Eigen::SparseMatrix<double>  y;
  //   y = 1 / pow(delta_t,2) * M * (2 * Ucur - Uprev) + fext + w /2 *  (C.transpose() * C) * V;


  //   std::vector<Eigen::Triplet<double> > ijv;
  //   for(int i = 0;i < 50;i++) //50 iterations 
  //   { 
      

  //     for (int m = 0; m < E.rows(); m++){
  //       Eigen::VectorXd dij = (Unext.row(E(m, 0)) - Unext.row(E(m, 1)) ).normalized() * r[m];
  //       ijv.emplace_back(Eigen::Triplet<double>(m, 0, dij[0]));
  //       ijv.emplace_back(Eigen::Triplet<double>(m, 1, dij[1]));
  //       ijv.emplace_back(Eigen::Triplet<double>(m, 2, dij[2]));
  //     }
  //     Eigen::SparseMatrix<double> d(E.rows(), 3);
  //     d.setFromTriplets(ijv.begin(),ijv.end());
  //     Eigen::SparseMatrix<double> l = k * A.transpose() * d + y;
  //     Unext = prefactorization.solve(l);




  // }
    
  // }

  // I tried to convert dense matrix into sparse above
  // but there are two errors that I can't solve
  // 1. can't get y through C.transpose() * C
  // 2. Unext cannot be solved if l is sparse matrix



    // Replace with your code
  const Eigen::MatrixXd l = Ucur;
  // pinned vertex penalty
  // equation from handout
  double w=1e10;
  Eigen::MatrixXd term = w /2 *  (C.transpose() * C) * V; //linear term with V

  // equation from the tut slide
  Eigen::MatrixXd y = 1 / pow(delta_t,2) * M * (2 * Ucur - Uprev) + fext + term;



  for(int i = 0;i < 50;i++) //50 iterations 
  { 
    Eigen::MatrixXd d = Eigen::MatrixXd::Zero(E.rows(), 3); // R ^ n*3

    for (int m = 0; m < d.rows(); m++){
      // since normalization will have magnitude of 1, r[m], which stores length in rest postion
      // will help get correct d

      d.row(m) = (Unext.row(E(m, 0)) - Unext.row(E(m, 1)) ).normalized() * r[m];

    }
    Eigen::MatrixXd l = k * A.transpose() * d + y;
    Unext = prefactorization.solve(l); // I thinkn it can't solve sparsematrix l into matrix unext?


  }
  
  //////////////////////////////////////////////////////////////////////////////
}
