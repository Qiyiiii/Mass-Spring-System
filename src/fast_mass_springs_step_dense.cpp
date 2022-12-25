#include "fast_mass_springs_step_dense.h"
#include <igl/matlab_format.h>

void fast_mass_springs_step_dense(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & E,
  const double k,
  const Eigen::VectorXi & b,
  const double delta_t,
  const Eigen::MatrixXd & fext,
  const Eigen::VectorXd & r,
  const Eigen::MatrixXd & M,
  const Eigen::MatrixXd & A,
  const Eigen::MatrixXd & C,
  const Eigen::LLT<Eigen::MatrixXd> & prefactorization,
  const Eigen::MatrixXd & Uprev,
  const Eigen::MatrixXd & Ucur,
  Eigen::MatrixXd & Unext)
{
  //////////////////////////////////////////////////////////////////////////////
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
    Unext = prefactorization.solve(l);



  }
  
  //////////////////////////////////////////////////////////////////////////////
}
