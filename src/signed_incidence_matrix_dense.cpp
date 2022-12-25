#include "signed_incidence_matrix_dense.h"

void signed_incidence_matrix_dense(
  const int n,
  const Eigen::MatrixXi & E,
  Eigen::MatrixXd & A)
{
  //////////////////////////////////////////////////////////////////////////////
  // Replace with your code

  
  // equation from the tut slide and under Matrices

  A = Eigen::MatrixXd::Zero(E.rows(),n);

  for (int e = 0; e < E.rows(); e++){
    // if (E(e, 0))
    A(e, E(e, 0)) = 1;  // if i
  	A(e, E(e, 1)) = -1; // if j
  }
  //////////////////////////////////////////////////////////////////////////////
}
