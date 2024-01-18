#ifndef GNN_HPP
#define GNN_HPP

#include "matrix.hpp"

double MSE_loss(const MatrixD& A, const MatrixD& B);
MatrixD MSE_loss_prime(const MatrixD& A, const MatrixD& B);


#endif
