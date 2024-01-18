
#include <cmath>
#include <cassert>

#include "matrix.hpp"

double MSE_loss(const MatrixD &Y, const MatrixD &Y_star)
{
    // TODO std parallel
    assert(Y.m == Y_star.m);
    assert(Y.n == Y_star.n);
    double sum{0};
    for (size_t i{}; i < Y.m; ++i)
    {
        for (size_t j{}; j < Y.n; ++j)
        {
            sum += std::pow(Y(i, j) - Y_star(i, j), 2);
        }
    }
    return sum/(Y.m*Y.n);
}
MatrixD MSE_loss_prime(const MatrixD &Y, const MatrixD &Y_star)
{
    return 2*(Y-Y_star)/(Y.m*Y.n);
}
