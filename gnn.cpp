
#include <cmath>
#include <cassert>

#include "matrix.hpp"

double MSELoss(const Matrix& A, const Matrix& B)
{
    // TODO std parallel
    assert(A.m == B.m);
    assert(A.n == B.n);
    double sum{0};
    for (size_t i{}; i < A.m; ++i)
    {
        for (size_t j{}; j < A.n; ++j)
        {   
            sum += std::pow(A(i, j) - B(i, j), 2);
        }
    }
    return std::sqrt(sum);
}
