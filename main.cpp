#include <iostream>
#include <algorithm>
#include <cassert>
#include <list>
#include <vector>
#include <map>
#include <utility>
#include <numeric>
#include <random>

#include "matrix.hpp"

int main()
{
    std::random_device seed;
    std::mt19937 generator(seed());
    std::uniform_real_distribution UniRe(-0.5, 0.5);
    MatrixCOO A{2, 2, [&](auto i, auto j)
                { return UniRe(generator); }};

    std::cout << A << std::endl;
    std::cout << relu(A) << std::endl;
    MatrixD B(A);
    MatrixD M(3, 2, 1);

    MatrixCOO m(2, 3);
    m(0, 0) = 1;
    m(0, 1) = 1;
    m(1, 1) = -2;
    std::cout << M << std::endl;
    std::cout << m << std::endl;
    std::cout << m * M << std::endl;
    relu(M);

    std::cout << MatrixCOO(MatrixD(m)) << std::endl;

    return 0;
}
