#include <iostream>
#include <algorithm>
#include <cassert>
#include <list>
#include <vector>
#include <map>
#include <utility>
#include <numeric>
#include <random>
#include <cmath>

#include "matrix.hpp"
#include "gnn.hpp"

class Foo
{
private:
    MatrixD W, B;
};

int main()
{
    std::random_device seed;
    std::mt19937 generator(seed());
    std::uniform_real_distribution UniRe(-0.5, 0.5);
    MatrixCOO A{3, 3, [&](auto i, auto j)
                { return (i < j) ? 1 : 0; }};

    std::cout << A << std::endl;

    MatrixCOO W(2, 2);

    W(0, 0) = 1;
    W(0, 1) = 2;
    W(1, 0) = 3;
    W(1, 1) = 4;

    MatrixD H(2, 3, 1.0);

    std::cout << W << std::endl;
    std::cout << H << std::endl;

    std::cout << W * H * A << std::endl;

    EdgeList edge_list{
        {{1, 1}, {2, 2}},
        {},
        {{1, 3}, {0, 4}}};

    MatrixD C{edge_list};
    std::cout << C << std::endl;

    MatrixCOO D{{{1, 2, 3},
                 {3, 4, 5}}};
    std::cout << D << std::endl;
    MatrixCOO E{{{1, 2, 3},
                 {4, 4, 6}}};
    std::cout << MSELoss(D, E) << std::endl;
    return 0;
}
