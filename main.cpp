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

std::cout<< edge_list.size() <<std::endl;
    MatrixD B{edge_list};
    
    std::cout << B << std::endl;
    return 0;
}
