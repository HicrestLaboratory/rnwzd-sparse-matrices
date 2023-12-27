#include <iostream>
#include <algorithm>
#include <cassert>
class Matrix
{
public:
};

//////////////////////////////////////////////

#include <list>
#include <vector>
#include <map>
#include <utility>
#include <numeric>
class MatrixCSR : public Matrix
{
private:
    std::vector<size_t> valperrow{};
    std::vector<size_t> cols{};
    std::vector<double> values;

public:
    // TODO public?
    size_t m, n;

    MatrixCSR(size_t m, size_t n) : m{m}, n{n}
    {
    }
    ~MatrixCSR()
    {
    }

};

int main()
{
    MatrixCSR m(2, 3);
    m(0, 0) = 1;
    m(1, 1) = 2;
    std::cout << m << std::endl;
    std::cout << m + m << std::endl;

    return 0;
}
