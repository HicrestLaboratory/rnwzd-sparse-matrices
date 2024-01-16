#include <iostream>
#include <algorithm>
#include <cassert>
#include <list>
#include <vector>
#include <map>
#include <utility>
#include <numeric>

#include "matrix.hpp"

class MatrixCSR : public Matrix
{
private:
    std::vector<size_t> valperrow{2, 2, 2};
    std::vector<size_t> cols{0, 1, 0, 1, 0, 1};
    std::vector<double> values{1, 2, 3, 4, 5, 6};

public:
    // TODO public?
    size_t m, n;

    MatrixCSR(size_t m, size_t n) : m{m}, n{n}
    {
    }
    ~MatrixCSR()
    {
    }
    friend MatrixCOO operator*(const MatrixCSR &A, const MatrixD &B)
    {
        assert(A.n == B.m);
        MatrixCOO C{A.m, B.n};
        size_t s{0};
        for (size_t i{0}; i < A.valperrow.size(); ++i)
        {
            for (size_t h{0}; h < A.valperrow[i]; ++h)
            {
                auto j = A.cols[s];
                auto A_ij = A.values[s];
                for (int k = 0; k < B.n; ++k)
                    C(i, k) += A_ij * B(j, k);
                s++;
            }
        }

        return C;
    }
    // MatrixCSR operator*(const MatrixD &B)
    // {
    //     assert(n == B.m);
    //     MatrixCSR C{m, B.n};
    //     size_t s{0};
    //     double C_ik{};
    //     for (size_t i{0}; i < valperrow.size(); ++i)
    //     {
    //         C.valperrow.push_back(0);
    //         for (size_t h{0}; h < valperrow[i]; ++h)
    //         {
    //             C_ik = 0;
    //             auto j = cols[s];
    //             auto A_ij = values[s];
    //             for (int k = 0; k < B.n; ++k)
    //                 C_ik += A_ij * B(j, k);

    //             C.valperrow[i] += 1;
    //             C.cols.push_back(k);
    //             C.values.push_back(C_ik);
    //             s++;
    //         }
    //     }

    //     return C;
    // }
    friend std::ostream &operator<<(std::ostream &out, const MatrixCSR &A)
    {
        size_t s{0};
        for (size_t r{0}; r < A.valperrow.size(); ++r)
        {
            for (size_t i{0}; i < A.valperrow[r]; ++i)
            {
                auto c = A.cols[s];
                auto value = A.values[s];
                out << "(" << r << ", " << c << "): " << value << std::endl;
                s++;
            }
        }

        return out;
    }
};

#include <fstream>
MatrixCOO MatrixCOOFromFile(std::string file_name)
{
    std::ifstream file(file_name);

    if (!file)
    {
        std::cerr << "Could not open the file " << file_name << ".";
        exit(EXIT_FAILURE);
    }

    std::vector<std::string> lines;
    std::string line;
    std::getline(file, line);
    // TODO split !!!

    while (std::getline(file, line))
        lines.push_back(line);
}
int main()
{
    MatrixCSR m(3, 2);
    MatrixD d(2, 2);
    d(0, 0) = 0;
    d(0, 1) = 1;
    d(1, 1) = 0;
    d(1, 0) = 1;
    std::cout << m << std::endl;
    std::cout << d << std::endl;
    std::cout << m * d << std::endl;

    return 0;
}
