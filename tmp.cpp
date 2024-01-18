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
    std::vector<size_t> valperrow{}; //{0, 2, 4, 6};
    std::vector<size_t> cols{};      //{0, 1, 0, 1, 0, 1};
    std::vector<double> values{};    //{1, 2, 3, 4, 5, 6};

public:
    MatrixCSR(size_t m, size_t n)
        : Matrix(m, n)
    {
    }
    MatrixCSR(const MatrixCOO &A)
        : Matrix(A.m, A.n)
    {
        std::vector<size_t> rows{};
        std::for_each(A.data.begin(), A.data.end(), [&](const auto &e)
                      {
            auto [i,j] {e.first};
            auto A_ij { e.second};
            rows.push_back(i);
            cols.push_back(j);
            values.push_back(A_ij); });
        size_t nnz = values.size();
        valperrow = std::vector<size_t>(m + 1, 0);
        for (size_t i = 0; i < nnz; i++)
            valperrow[rows[i] + 1]++;
        for (size_t i = 0; i < m; i++)
            valperrow[i + 1] += valperrow[i];
    }
    ~MatrixCSR()
    {
    }
    friend MatrixCSR operator*(const MatrixCSR &A, const MatrixD &B)
    {
        assert(A.n == B.m);
        MatrixCOO C{A.m, B.n};
        size_t s{0};
        for (size_t i{0}; i < A.m ; ++i)
        {
            for (size_t h{0}; h < A.valperrow[i + 1] - A.valperrow[i]; ++h)
            {
                auto j = A.cols[s];
                auto A_ij = A.values[s];
                for (int k = 0; k < B.n; ++k)
                    C(i, k) += A_ij * B(j, k);
                s++;
            }
        }

        return MatrixCSR(C);
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
        for (size_t r{0}; r < A.valperrow.size() - 1; ++r)
        {
            for (size_t i{0}; i < A.valperrow[r + 1] - A.valperrow[r]; ++i)
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

// #include <fstream>
// MatrixCOO MatrixCOOFromFile(std::string file_name)
// {
//     std::ifstream file(file_name);

//     if (!file)
//     {
//         std::cerr << "Could not open the file " << file_name << ".";
//         exit(EXIT_FAILURE);
//     }

//     std::vector<std::string> lines;
//     std::string line;
//     std::getline(file, line);
//     // TODO split !!!

//     while (std::getline(file, line))
//         lines.push_back(line);
// }
int main()
{

    MatrixCOO d(2, 2);
    // d(0, 0) = 0;
    d(0, 1) = 1;
    // d(1, 1) = 0;
    d(1, 0) = 1;
    std::cout << (d) << std::endl;
    std::cout << MatrixCSR(d) << std::endl;
    // 2 1
    // 4 3
    // 6 5
    return 0;
}
