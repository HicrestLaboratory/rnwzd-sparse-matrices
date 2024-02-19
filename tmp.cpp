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

class MatrixCSR : public Matrix
{
public: // TODO !!!!
//private 
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
        for (size_t i{0}; i < A.m; ++i)
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

int main()
{
    const double p = 0.3; // prob of nonzero entry
    const size_t N = 1 << 2;

    std::random_device seed;
    std::mt19937 rand_gen(seed());

    std::uniform_real_distribution Unif(0.0, 1.0);
    auto sp_gen = [&](auto i, auto j)
    { return (Unif(rand_gen) < p) ? 1 : 0; };
    auto A{MatrixCSR(MatrixCOO(N, N, sp_gen))};

    std::uniform_int_distribution UnifInt(0, 10);
    auto int_gen = [&](auto i, auto j)
    { return UnifInt(rand_gen); };
    auto v{MatrixD{N, 1, int_gen}};
    ///////////////////////////////////////////////

    
    return 0;
}
