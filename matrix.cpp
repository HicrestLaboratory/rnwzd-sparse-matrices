#include <iostream>
#include <algorithm>
#include <cassert>
#include <vector>

#include "matrix.hpp"

double relu(double x)
{
    return (x > 0) ? x : 0;
}

/////////////////////////////////

MatrixD::MatrixD(size_t m, size_t n, double value = 0) : m{m}, n{n}
{
    data.resize(m * n);
    std::fill(data.begin(), data.end(), value);
}
MatrixD::MatrixD(size_t m, size_t n, std::function<double(size_t, size_t)> generator) : m{m}, n{n}
{ // TODO std generate parallel
    data.resize(m * n);
    for (size_t i{}; i < m; ++i)
    {
        for (size_t j{}; j < n; ++j)
        {
            data[i * n + j] = generator(i, j);
        }
    }
}
MatrixD::MatrixD(const MatrixCOO &A) : m{A.m}, n{A.n}
{
    data.resize(m * n);
    std::for_each(A.data.begin(), A.data.end(),
                  [&](const auto &e)
                  {
                    auto [i,j] {e.first};
                    auto A_ij { e.second};
                    data[i * n + j] = A_ij; });
}
// TODO rvalue constructor
MatrixD::~MatrixD()
{
}
double &MatrixD::operator()(size_t row, size_t col)
{
    assert(row >= 0 && row < m);
    assert(col >= 0 && col < n);

    return data[row * n + col];
}

double MatrixD::operator()(size_t row, size_t col) const
{
    assert(row >= 0 && row < m);
    assert(col >= 0 && col < n);

    return data[row * n + col];
}
double &MatrixD::operator[](size_t index)
{
    assert(index >= 0 && index < m * n);
    return data[index];
}
double MatrixD::operator[](size_t index) const
{
    assert(index >= 0 && index < m * n);
    return data[index];
}
MatrixD MatrixD::operator+() const
{ //?? is this right? gets a copy?
    return *this;
}
MatrixD MatrixD::operator-() const
{
    MatrixD B{m, n, 0};
    std::transform(data.begin(), data.end(), B.data.begin(), [](const double e)
                   { return -e; });
    return B;
}
MatrixD operator+(const MatrixD &A, const MatrixD &B)
{
    assert(A.m == B.m && A.n == B.n);
    MatrixD C{A.m, A.n, 0};
    std::transform(A.data.begin(), A.data.end(), B.data.begin(), C.data.begin(), std::plus<double>());
    return C;
}

MatrixD operator-(const MatrixD &A, const MatrixD &B)
{
    return A + (-B);
}
MatrixD operator+=(MatrixD &A, const MatrixD &B)
{
    assert(A.m == B.m && A.n == B.n);
    MatrixD C{MatrixD(A.m, A.n, 0)};
    std::transform(A.data.begin(), A.data.end(), B.data.begin(), A.data.begin(), std::plus<double>());
    return A;
}
MatrixD operator-=(MatrixD &A, const MatrixD &B)
{
    return A += (-B);
}

MatrixD operator*(const MatrixD &A, double value)
{
    MatrixD B{MatrixD(A.m, A.n, 0)};
    std::transform(A.data.begin(), A.data.end(), B.data.begin(), [value](const double e)
                   { return e * value; });
    return B;
}
MatrixD operator*(double value, const MatrixD &A)
{
    return A * value;
}

MatrixD operator/(const MatrixD &A, double value)
{
    return A * (1 / value);
}

MatrixD operator*(const MatrixD &A, const MatrixD &B)
{ // TODO use std algorithms
    assert(A.n == B.m);

    MatrixD C{A.m, B.n, 0};
    for (int i = 0; i < A.m; ++i)
    {
        for (int k = 0; k < B.n; ++k)
        {
            for (int j = 0; j < A.n; ++j)
            {
                C(i, k) += A(i, j) * B(j, k);
            }
        }
    }
    return C;
}
MatrixD relu(const MatrixD &A)
{
    MatrixD C{A.m, A.n};
    std::transform(
        A.data.begin(), A.data.end(), C.data.begin(), [](const auto &e)
        { return relu(e); });
    return C;
}
std::ostream &operator<<(std::ostream &out, const MatrixD &matrix)
{
    for (int x = 0; x < matrix.m; x++)
    {
        for (int y = 0; y < matrix.n; y++)
        {
            out << matrix(x, y) << " ";
        }
        out << std::endl;
    };
    return out;
}

/////////////////////////////////////////

MatrixCOO::MatrixCOO(size_t m, size_t n) : m{m}, n{n}
{
}
// MatrixCOO(const MatrixCOO &A) : m{A.m}, n{A.n}
// {
//     std::copy(A.data.begin(),A.data.end(),data.begin() );
// }
MatrixCOO::MatrixCOO(size_t m, size_t n, std::function<double(size_t, size_t)> generator) : m{m}, n{n}
{ // TODO std generate parallel
    double A_ij;
    for (size_t i{}; i < m; ++i)
    {
        for (size_t j{}; j < n; ++j)
        {
            A_ij = generator(i, j);
            if (A_ij != 0.0)
            {
                data[{i, j}] = A_ij;
            }
        }
    }
}
// TODO rvalue constructor
MatrixCOO::~MatrixCOO()
{
}
// MatrixCOO operator=(const MatrixCOO &A)
// {
//     std::copy(A.data.begin(),A.data.end(),data.begin() );
//     return *this;
// }
double &MatrixCOO::operator()(size_t row, size_t col)
{
    assert(row >= 0 && row < m);
    assert(col >= 0 && col < n);

    return data[{row, col}];
}
double MatrixCOO::operator()(size_t row, size_t col) const
{
    assert(row >= 0 && row < m);
    assert(col >= 0 && col < n);
    double value = 0;
    try
    {
        value = data.at({row, col});
    }
    catch (std::out_of_range oor)
    {
    }
    return value;
}
double &MatrixCOO::operator[](size_t index)
{ // returns with liear index row major
    assert(index >= 0 && index < m * n);
    size_t row = index / n;
    size_t col = index % n;
    return data[{row, col}];
}
double MatrixCOO::operator[](size_t index) const
{ // returns with liear index row major
    assert(index >= 0 && index < m * n);
    size_t row = index / n;
    size_t col = index % n;
    double value = 0;
    try
    {
        value = data.at({row, col});
    }
    catch (std::out_of_range oor)
    {
    }
    return value;
}
MatrixCOO MatrixCOO::operator+() const
{ //?? is this right? gets a copy?
    return *this;
}
MatrixCOO MatrixCOO::operator-() const
{
    MatrixCOO B{MatrixCOO(m, n)};

    std::transform(data.begin(), data.end(), std::inserter(B.data, B.data.begin()), [](auto &e)
                   { return std::make_pair(e.first, -e.second); });
    return B;
}

MatrixCOO operator+(const MatrixCOO &A, const MatrixCOO &B)
{
    assert(A.m == B.m && A.n == B.n);
    MatrixCOO C(A);
    std::for_each(B.data.begin(), B.data.end(), [&C](auto &e)
                  { C(e.first.first, e.first.second) += e.second; });
    return C;
}

MatrixD operator+(const MatrixD &A, const MatrixCOO &B)
{

    assert(A.m == B.m && A.n == B.n);
    MatrixD C(A);
    std::for_each(B.data.begin(), B.data.end(), [&C](auto &e)
                  { C(e.first.first, e.first.second) += e.second; });

    return C;
}
MatrixD operator+(const MatrixCOO &A, const MatrixD &B)
{

    return B + A;
}
MatrixCOO operator-(const MatrixCOO &A, const MatrixCOO &B)
{
    return A + (-B);
}
MatrixCOO operator+=(MatrixCOO &A, const MatrixCOO &B)
{
    A = A + B;
    return A;
}
MatrixCOO operator-=(MatrixCOO &A, const MatrixCOO &B)
{
    return A += (-B);
}

MatrixCOO operator*(const MatrixCOO &A, double value)
{
    MatrixCOO B{A.m, A.n};
    std::for_each(A.data.begin(), A.data.end(), [value, &B](const auto e)
                  { B(e.first.first, e.first.second) = e.second * value; });
    return B;
}
MatrixCOO operator*(double value, const MatrixCOO &A)
{
    return A * value;
}

MatrixCOO operator/(const MatrixCOO &A, double value)
{
    return A * (1 / value);
}
MatrixCOO operator*(const MatrixCOO &A, const MatrixD &B)
{
    assert(A.n == B.m);
    MatrixCOO C{A.m, B.n};
    std::for_each(A.data.begin(), A.data.end(), [&B, &C](const auto &e)
                  {
            auto [i,j] {e.first};
            auto A_ij { e.second};
            for (size_t k {0}; k < B.n; ++k)
                 C(i,k) += A_ij * B(j,k); });
    return C;
}
MatrixCOO operator*(const MatrixCOO &A, const MatrixCOO &B)
{
    assert(A.n == B.m);
    MatrixCOO C{A.m, B.n};
    std::for_each(A.data.begin(), A.data.end(), [&B, &C](const auto &e)
                  {
            auto [i,j] {e.first};
            auto A_ij { e.second};
            for (int k = 0; k < B.n; ++k)
                 C(i,k) += A_ij * B(j,k); });
    return C;
}
MatrixCOO relu(const MatrixCOO &A)
{
    MatrixCOO C{A.m, A.n};
    std::for_each(A.data.begin(), A.data.end(), [&C](const auto &e)
                  {
            auto [i,j] {e.first};
            auto A_ij { e.second};
            auto C_ij = relu(A_ij); 
            if ( C_ij != 0)
                {C.data[{i,j}] = C_ij;}
            else
                {C.data.erase({i,j});} });
    return C;
}
std::ostream &operator<<(std::ostream &out, const MatrixCOO &matrix)
{
    for (size_t x = 0; x < matrix.m; x++)
    {
        for (size_t y = 0; y < matrix.n; y++)
        {
            out << matrix(x, y) << " ";
        }
        out << std::endl;
    };

    // for (const auto &[k, v] : matrix.data)
    //     std::cout << "m[" << k.first << "," << k.second << "] = " << v << std::endl;

    return out;
}
