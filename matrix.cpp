#include <iostream>
#include <algorithm>
#include <cassert>
#include <vector>

#include "matrix.hpp"

constexpr double relu(double x)
{
    return (x > 0) ? x : 0;
}
constexpr double relu_prime(double x)
{
    return (x > 0) ? 1 : 0;
}

//////////////////////////////

Matrix::Matrix()
    : m{0}, n{0}
{
}
Matrix::Matrix(size_t m, size_t n)
    : m{m}, n{n}
{
}
double &Matrix::operator()(size_t row, size_t col)
{
    return value;
}
double Matrix::operator()(size_t row, size_t col) const
{
    return value;
}
/////////////////////////////////
MatrixD::MatrixD()
    : Matrix(0,0), data{}
{
}
MatrixD::MatrixD(size_t m, size_t n, double value)
    : Matrix(m, n), data{std::vector<double>(m * n, value)}
{
}
MatrixD::MatrixD(size_t m, size_t n, std::function<double(size_t, size_t)> generator)
    : Matrix(m, n), data{std::vector<double>(m * n)}
{ // TODO std generate parallel
    for (size_t i{}; i < m; ++i)
    {
        for (size_t j{}; j < n; ++j)
        {
            data[i * n + j] = generator(i, j);
        }
    }
}
MatrixD::MatrixD(const std::vector<std::vector<double>> &vv)
    : Matrix(vv.size(), vv[0].size()), data{std::vector<double>(m * n)}
{ // TODO std generate parallel
    for (size_t i{}; i < m; ++i)
    {
        assert(vv[i].size() == n);
        for (size_t j{}; j < n; ++j)
        {
            data[i * n + j] = vv[i][j];
        }
    }
}
MatrixD::MatrixD(const MatrixCOO &A)
    : Matrix(A.m, A.n), data{std::vector<double>(A.m * A.n)}
{
    std::for_each(A.data.begin(), A.data.end(),
                  [&](const auto &e)
                  {
                    auto [i,j] {e.first};
                    auto A_ij { e.second};
                    data[i * n + j] = A_ij; });
}
MatrixD::MatrixD(const EdgeList &EL)
    : Matrix(EL.size(), EL.size()), data{std::vector<double>(EL.size() * EL.size())}
{
    for (size_t v{}; v < n; ++v)
    {
        for (auto [u, w] : EL[v])
        {
            data[v * n + u] = w;
        }
    }
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
MatrixD ewprod(const MatrixD &A, const MatrixD &B)
{
    assert(A.m == B.m);
    assert(A.n == B.n);
    MatrixD C{A.m, B.n, 0};
    for (size_t i = 0; i < A.m; i++)
    {
        for (size_t j = 0; j < A.n; j++)
        {
            C(i, j) = A(i, j) * B(i, j);
        }
    }
    return C;
}
std::ostream &operator<<(std::ostream &out, const MatrixD &matrix)
{
    for (size_t i = 0; i < matrix.m; i++)
    {
        for (size_t j = 0; j < matrix.n; j++)
        {
            out << matrix(i, j) << " ";
        }
        out << std::endl;
    };
    return out;
}
MatrixD relu(const MatrixD &A)
{
    MatrixD C{A.m, A.n};
    std::transform(
        A.data.begin(), A.data.end(), C.data.begin(), [](const auto &e)
        { return relu(e); });
    return C;
}
MatrixD relu_prime(const MatrixD &A)
{
    MatrixD C{A.m, A.n};
    std::transform(
        A.data.begin(), A.data.end(), C.data.begin(), [](const auto &e)
        { return relu_prime(e); });
    return C;
}
double MatrixD::outdegree(size_t v)
{
    // TODO unweighted ??
    // TODO std algorithm
    double sum{0};
    auto i{v};
    for (size_t j{}; j < n; ++j)
    {
        sum += data[i * n + j];
    }
    return sum;
}
double MatrixD::indegree(size_t v)
{
    // TODO unweighted ??
    // TODO std algorithm
    double sum{0};
    auto j{v};
    for (size_t i{}; i < m; ++i)
    {
        sum += data[i * n + j];
    }
    return sum;
}
MatrixD MatrixD::t()
{
    MatrixD T{
        n, m, [&](size_t j, size_t i)
        {
            return (*this)(i, j);
        }};
    return T;
}
/////////////////////////////////////////

MatrixCOO::MatrixCOO(size_t m, size_t n)
    : Matrix(m, n)
{
}
MatrixCOO::MatrixCOO(size_t m, size_t n, std::function<double(size_t, size_t)> generator)
    : Matrix(m, n)
{ // TODO std algorithm generate parallel
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
MatrixCOO::MatrixCOO(const std::vector<std::vector<double>> &vv)
    : Matrix(vv.size(), vv[0].size())
{ // TODO std generate parallel
    double A_ij;
    for (size_t i{}; i < m; ++i)
    {
        assert(vv[i].size() == n);
        for (size_t j{}; j < n; ++j)
        {
            A_ij = vv[i][j];
            if (A_ij != 0.0)
            {
                data[{i, j}] = A_ij;
            }
        }
    }
}
MatrixCOO::MatrixCOO(const MatrixD &A)
    : Matrix(A.m, A.n)
{
    // TODO std algorithm parallel
    double A_ij;
    for (size_t i{}; i < m; ++i)
    {
        for (size_t j{}; j < n; ++j)
        {
            A_ij = A(i, j);
            if (A_ij != 0.0)
            {
                data[{i, j}] = A_ij;
            }
        }
    }
}
MatrixCOO::MatrixCOO(const EdgeList &EL)
    : Matrix(EL.size(), EL.size())
{
    for (size_t v{}; v < n; ++v)
    {
        for (auto [u, w] : EL[v])
        {
            data[{v, u}] = w;
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
    catch (const std::out_of_range &e)
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
    catch (const std::out_of_range &e)
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
MatrixCOO operator*(const MatrixD &A, const MatrixCOO &B)
{
    assert(A.n == B.m);
    MatrixCOO C{A.m, B.n};
    std::for_each(B.data.begin(), B.data.end(), [&](const auto &e)
                  {
            auto [j,k] {e.first};
            auto B_jk { e.second};
            for (int i = 0; i < A.m; ++i)
                 C(i,k) += A(i,j) * B_jk; });
    return C;
}
MatrixD ewprod(const MatrixCOO &A, const MatrixCOO &B)
{   
    
    assert(A.m == B.m);
    assert(A.n == B.n);
    MatrixCOO C{A.m, A.n};
    std::for_each(A.data.begin(), A.data.end(), [&B, &C](const auto &e)
                  {
            auto [i,j] {e.first};
            auto A_ij { e.second};
            C(i,j) = A_ij * B(i,j); });
    return C;
}
MatrixD ewprod(const MatrixD &A, const MatrixCOO &B)
{   
    
    assert(A.m == B.m);
    assert(A.n == B.n);
    MatrixCOO C{A.m, A.n};
    std::for_each(B.data.begin(), B.data.end(), [&A, &C](const auto &e)
                  {
            auto [i,j] {e.first};
            auto B_ij { e.second};
            C(i,j) = B_ij * A(i,j); });
    return C;
}
MatrixD ewprod(const MatrixCOO &A, const MatrixD &B)
{   
    
    assert(A.m == B.m);
    assert(A.n == B.n);
    MatrixCOO C{A.m, A.n};
    std::for_each(A.data.begin(), A.data.end(), [&B, &C](const auto &e)
                  {
            auto [i,j] {e.first};
            auto A_ij { e.second};
            C(i,j) = A_ij * B(i,j); });
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
MatrixCOO relu_prime(const MatrixCOO &A)
{
    MatrixCOO C{A.m, A.n};
    std::for_each(A.data.begin(), A.data.end(), [&C](const auto &e)
                  {
            auto [i,j] {e.first};
            auto A_ij { e.second};
            auto C_ij = relu_prime(A_ij); 
            if ( C_ij != 0)
                {C.data[{i,j}] = C_ij;}
            else
                {C.data.erase({i,j});} });
    return C;
}
double MatrixCOO::outdegree(size_t v)
{
    // TODO unweighted ??
    // TODO std algorithm
    double sum{0};
    std::for_each(data.begin(), data.end(), [&](const auto &e)
                  {
            auto [i,j] {e.first};
            auto A_ij { e.second};
            if (i == v) sum += A_ij; });
    return sum;
}
double MatrixCOO::indegree(size_t v)
{
    // TODO unweighted ??
    // TODO std algorithm
    double sum{0};
    std::for_each(data.begin(), data.end(), [&](const auto &e)
                  {
            auto [i,j] {e.first};
            auto A_ij { e.second};
            if (j == v) sum += A_ij; });
    return sum;
}
MatrixCOO MatrixCOO::t()
{
    MatrixCOO T{n, m};
    std::for_each(data.begin(), data.end(), [&](const auto &e)
                  {
            auto [i,j] {e.first};
            auto A_ij { e.second};
            T(j,i) = A_ij; });
    return T;
}
