#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <iostream>
#include <functional>
#include <vector>
#include <map>

constexpr double relu(double x);

class Matrix;
class MatrixD;
class MatrixCOO;

typedef std::vector<std::vector<std::pair<size_t, double>>> EdgeList;

class Matrix
{ // TODO
};

class MatrixD : public Matrix
{
public:
    // TODO public?
    size_t m, n;

private:
    std::vector<double> data{};

public:
    MatrixD(size_t m, size_t n, double value = 0);
    MatrixD(size_t m, size_t n, std::function<double(size_t, size_t)> generator);
    MatrixD(const std::vector<std::vector<double>> &vv);
    MatrixD(const MatrixCOO &A);
    MatrixD(const EdgeList &EL);
    ~MatrixD();

    friend class MatrixCOO;

    double &operator()(size_t row, size_t col);
    double operator()(size_t row, size_t col) const;
    double &operator[](size_t index);
    double operator[](size_t index) const;
    MatrixD operator+() const;
    MatrixD operator-() const;
    friend MatrixD operator+(const MatrixD &A, const MatrixD &B);
    friend MatrixD operator-(const MatrixD &A, const MatrixD &B);
    friend MatrixD operator+=(MatrixD &A, const MatrixD &B);
    friend MatrixD operator-=(MatrixD &A, const MatrixD &B);
    friend MatrixD operator*(const MatrixD &A, double value);
    friend MatrixD operator*(double value, const MatrixD &A);
    friend MatrixD operator/(const MatrixD &A, double value);
    friend MatrixD operator*(const MatrixD &A, const MatrixD &B);

    friend std::ostream &operator<<(std::ostream &out, const MatrixD &matrix);

    friend MatrixD relu(const MatrixD &A);
    double outdegree(size_t v);
    double indegree(size_t v);

    MatrixD t();
};

/////////////////////////////////////

class MatrixCOO : public Matrix
{
public:
    // TODO public?
    size_t m, n;

private:
    std::map<std::pair<size_t, size_t>, double> data{};

public:
    MatrixCOO(size_t m, size_t n);
    MatrixCOO(size_t m, size_t n, std::function<double(size_t, size_t)> generator);
    MatrixCOO(const std::vector<std::vector<double>> &vv);
    MatrixCOO(const MatrixD &A);
    MatrixCOO(const EdgeList &EL);
    ~MatrixCOO();

    friend class MatrixD;

    double &operator()(size_t row, size_t col);
    double operator()(size_t row, size_t col) const;
    double &operator[](size_t index);
    double operator[](size_t index) const;
    MatrixCOO operator+() const;
    MatrixCOO operator-() const;
    friend MatrixCOO operator+(const MatrixCOO &A, const MatrixCOO &B);
    friend MatrixD operator+(const MatrixD &A, const MatrixCOO &B);
    friend MatrixD operator+(const MatrixCOO &A, const MatrixD &B);
    friend MatrixCOO operator-(const MatrixCOO &A, const MatrixCOO &B);
    friend MatrixCOO operator+=(MatrixCOO &A, const MatrixCOO &B);
    friend MatrixCOO operator-=(MatrixCOO &A, const MatrixCOO &B);
    friend MatrixCOO operator*(const MatrixCOO &A, double value);
    friend MatrixCOO operator*(double value, const MatrixCOO &A);
    friend MatrixCOO operator/(const MatrixCOO &A, double value);
    friend MatrixCOO operator*(const MatrixCOO &A, const MatrixD &B);
    friend MatrixCOO operator*(const MatrixCOO &A, const MatrixCOO &B);
    friend MatrixCOO operator*(const MatrixD &A, const MatrixCOO &B);
    friend std::ostream &operator<<(std::ostream &out, const MatrixCOO &matrix);

    friend MatrixCOO relu(const MatrixCOO &A);
    double outdegree(size_t v);
    double indegree(size_t v);

    MatrixCOO t();
};
///////////////////////////////////
#endif
