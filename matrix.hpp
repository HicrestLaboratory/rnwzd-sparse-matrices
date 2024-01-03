#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <iostream>
#include <functional>
#include <vector>
#include <map>

double relu(double x);

class Matrix;
class MatrixD;
class MatrixCOO;

class Matrix
{ // TODO
};

class MatrixD : public Matrix
{
private:
    std::vector<double> data{};

public:
    // TODO public?
    size_t m, n;

    MatrixD(size_t m, size_t n, double value);
    MatrixD(size_t m, size_t n, std::function<double(size_t, size_t)> generator);
    MatrixD(const MatrixCOO &A);
    // TODO rvalue constructor
    ~MatrixD();
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
    friend MatrixD relu(const MatrixD &A);
    friend std::ostream &operator<<(std::ostream &out, const MatrixD &matrix);
};

/////////////////////////////////////

class MatrixCOO : public Matrix
{
private:
    std::map<std::pair<size_t, size_t>, double> data{};

public:
    // TODO public?
    size_t m, n;

    MatrixCOO(size_t m, size_t n);
    MatrixCOO(size_t m, size_t n, std::function<double(size_t, size_t)> generator);
    ~MatrixCOO();

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
    friend MatrixCOO relu(const MatrixCOO &A);
    friend std::ostream &operator<<(std::ostream &out, const MatrixCOO &matrix);
    friend MatrixD::MatrixD(const MatrixCOO &A);
};
///////////////////////////////////
#endif
