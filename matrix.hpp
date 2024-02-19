#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <iostream>
#include <functional>
#include <vector>
#include <map>

constexpr double relu(double x);
constexpr double relu_prime(double x);

class Matrix;
class MatrixD;
class MatrixCOO;

typedef std::vector<std::vector<std::pair<size_t, double>>> AdjList;
typedef std::vector<std::tuple<size_t, size_t, double>> EdgeList;

MatrixCOO mcoo_from_el_file(std::string filename,
                      bool weighted = true,
                      bool directed = true);

class Matrix
{ // TODO
private:
    double value{0};

public:
    // TODO public?
    size_t m, n;
    Matrix();
    Matrix(size_t m, size_t n);
    virtual double &operator()(size_t row, size_t col);
    virtual double operator()(size_t row, size_t col) const;
};

class MatrixD : public Matrix
{
private:
    std::vector<double> data{};

public:
    MatrixD();
    MatrixD(size_t m, size_t n, double value = 0);
    MatrixD(size_t m, size_t n, std::function<double(size_t, size_t)> generator);
    MatrixD(const std::vector<std::vector<double>> &vv);
    MatrixD(const MatrixCOO &A);
    MatrixD(const AdjList &AL);
    ~MatrixD();

    friend class MatrixCOO;
    friend class MatrixCSR;

    double &operator()(size_t row, size_t col) override;
    double operator()(size_t row, size_t col) const override;
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

    friend MatrixD ewprod(const MatrixD &A, const MatrixD &B);
    friend std::ostream &operator<<(std::ostream &out, const MatrixD &matrix);

    friend MatrixD relu(const MatrixD &A);
    friend MatrixD relu_prime(const MatrixD &A);
    double outdegree(size_t v);
    double indegree(size_t v);

    MatrixD t();
};

/////////////////////////////////////

class MatrixCOO : public Matrix
{
private:
    std::map<std::pair<size_t, size_t>, double> data{};

public:
    MatrixCOO(size_t m, size_t n);
    MatrixCOO(size_t m, size_t n, std::function<double(size_t, size_t)> generator);
    MatrixCOO(const std::vector<std::vector<double>> &vv);
    MatrixCOO(const MatrixD &A);
    MatrixCOO(const AdjList &AL, bool directed = true);
    ~MatrixCOO();

    friend class MatrixD;
    friend class MatrixCSR;

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

    friend MatrixD ewprod(const MatrixCOO &A, const MatrixCOO &B);
    friend MatrixD ewprod(const MatrixD &A, const MatrixCOO &B);
    friend MatrixD ewprod(const MatrixCOO &A, const MatrixD &B);
    friend MatrixCOO relu(const MatrixCOO &A);
    friend MatrixCOO relu_prime(const MatrixCOO &A);
    double outdegree(size_t v);
    double indegree(size_t v);

    MatrixCOO t();
};
///////////////////////////////////
#endif
