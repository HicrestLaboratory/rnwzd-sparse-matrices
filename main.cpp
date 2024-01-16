#include <iostream>
#include <algorithm>
#include <cassert>
#include <list>
#include <vector>
#include <map>
#include <random>
#include <cmath>

#include "matrix.hpp"
#include "gnn.hpp"

// Base class
class Layer
{
protected:
    MatrixD m_input;
    MatrixD m_output;

public:
    Layer()
        : m_input{0, 0}, m_output{0, 0} {};
    virtual MatrixD forward_propagation(const MatrixD &input, const MatrixCOO &A) = 0;
    virtual MatrixD backward_propagation(const MatrixD &output_error, double learning_rate) = 0;
};

class GNNLayer : public Layer
{
private:
    size_t m_input_size, m_output_size;
    MatrixD m_W, m_B;
    MatrixCOO m_A{0, 0};

public:
    GNNLayer(size_t input_size, size_t output_size, std::function<double(size_t, size_t)> generator)
        : m_input_size{input_size},
          m_output_size{output_size},
          m_W{output_size, input_size, generator},
          m_B{output_size, input_size, generator}
    {
    }
    MatrixD forward_propagation(const MatrixD &input, const MatrixCOO &A)
    {
        m_input = input;
        m_A = A;
        m_output = m_W * input * A + m_B * input;
        return m_output;
    }
    MatrixD backward_propagation(const MatrixD &output_error, double learning_rate)
    {
        MatrixD W_error{output_error * m_A.t() * m_input.t()};
        MatrixD B_error{output_error * m_input.t()};
        MatrixD input_error{m_W.t() * output_error * m_A.t() + m_B.t() * output_error};

        m_W -= learning_rate * W_error;
        m_B -= learning_rate * B_error;

        return input_error;
    }
};

class ReLULayer : public Layer
{
public:
    ReLULayer()
    {
    }
    MatrixD forward_propagation(const MatrixD &input, const MatrixCOO &A)
    {
        m_input = input;
        m_output = relu(input);
        return m_output;
    }
    MatrixD backward_propagation(const MatrixD &output_error, double learning_rate)
    {
        return ewprod(output_error, relu_prime(m_input));
    }
};

int main()
{
    std::random_device seed;
    std::mt19937 rand_gen(seed());
    std::normal_distribution Normal(0.0, 2.0);
    auto normal_gen = [&](auto i, auto j)
    { return Normal(rand_gen); };
    MatrixD In{2, 4, normal_gen};

    std::cout << In << std::endl;

    MatrixCOO A{4, 4, [](auto i, auto j)
                {
                    return (i < j) ? 1 : 0;
                }};
    GNNLayer ggnl{2, 3, normal_gen};

    MatrixD Out = ggnl.forward_propagation(In, A);

    std::cout << Out << std::endl;

    MatrixD a{{{1, 2, 3},
                 {1, 2, 3}}};
    std::cout << a.m << std::endl;
    std::cout << a.n << std::endl;

    MatrixCOO b{{{0, 1, 0},
                 {1, 0, 1}}};
    std::cout << b.m << std::endl;
    std::cout << b.n << std::endl;
    std::cout << ewprod(a, b) << std::endl;
    return 0;
}
