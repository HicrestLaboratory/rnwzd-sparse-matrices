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
    MatrixD W(){return m_W;}
    MatrixD B(){return m_B;}
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

class Network
{
private:
    std::vector<Layer *> m_layerps{};
    std::function<double(const MatrixD &, const MatrixD &)> m_loss{MSE_loss};
    std::function<MatrixD(const MatrixD &, const MatrixD &)> m_loss_prime{MSE_loss_prime};

public:
    Network()
    {
    }
    Network(std::vector<Layer *> layerps)
        : m_layerps{layerps}
    {
    }
    void add_layerp(Layer *layer)
    {
        m_layerps.push_back(layer);
    }
    MatrixD forward(MatrixD input, const MatrixCOO &A)
    {
        MatrixD output{input};
        for (auto &layerp : m_layerps)
        {
            output = layerp->forward_propagation(output, A);
        }
        return output;
    }
    void fit(std::vector<MatrixD> xs,
             std::vector<MatrixCOO> As,
             std::vector<MatrixD> ys,
             size_t n_epochs, double learning_rate)
    {
        assert(xs.size() == As.size() && xs.size() == ys.size());

        auto n_samples{xs.size()};
        double loss_value{};
        MatrixD output;
        MatrixD error;
        for (size_t epoch{}; epoch < n_epochs; ++epoch)
        {
            loss_value = 0;
            for (size_t j{}; j < n_samples; ++j)
            {
                output = forward(xs[j], As[j]);
                loss_value += m_loss(ys[j], output);

                error = m_loss_prime(ys[j], output);
                for (auto it = m_layerps.rbegin(); it != m_layerps.rend(); ++it)
                {
                    auto layerp = *it;
                    error = layerp->backward_propagation(error, learning_rate);
                }
            }
            loss_value / n_samples;
            std::cout << "Epoch " << epoch+1 << "/" << n_epochs << " loss = " << loss_value << std::endl;
        }
    };
};
int main()
{

    const auto n_samples{100};
    const auto n_v{10};     // number of vertices
    const auto n_h_in{10};  // dimension of input latent vector
    const auto n_h_out{10}; // dimension of output latent vector

    const auto n_epochs{100};
    const double learning_rate{-0.01};

    std::random_device seed;
    std::mt19937 rand_gen(seed());
    std::normal_distribution Normal(0.0, 2.0);
    auto normal_gen = [&](auto i, auto j)
    { return Normal(rand_gen); };

    GNNLayer gnn1{n_h_in, n_h_out, normal_gen};
    ReLULayer relu1{};
    Network net1{{&gnn1, &relu1}};

    std::vector<MatrixD> xs{};
    std::vector<MatrixCOO> As{};
    std::vector<MatrixD> ys{};
    for (size_t i{}; i < n_samples; ++i)
    {
        MatrixD x{n_h_in, n_v, normal_gen};
        MatrixCOO A{n_v, n_v, normal_gen};
        xs.push_back(x);
        As.push_back(A);
        ys.push_back(net1.forward(x, A));
    }

    GNNLayer gnn2{n_h_in, n_h_out, normal_gen};
    ReLULayer relu2{};
    Network net2{{&gnn2, &relu2}};

    net2.fit(xs, As, ys, n_epochs, learning_rate);


    return 0;
}
