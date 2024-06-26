
#include <cmath>
#include <cassert>

#include "matrix.hpp"
#include "gnn.hpp"
double MSE_loss(const MatrixD &Y, const MatrixD &Y_star)
{
    // TODO std parallel
    assert(Y.m == Y_star.m);
    assert(Y.n == Y_star.n);
    double sum{0};
    for (size_t i{}; i < Y.m; ++i)
    {
        for (size_t j{}; j < Y.n; ++j)
        {
            sum += std::pow(Y(i, j) - Y_star(i, j), 2);
        }
    }
    return sum / (Y.m * Y.n);
}
MatrixD MSE_loss_prime(const MatrixD &Y, const MatrixD &Y_star)
{
    assert(Y.m == Y_star.m);
    assert(Y.n == Y_star.n);
    return 2 * (Y - Y_star) / (Y.m * Y.n);
}
double BCE_loss(const MatrixD &output, const MatrixD &target)
{
    // TODO std parallel
    assert(output.n == 1 && target.n == 1); // TODO
    assert(output.m == target.m);
    double sum{0}, t, o;
    size_t i = 0;
    for (size_t j{}; j < target.n; ++j)
    {
        t = target(i, j);
        o = output(i, j);
        sum += -t * std::log(o) - (1 - t) * std::log(1 - o);
    }

    return sum;
}
MatrixD BCE_loss_prime(const MatrixD &output, const MatrixD &target)
{
    // TODO std parallel
    assert(output.n == 1 && target.n == 1); // TODO
    assert(output.m == target.m);

    // TODO!!!!!!
    MatrixD C{target.m, target.n, 0};
    for (size_t i = 0; i < target.m; i++)
    {
        for (size_t j = 0; j < target.n; j++)
        {
            C(i, j) = -target(i, j) / output(i, j) + (1 - target(i, j)) / (1 - output(i, j));
        }
    }
    return C;
    // return -ewdiv(target, output)+ewdiv(1-target,1-output);
}

Layer::Layer()
    : m_input{0, 0}, m_output{0, 0} {};

GNNLayer::GNNLayer(size_t input_size, size_t output_size, std::function<double(size_t, size_t)> generator)
    : m_input_size{input_size},
      m_output_size{output_size},
      m_W{input_size, output_size, generator},
      m_B{input_size, output_size, generator}
{
}
MatrixD GNNLayer::forward_propagation(const MatrixD &input, const MatrixCOO &A)
{
    m_input = input;
    m_A = A;
    // m_output = m_W * input * A + m_B * input;
    m_output = A * input * m_W + input * m_B;
    return m_output;
}
MatrixD GNNLayer::backward_propagation(const MatrixD &output_error, double learning_rate)
{
    // MatrixD W_error{output_error * m_A.t() * m_input.t()};
    // MatrixD B_error{output_error * m_input.t()};
    // MatrixD input_error{m_W.t() * output_error * m_A.t() + m_B.t() * output_error};

    MatrixD W_error{m_input.t() * m_A * output_error};
    MatrixD B_error{m_input.t() * output_error};
    MatrixD input_error{m_A * output_error * m_W.t() + output_error * m_B.t()};

    m_W -= learning_rate * W_error;
    m_B -= learning_rate * B_error;

    return input_error;
}
MatrixD GNNLayer::W() { return m_W; }
MatrixD GNNLayer::B() { return m_B; }

ActLayer::ActLayer()
{
}
MatrixD ActLayer::forward_propagation(const MatrixD &input, const MatrixCOO &A)
{
    m_input = input;
    m_output = act(input);
    return m_output;
}
MatrixD ActLayer::backward_propagation(const MatrixD &output_error, double learning_rate)
{
    return ewprod(output_error, act_prime(m_input));
}

Network::Network(std::function<double(const MatrixD &, const MatrixD &)> loss,
                 std::function<MatrixD(const MatrixD &, const MatrixD &)> loss_prime)
    : m_loss{loss}, m_loss_prime{loss_prime}
{
}
Network::Network(std::vector<Layer *> layerps,
                 std::function<double(const MatrixD &, const MatrixD &)> loss,
                 std::function<MatrixD(const MatrixD &, const MatrixD &)> loss_prime)
    : m_layerps{layerps}, m_loss{loss}, m_loss_prime{loss_prime}
{
}
void Network::add_layerp(Layer *layer)
{
    m_layerps.push_back(layer);
}
MatrixD Network::forward(MatrixD input, const MatrixCOO &A)
{
    MatrixD output{input};
    for (auto &layerp : m_layerps)
    {
        output = layerp->forward_propagation(output, A);
    }
    return output;
}
void Network::fit(std::vector<MatrixD> xs,
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
            loss_value += m_loss(output, ys[j]);

            error = m_loss_prime(output, ys[j]);
            for (auto it = m_layerps.rbegin(); it != m_layerps.rend(); ++it)
            {
                auto layerp = *it;
                error = layerp->backward_propagation(error, learning_rate);
            }
        }
        loss_value / n_samples;
        std::cout << "Epoch " << epoch + 1 << "/" << n_epochs << " loss = " << loss_value << std::endl;
    }
};
