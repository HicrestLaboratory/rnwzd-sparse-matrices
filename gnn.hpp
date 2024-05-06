#ifndef GNN_HPP
#define GNN_HPP

#include "matrix.hpp"

double MSE_loss(const MatrixD &A, const MatrixD &B);
MatrixD MSE_loss_prime(const MatrixD &A, const MatrixD &B);

double BCE_loss(const MatrixD &output, const MatrixD &target);
MatrixD BCE_loss_prime(const MatrixD &output, const MatrixD &target);

class Layer
{
protected:
    MatrixD m_input;
    MatrixD m_output;

public:
    Layer();
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
    GNNLayer(size_t input_size, size_t output_size, std::function<double(size_t, size_t)> generator);
    MatrixD forward_propagation(const MatrixD &input, const MatrixCOO &A);
    MatrixD backward_propagation(const MatrixD &output_error, double learning_rate);
    MatrixD W();
    MatrixD B();
};

class ActLayer : public Layer
{
public:
    ActLayer();
    MatrixD forward_propagation(const MatrixD &input, const MatrixCOO &A);
    MatrixD backward_propagation(const MatrixD &output_error, double learning_rate);
};

class Network
{
private:
    std::vector<Layer *> m_layerps{};
    std::function<double(const MatrixD &, const MatrixD &)> m_loss{BCE_loss};
    std::function<MatrixD(const MatrixD &, const MatrixD &)> m_loss_prime{BCE_loss_prime};

public:
    Network(
        std::function<double(const MatrixD &, const MatrixD &)> loss,
        std::function<MatrixD(const MatrixD &, const MatrixD &)> loss_prime);
    Network(std::vector<Layer *> layerps,
            std::function<double(const MatrixD &, const MatrixD &)> loss,
            std::function<MatrixD(const MatrixD &, const MatrixD &)> loss_prime);
    void add_layerp(Layer *layer);
    MatrixD forward(MatrixD input, const MatrixCOO &A);
    void fit(std::vector<MatrixD> xs,
             std::vector<MatrixCOO> As,
             std::vector<MatrixD> ys,
             size_t n_epochs, double learning_rate);
};
#endif
