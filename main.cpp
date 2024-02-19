#include <iostream>
#include <algorithm>
#include <cassert>
#include <list>
#include <vector>
#include <map>
#include <random>
#include <cmath>
#include <fstream>

#include "matrix.hpp"
#include "gnn.hpp"

int main()
{

    // const auto n_samples{1};
    // const auto n_v{1};  // number of vertices
    // const auto n_h_in{1};  // dimension of input latent vector
    // const auto n_h_out{1}; // dimension of output latent vector

    // const auto n_epochs{100};
    // const double learning_rate{-0.01};

    // std::random_device seed;
    // std::mt19937 rand_gen(seed());
    // std::normal_distribution Normal(0.0, 2.0);
    // auto normal_gen = [&](auto i, auto j)
    // { return Normal(rand_gen); };

    // GNNLayer gnn1{n_h_in, n_h_out, normal_gen};
    // ReLULayer relu1{};
    // Network net1{{&gnn1, &relu1}};

    // std::vector<MatrixD> xs{};
    // std::vector<MatrixCOO> As{};
    // std::vector<MatrixD> ys{};
    // for (size_t i{}; i < n_samples; ++i)
    // {
    //     MatrixD x{n_h_in, n_v, normal_gen};
    //     MatrixCOO A{n_v, n_v, normal_gen};
    //     xs.push_back(x);
    //     As.push_back(A);
    //     ys.push_back(net1.forward(x, A));
    // }

    // GNNLayer gnn2{n_h_in, n_h_out, normal_gen};
    // ReLULayer relu2{};
    // Network net2{{&gnn2, &relu2}};

    // net2.fit(xs, As, ys, n_epochs, learning_rate);

    std::string filename = "data/karate.edgelist";
    bool weighted = false;
    bool directed = false;
    auto A = mcoo_from_el_file(filename,weighted,directed);
    
    std::cout << A;
    return 0;
}
