#include <iostream>
#include <algorithm>
#include <cassert>
#include <list>
#include <vector>
#include <map>
#include <random>
#include <cmath>
#include <fstream>

#include "csv.h"

#include "matrix.hpp"
#include "gnn.hpp"

int main()
{

  std::string filename = "data/dummy.edgelist";
  bool weighted = false;
  bool directed = false;
  auto A = mcoo_from_el_file(filename, weighted, directed);

  // std::cout << A;

  io::CSVReader<5> in("data/dummy.csv");
  in.read_header(io::ignore_extra_column, "node", "member", "instructor", "administrator", "community");
  int node;
  double member, instructor, administrator, community;
  MatrixCOO x{0, 0}, y{0, 0};
  while (in.read_row(node, member, instructor, administrator, community))
  {
    x.add_row({member, instructor, administrator});
    y.add_row({community, 1 - community});
  }
  // x = x.t(); // TODO
  // y = y.t(); // TODO

  const auto n_samples{1};
  const auto n_v{x.m};    // number of vertices
  const auto n_h_in{x.n}; // dimension of input latent vector
  const auto n_h_hidden{4};
  const auto n_h_out{y.n}; // dimension of output latent vector

  const auto n_epochs{10000};
  const double learning_rate{0.1}; // TODO

  std::random_device seed;
  std::mt19937 rand_gen(seed());
  std::normal_distribution Normal(0.0, 2.0);
  auto normal_gen = [&](auto i, auto j)
  { return Normal(rand_gen); };

  // GNNLayer gnni{n_h_in, n_h_hidden, normal_gen};
  // ActLayer acti{};

  // // GNNLayer gnn1{n_h_hidden, n_h_hidden, normal_gen};
  // //ActLayer act1{};

  // GNNLayer gnnf{n_h_hidden, n_h_out, normal_gen};
  // ActLayer actf{};
  // Network net{{&gnni, &acti,
  // //   &gnn1, &act1,
  //   &gnnf, &actf}};

  GNNLayer gnn{n_h_in, n_h_out, normal_gen};
  ActLayer act{};
  Network net{{&gnn, &act}};

  net.fit({x}, {A}, {y}, // TODO
          n_epochs, learning_rate);

  std::cout << net.forward(x, A);
  std::cout << MatrixD{y};

  return 0;
}
