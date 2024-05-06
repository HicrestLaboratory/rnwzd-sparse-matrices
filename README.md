# Sparse Matrices and GNN from scratch

This repo contains code that implements common functions on dense matrices and sparse matrices in COO format. (matrix.cpp and matrix.hpp)

On top of that the repo contains code that impements end-to-end GNN training and inference. (gnn.cpp and gnn.hpp)

To compile run (on linux)
```
./make.sh
```
To run training and inference run
```
bin/main
```
## Data: Karate Club

I have use the famous Zachary's Karate Club as an example. The data is for binary classification.

To use another dataset you need to change the main file.

## main

The main file reads data, creates a simple gnn, trains it and then evaluates it on train data.

## Limits

- Only Gradient descent is implemented (no SGD, no ADAM etc.)
- No traint/validation/test splitting
- Only two types of loss (MSE, BCE)
- Only one type of layer and activation
- Somewhat incomplete interface of Matrix classes
- Lack of documentation
- Needs cleaning up

# WARNING

This code was made as code assignment and it is not user-friendly, however the code is well organized and quite self explanatory.
