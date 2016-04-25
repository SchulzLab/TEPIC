//
// Created by Fabian Kern on 25.04.16.
//

#ifndef TEPIC_HIC_READMATRIX_HPP
#define TEPIC_HIC_READMATRIX_HPP

#include <cstdio>
#include <iostream>
#include <Eigen/Sparse>

using namespace std;
using Eigen::SparseMatrix;

SparseMatrix<float> readMatrix(FILE* file);

#endif //TEPIC_HIC_READMATRIX_HPP
