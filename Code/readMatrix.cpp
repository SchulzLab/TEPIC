//
// Created by Fabian Kern on 25.04.16.
//

#include "readMatrix.hpp"

typedef Eigen::Triplet<double> T;

Eigen::SparseMatrix<float> readMatrix(FILE* file){
    if(file){
        int i,j = 0;
        float v_ij = 0.0;
        std::vector<T> tripletList;
        int maxi, maxj = 0;

        while(fscanf (file, "%i %i %f", &i, &j, &v_ij) != EOF){
            // read file content into memory using EIGEN
            tripletList.push_back(T(i,j,v_ij));
            maxi = i > maxi ? i : maxi;
            maxj = j > maxj ? j : maxj;
        }

        Eigen::SparseMatrix<float> mat(maxi+1,maxj+1); // +1 because EIGEN may throwing an assertion if the matrix is dense
        mat.setFromTriplets(tripletList.begin(), tripletList.end());

        return mat;
    }
    else{
        cout << "Could not open file!" << endl;
        Eigen::SparseMatrix<float> mat;
        return mat;
    }
}