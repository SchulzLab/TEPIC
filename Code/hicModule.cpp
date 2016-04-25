#include <iostream>
#include "readMatrix.hpp"
#include <Eigen/Dense>

using namespace std;
using Eigen::MatrixXd;

void testFileRead(char* file){

    FILE* pFile = fopen (file, "r");
    SparseMatrix<float> mat = readMatrix(pFile);

    cout << "Width of x-dimension: " << mat.cols() << endl;
    cout << "Height of y-dimension: " << mat.rows() << endl;

    fclose(pFile);
}

// Use this to test Eigen on deployment machines
void runEigenExample(){
    cout << "Testing EIGEN..." << endl;
    MatrixXd m(2,2);
    m(0,0) = 3;
    m(1,0) = 2.5;
    m(0,1) = -1;
    m(1,1) = m(1,0) + m(0,1);
    cout << m << std::endl;
    cout << "Finished." << endl;
}

int main(int argc, char* argv[]) {

    cout << "Running the TEPIC-HIC module..." << endl;

    testFileRead(argv[1]);
    runEigenExample();

    return 0;
}

