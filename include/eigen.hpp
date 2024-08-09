/**
 * @file eigen.hpp
 * @brief A class to organize the eigenvalue and eigenvector calculations for 3x3 matrices.
 */

#ifndef MY_EIGEN_HEADER
#define MY_EIGEN_HEADER
/**
 * @class eigen
 * @brief Wrapper class of the eigen-system.
 *
 */
class eigen
{
public:
    // calculate the eigenvalues and eigenvectors of a given 3x3 symmetric matrix.
    static void eigens_sym_33( double matrixData[ 9 ], double eigenValues[ 3 ],
                               double eigenVectors[ 9 ] );
};

#endif
