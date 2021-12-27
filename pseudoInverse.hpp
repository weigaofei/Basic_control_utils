#include <Eigen/LU>
#include <Eigen/SVD>
#include "cppTypes.hpp"

template <typename T>
void pseudoInverse(templateDMat<T> &matrix, T sigmaThreshold, templateDMat<T> &invMatrix)
{
    if ((1 == matrix.rows()) && (1 == matrix.cols()))
    {
        invMatrix.resize(1, 1);
        if (matrix.coeff(0, 0) > sigmaThreshold)
        {
            invMatrix.coeffRef(0, 0) = 1.0 / matrix.coeff(0, 0);
        }
        else
        {
            invMatrix.coeffRef(0, 0) = 0.0;
        }
        return;
    }

    Eigen::JacobiSVD<templateDMat<T>> svd(matrix, Eigen::ComputeThinU | Eigen::ComputeThinV);
    // not sure if we need to svd.sort()... probably not
    int const nrows(svd.singularValues().rows());
    templateDMat<T> invS;
    invS = templateDMat<T>::Zero(nrows, nrows);
    for (int ii(0); ii < nrows; ++ii)
    {
        if (svd.singularValues().coeff(ii) > sigmaThreshold)
        {
            invS.coeffRef(ii, ii) = 1.0 / svd.singularValues().coeff(ii);
        }
        else
        {
            invS.coeffRef(ii, ii) = 1.0 / sigmaThreshold;
            printf("sigular value is too small: %f\n", svd.singularValues().coeff(ii));
        }
    }
    invMatrix = svd.matrixV() * invS * svd.matrixU().transpose();
}