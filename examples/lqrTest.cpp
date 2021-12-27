#include <iostream>
#include "LQR/infiniteContinousLqr.hpp"
#include "LQR/infiniteDiscreteLqr.hpp"
#include "LQR/finiteDisLqr.hpp"

using namespace std;
using Eigen::Matrix;
using Eigen::MatrixXf;

int main(void)
{
    const size_t stateDim = 2;
    const size_t controlDim = 1;

    Eigen::Matrix<float, stateDim, stateDim> A, Q;
    Eigen::Matrix<float, stateDim, controlDim> B;
    Eigen::Matrix<float, controlDim, controlDim> R;
    Eigen::Matrix<float, controlDim, stateDim> K_dis, K_con;
    Eigen::Matrix<float, stateDim, stateDim> Pn_dis, Pn_con;

    A.setZero();
    Q.setZero();
    B.setZero();
    R.setZero();

    A << 1, 2, 3, 4;
    B << 1, 10;
    Q.diagonal() << 100, 1;
    R.diagonal() << 0.001;
    InfiniteDisLqr<stateDim, controlDim, float> indislqr;
    Pn_dis = indislqr.compute(Q, R, A, B, K_dis);

    InfiniteConLqr<stateDim, controlDim, float> inconlqr;
    Pn_con = inconlqr.compute(Q, R, A, B, K_con);

    cout << "discrete lqr\n"
         << "K\n"
         << K_dis << "\n"
         << "Pn\n"
         << Pn_dis << "\n\n"
         << "continous lqr\n"
         << "K\n"
         << K_con << "\n"
         << "Pn\n"
         << Pn_con << "\n\n";
}