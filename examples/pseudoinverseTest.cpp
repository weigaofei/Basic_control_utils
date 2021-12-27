#include <iostream>
#include "pseudoInverse.hpp"
#include "Timer.hpp"

using namespace std;
using Eigen::Matrix;
using Eigen::MatrixXd;

int main()
{
    MatrixXd a = MatrixXd::Random(20, 24);
    MatrixXd a_inv_svd, a_inv_common;
    Timer t1;
    pseudoInverse(a, 1e-3, a_inv_svd);
    cout << "SVD's pseudo inverse calcu time  " << t1.getMs() << "  ms" << endl;
    Timer t2;
    a_inv_common = a.transpose() * (a * a.transpose()).inverse();
    cout << "common's pseudo inverse calcu time  " << t2.getMs() << "  ms" << endl;
    double diff = (a_inv_svd - a_inv_common).norm();
    cout << "two method's pseudo inverse matrix's diff norm  " << diff << endl;

    return 0;
}