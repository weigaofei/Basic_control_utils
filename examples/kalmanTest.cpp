#include <iostream>
#include "Filter/kalmanFilter.hpp"

using namespace std;

int main(void)
{
    const size_t stateDim = 2;
    const size_t controlDim = 1;
    const size_t outputDim = 1;

    Eigen::Matrix<double, stateDim, stateDim> A, Q;
    Eigen::Matrix<double, stateDim, controlDim> B;
    Eigen::Matrix<double, outputDim, stateDim> C;
    Eigen::Matrix<double, controlDim, controlDim> R;
    Eigen::Matrix<double, stateDim, 1> curr_x, final_x, process_x;
    Eigen::Matrix<double, controlDim, 1> curr_u;
    Eigen::Matrix<double, outputDim, 1> measured_y;
    A << 1, 2, 4, 5;
    B << 1, 0;
    C << 1, 0;
    Q << 10, 0, 0, 100;
    R << 0.003;
    curr_x << 1, 4;
    curr_u << 3;
    measured_y << 13;

    kalmanFilter<stateDim, controlDim, outputDim, double> kafilter(A, B, C);
    kafilter.setNoiseMatrix(Q, R);
    kafilter.predict(curr_x, curr_u, 0);
    kafilter.measureUpdate(measured_y);
    kafilter.getFinalPreValue(final_x);
    kafilter.getProcessPreValue(process_x);
    cout << "Ax+Bu esimate\n"
         << process_x << "\n\n";
    cout << "first filter\n"
         << final_x << "\n\n";
    kafilter.predict(curr_x, curr_u, 0);
    kafilter.measureUpdate(measured_y);
    kafilter.getFinalPreValue(final_x);
    cout << "second filter\n"
         << final_x << "\n";
}