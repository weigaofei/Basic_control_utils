#include "cppTypes.hpp"
#include "orientation.hpp"

using Eigen::Matrix;
using Eigen::Quaternion;
using std::cout;
using std::endl;

int main(void)
{
    Matrix<double, 3, 1> rpy_d, rpy_slope, rpy_result;
    Matrix<double, 3, 3> R0, R_slope;
    Quaternion<double> quat;
    Matrix<double, 4, 1> quaVec, quaVec2;
    rpy_d << 0, -0.2, 0;
    rpy_slope << 0, -0.4, 0;
    R0 = coordinateRotation<double>(CoordinateAxis::Z, rpy_d[2]) *
         coordinateRotation<double>(CoordinateAxis::Y, rpy_d[1]) *
         coordinateRotation<double>(CoordinateAxis::X, rpy_d[0]);
    R_slope = coordinateRotation<double>(CoordinateAxis::Z, rpy_slope[2]) *
              coordinateRotation<double>(CoordinateAxis::Y, rpy_slope[1]) *
              coordinateRotation<double>(CoordinateAxis::X, rpy_slope[0]);
    rpy_result = rotationMatrixToRPY(R_slope * R0);

    cout << rpy_result.transpose() << endl;

    return 0;
}