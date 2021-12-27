
#include <iostream>
#include "linearizer.hpp"
#include "LTISystem.hpp"

using namespace std;

int main(void)
{
    const size_t stateDim = 2;
    const size_t controlDim = 1;
    const size_t outputDim = 2;
    Eigen::Matrix<double, stateDim, stateDim> A, dfdx;
    Eigen::Matrix<double, stateDim, controlDim> B, dfdu;
    Eigen::Matrix<double, outputDim, stateDim> C;
    Eigen::Matrix<double, stateDim, 1> curr_x, new_x, constTermC;
    Eigen::Matrix<double, controlDim, 1> curr_u;
    A.setZero();
    B.setZero();
    C.setZero();
    A << 4, 1, 2, 0;
    B << 0, 1;
    C(0, 0) = 1;
    C(0, 1) = 3;
    constTermC << 0, 2;
    curr_x << 1, 2;
    curr_u << 1;
    LTISystem<stateDim, controlDim, outputDim, double> lti_sys(A, B, C);
    bool is_stable = lti_sys.isStable(lti_sys.A(), true);
    bool is_controllable = lti_sys.isControllable(lti_sys.A(), lti_sys.B());
    bool is_oberser = lti_sys.isOberservable(lti_sys.A(), lti_sys.C());
    double dis_timestep = 0.025;
    lti_sys.discreteSystem(dis_timestep, DiscreteMethod::FORWARD_EULER, false);
    cout << "stable:  " << is_stable << "\n"
         << "controllable:  " << is_controllable << "\n"
         << "observe:  " << is_oberser << "\n"
         << "discrete A matrix:\n"
         << lti_sys.disA() << "\n"
         << "discrete B:\n " << lti_sys.disB() << endl;
    lti_sys.updateSystem(A, B, C, constTermC);
    lti_sys.discreteSystem(dis_timestep, DiscreteMethod::EIGEN_EXP, true);
    cout << lti_sys.disA() << "\n"
         << lti_sys.disB() << "\n"
         << lti_sys.disConstTermC() << endl;

    LinerizerNumDiff<stateDim, controlDim, double> linear(
        std::bind(&LTISystem<stateDim, controlDim, outputDim, double>::computeConDynamics,
                  &lti_sys,
                  std::placeholders::_1,
                  std::placeholders::_2,
                  std::placeholders::_3,
                  std::placeholders::_4),
        true);
    dfdx = linear.getDerivativeState(curr_x, curr_u);
    dfdu = linear.getDerivativeControl(curr_x, curr_u);

    cout << "df/dx\n"
         << dfdx << "\ndf/du\n"
         << dfdu << std::endl;

    return 0;
}