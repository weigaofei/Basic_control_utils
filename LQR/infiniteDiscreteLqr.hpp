#pragma once

#include "riccati/DARE.hpp"
#include "LTISystem.hpp"

template <size_t STATE_DIM, size_t CONTROL_DIM, typename SCALAR>
class InfiniteDisLqr
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    typedef Eigen::Matrix<SCALAR, STATE_DIM, STATE_DIM> state_matrix_t;
    typedef Eigen::Matrix<SCALAR, CONTROL_DIM, CONTROL_DIM> control_matrix_t;
    typedef Eigen::Matrix<SCALAR, CONTROL_DIM, STATE_DIM> control_state_matrix_t;
    typedef Eigen::Matrix<SCALAR, STATE_DIM, CONTROL_DIM> control_gain_matrix_t;
    typedef Eigen::Matrix<SCALAR, CONTROL_DIM, STATE_DIM> control_feedback_t;

    state_matrix_t compute(const state_matrix_t &Q,
                           const control_matrix_t &R,
                           const state_matrix_t &A,
                           const control_gain_matrix_t &B,
                           control_feedback_t &K,
                           bool verbose = false,
                           const SCALAR eps = 1e-2,
                           size_t maxIter = 2000)
    {
        LTISystem<STATE_DIM, CONTROL_DIM, STATE_DIM, SCALAR> lti_sys(A, B, Eigen::Matrix<SCALAR, STATE_DIM, STATE_DIM>::Zero());
        bool stabilizable = lti_sys.isStabilizable(lti_sys.A(), lti_sys.B());
        if (!stabilizable)
            std::cout << "System is not stabilizable, can not calculate K" << std::endl;
        return dare_.computeSteadyStateRiccatiMatrix(Q, R, A, B, K, verbose, eps, maxIter);
    }

private:
    DARE<STATE_DIM, CONTROL_DIM, SCALAR> dare_; //discrete time riccati equations
};
