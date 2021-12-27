/**********************************************************************************************************************
This file is part of the Control Toolbox (https://github.com/ethz-adrl/control-toolbox), copyright by ETH Zurich.
Licensed under the BSD-2 license (see LICENSE file in main directory)
**********************************************************************************************************************/

#pragma once

#include <iostream>
#include "cppTypes.hpp"
/*!
 * \ingroup LQR
 *+
 * \brief Continuous-Time Algebraic Riccati Equation
 *
 * solves the continuous-time Infinite-Horizon Algebraic Riccati Equation
 *
 * @tparam STATE_DIM system state dimension
 * @tparam CONTROL_DIM system control input dimension
 */
template <size_t STATE_DIM, size_t CONTROL_DIM, typename SCALAR>
class CARE
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    typedef Eigen::Matrix<SCALAR, STATE_DIM, STATE_DIM> state_matrix_t;
    typedef Eigen::Matrix<SCALAR, CONTROL_DIM, CONTROL_DIM> control_matrix_t;
    typedef Eigen::Matrix<SCALAR, CONTROL_DIM, STATE_DIM> control_state_matrix_t;
    typedef Eigen::Matrix<SCALAR, STATE_DIM, CONTROL_DIM> control_gain_matrix_t;
    typedef Eigen::Matrix<SCALAR, CONTROL_DIM, STATE_DIM> control_feedback_t;

    typedef Eigen::Matrix<SCALAR, 2 * STATE_DIM, 2 * STATE_DIM> schur_matrix_t;
    typedef Eigen::Matrix<SCALAR, 2 * STATE_DIM, STATE_DIM> factor_matrix_t;

    CARE(){};

    // Implementation using the Schur-Method
    // This is numerically more stable and should be preferred over the naive implementation
    bool solve(const state_matrix_t &Q,
               const control_matrix_t &R,
               const state_matrix_t &A,
               const control_gain_matrix_t &B,
               state_matrix_t &P,
               bool RisDiagonal)
    {
        control_matrix_t R_inverse;
        if (RisDiagonal)
        {
            R_inverse.setZero();
            R_inverse.diagonal().noalias() = R.diagonal().cwiseInverse();
        }
        else
        {
            R_inverse.noalias() = R.inverse();
        }

        schur_matrix_t M;
        M << A, -B * R_inverse * B.transpose(), -Q, -A.transpose();

        return solveSchurIterative(M, P);
    }

    state_matrix_t computeSteadyStateRiccatiMatrix(const state_matrix_t &Q,
                                                   const control_matrix_t &R,
                                                   const state_matrix_t &A,
                                                   const control_gain_matrix_t &B,
                                                   const bool RisDiagonal = false)
    {
        state_matrix_t P;
        control_matrix_t Rinv;

        solve(Q, R, A, B, P, RisDiagonal);

        return P;
    }

private:
    bool solveSchurIterative(const schur_matrix_t &M,
                             state_matrix_t &P,
                             SCALAR epsilon = 1e-6,
                             int maxIterations = 100000)
    {

        bool converged = false;

        schur_matrix_t Mlocal = M;

        int iterations = 0;
        while (!converged)
        {
            if (iterations > maxIterations)
                return false;

            schur_matrix_t Mdiff = Mlocal - Mlocal.inverse();

            schur_matrix_t Mnew = Mlocal - 0.5 * Mdiff;

            converged = Mnew.isApprox(Mlocal, epsilon);

            Mlocal = Mnew;

            iterations++;
        }

        /* break down W and extract W11 W12 W21 W22  (what is the size of these?) */
        state_matrix_t M11(Mlocal.template block<STATE_DIM, STATE_DIM>(0, 0));
        state_matrix_t M12(Mlocal.template block<STATE_DIM, STATE_DIM>(0, STATE_DIM));
        state_matrix_t M21(Mlocal.template block<STATE_DIM, STATE_DIM>(STATE_DIM, 0));
        state_matrix_t M22(Mlocal.template block<STATE_DIM, STATE_DIM>(STATE_DIM, STATE_DIM));

        /* find M and N using the elements of W	 */
        factor_matrix_t U;
        factor_matrix_t V;

        U.template block<STATE_DIM, STATE_DIM>(0, 0) = M12;
        U.template block<STATE_DIM, STATE_DIM>(STATE_DIM, 0) = M22 + state_matrix_t::Identity();

        V.template block<STATE_DIM, STATE_DIM>(0, 0) = M11 + state_matrix_t::Identity();
        V.template block<STATE_DIM, STATE_DIM>(STATE_DIM, 0) = M21;

        /* Solve for S from the equation   MS=N */
        FullPivLU_.compute(U);

        P = FullPivLU_.solve(-V);

        return true;
    }

    Eigen::FullPivLU<factor_matrix_t> FullPivLU_;
};
