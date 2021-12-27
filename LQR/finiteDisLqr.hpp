#pragma once

#include "riccati/DARE.hpp"
#include "cppTypes.hpp"
#include <vector>

//Finite discrete lqr often used for local trajectory stabilization for nonlinear systems
//see more http://underactuated.mit.edu/trajopt.html#section4

template <size_t STATE_DIM, size_t CONTROL_DIM, typename SCALAR>
class FinieDisLqr
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    using state_vector_t = typename Eigen::Matrix<SCALAR, STATE_DIM, 1>;
    using control_vector_t = typename Eigen::Matrix<SCALAR, CONTROL_DIM, 1>;
    using state_matrix_t = typename Eigen::Matrix<SCALAR, STATE_DIM, STATE_DIM>;
    using control_matrix_t = typename Eigen::Matrix<SCALAR, STATE_DIM, CONTROL_DIM>;
    using state_weight_t = typename Eigen::Matrix<SCALAR, STATE_DIM, STATE_DIM>;
    using control_weight_t = typename Eigen::Matrix<SCALAR, CONTROL_DIM, CONTROL_DIM>;
    using control_feedback_t = typename Eigen::Matrix<SCALAR, CONTROL_DIM, STATE_DIM>;
    using state_vector_array_t = std::vector<state_vector_t>;
    using control_vector_array_t = std::vector<control_vector_t>;
    using state_matrix_array_t = std::vector<state_matrix_t>;
    using control_matrix_array_t = std::vector<control_matrix_t>;
    using control_feedback_array_t = std::vector<control_feedback_t>;

    /*
    One of the most powerful applications of time-varying LQR involves linearizing around a reference 
    trajectory of a nonlinear system and using LQR to provide a trajectory controller. This is very similar 
    to using LQR to stablize a fixed-point, but for a trajectory, the linearization is time-varying.
    final controller: U(t) = Uref(t) - K(t)(X(t) - Xref(t));
    see more http://underactuated.mit.edu/lqr.html
    */
    void designController(const state_vector_array_t &x_trajectory,
                          const control_vector_array_t &u_trajectory,
                          const state_matrix_array_t &A,
                          const control_matrix_array_t &B,
                          const state_weight_t &Q,
                          const control_weight_t &R,
                          const state_weight_t &terminalWeight,
                          control_feedback_array_t &K)
    {
        size_t N = x_trajectory.size() - 1;
        state_matrix_t P = terminalWeight;
        for (int i = N - 1; i >= 0; i--)
        {
            dare_.iterateRobust(Q, R, A[i], B[i], P, K[i]);
        }
    }

    // convenient to opreate two trajectories, sunch as + - *
    // +  -  *  /  %
    enum TrajOpretion
    {
        ADD = 0,
        SUB,
        MUL,
        DIV,
        MODULUS
    };

    void trajCalcu(const state_vector_array_t &a, const state_vector_array_t &b, const TrajOpretion &op,
                   state_vector_array_t &result)
    {
        result.resize(a.size());
        switch (op)
        {
        case TrajOpretion::ADD:
            std::transform(a.begin(), a.end(), b.begin(), result.begin(), std::plus<state_matrix_t>());
            break;
        case TrajOpretion::SUB:
            std::transform(a.begin(), a.end(), b.begin(), result.begin(), std::minus<state_matrix_t>());
            break;
        case TrajOpretion::MUL:
            std::transform(a.begin(), a.end(), b.begin(), result.begin(), std::multiplies<state_matrix_t>());
            break;
        // case TrajOpretion::DIV:
        //     std::transform(a.begin(), a.end(), b.begin(), result.begin(), std::divides<state_matrix_t>());
        //     break;
        // case TrajOpretion::MODULUS:
        //     std::transform(a.begin(), a.end(), b.begin(), result.begin(), std::modulus<state_matrix_t>());
        //     break;
        default:
            throw std::runtime_error("Unknown opration method");
        }
    }

private:
    DARE<STATE_DIM, CONTROL_DIM, SCALAR> dare_;
};