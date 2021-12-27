#pragma once

#include "riccati/CARE.hpp"
#include "LTISystem.hpp"

/*!
 * \ingroup LQR
 *
 * \brief continuous-time infinite-horizon LQR
 *
 * Implements continous-time infinite-horizon LQR.
 * Resulting feedback law will take the form
 * \f[
 * u_{fb} = -K \cdot (x - x_{ref})
 * \f]
 *
 * @tparam STATE_DIM
 * @tparam CONTROL_DIM
 */
template <size_t STATE_DIM, size_t CONTROL_DIM, typename SCALAR>
class InfiniteConLqr
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    typedef Eigen::Matrix<SCALAR, STATE_DIM, STATE_DIM> state_matrix_t;
    typedef Eigen::Matrix<SCALAR, CONTROL_DIM, CONTROL_DIM> control_matrix_t;
    typedef Eigen::Matrix<SCALAR, CONTROL_DIM, STATE_DIM> control_state_matrix_t;
    typedef Eigen::Matrix<SCALAR, STATE_DIM, CONTROL_DIM> control_gain_matrix_t;
    typedef Eigen::Matrix<SCALAR, CONTROL_DIM, STATE_DIM> control_feedback_t;

    //! design the infinite-horizon LQR controller.
    /*!
	 * @param Q state-weighting matrix
	 * @param R control input weighting matrix
	 * @param A linear system dynamics matrix A
	 * @param B linear system dynamics matrix B
	 * @param K control feedback matrix K (to be designed)
	 * @param RisDiagonal set to true if R is a diagonal matrix (efficiency boost)
	 * @param solveRiccatiIteratively
	 * 	use closed-form solution of the infinite-horizon Riccati Equation
	 * @return success
	 */
    state_matrix_t compute(const state_matrix_t &Q,
                           const control_matrix_t &R,
                           const state_matrix_t &A,
                           const control_gain_matrix_t &B,
                           control_feedback_t &K,
                           bool RisDiagonal = false)
    {
        LTISystem<STATE_DIM, CONTROL_DIM, STATE_DIM, SCALAR> lti_sys(A, B, Eigen::Matrix<SCALAR, STATE_DIM, STATE_DIM>::Zero());
        bool stabilizable = lti_sys.isStabilizable(lti_sys.A(), lti_sys.B());
        if (!stabilizable)
            std::cout << "System is not stabilizable, can not calculate K" << std::endl;

        control_matrix_t R_inverse;
        state_matrix_t P;

        P = care_.computeSteadyStateRiccatiMatrix(Q, R, A, B);

        K = (R_inverse * (B.transpose() * P));

        return P;
    }

private:
    CARE<STATE_DIM, CONTROL_DIM, SCALAR> care_; // continuous-time algebraic riccati equation
};
