#pragma once

#include "cppTypes.hpp"
#include "LTISystem.hpp"

template <size_t STATE_DIM, size_t CONTROL_DIM, size_t OUTPUT_DIM, typename SCALAR>
class kalmanFilter
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    using state_vector_t = typename Eigen::Matrix<SCALAR, STATE_DIM, 1>;
    using control_vector_t = typename Eigen::Matrix<SCALAR, CONTROL_DIM, 1>;
    using output_vector_t = typename Eigen::Matrix<SCALAR, OUTPUT_DIM, 1>;
    using state_matrix_t = typename Eigen::Matrix<SCALAR, STATE_DIM, STATE_DIM>;
    using control_matrix_t = typename Eigen::Matrix<SCALAR, STATE_DIM, CONTROL_DIM>;
    using output_matrix_t = typename Eigen::Matrix<SCALAR, OUTPUT_DIM, STATE_DIM>;
    using process_noise_t = typename Eigen::Matrix<SCALAR, STATE_DIM, STATE_DIM>;
    using measure_noise_t = typename Eigen::Matrix<SCALAR, OUTPUT_DIM, OUTPUT_DIM>;

    kalmanFilter(const state_matrix_t &statematrix, const control_matrix_t &controlmatrix,
                 const output_matrix_t &outputmatrix) : lti_sys(statematrix, controlmatrix, outputmatrix)
    {
        Q_ = process_noise_t::Identity();
        R_ = measure_noise_t::Identity();
        P_ = state_matrix_t::Identity();
    }
    ~kalmanFilter() {}

    //refer to https://zhuanlan.zhihu.com/p/45238681
    void predict(const state_vector_t &curr_state, const control_vector_t &curr_control, const SCALAR &curr_time)
    {
        lti_sys.computeConDynamics(curr_state, curr_control, curr_time, estimate_x_);
        P_ = lti_sys.A() * P_ * lti_sys.A().transpose() + Q_;
        //std::cout << P_ <<"\n\n";
    }

    void measureUpdate(const output_vector_t &curr_measure_y)
    {
        estimate_y_ = lti_sys.C() * estimate_x_;
        delta_y_ = curr_measure_y - estimate_y_;
        S_ = lti_sys.C() * P_ * lti_sys.C().transpose() + R_;
        K_ = P_ * lti_sys.C().transpose() * S_.inverse();
        final_x_ = estimate_x_ + K_ * delta_y_;
        P_ -= (K_ * lti_sys.C() * P_).eval();
    }

    //in practice, mainly adjust Q and R to improve accuracy
    void setNoiseMatrix(const process_noise_t &Q, const measure_noise_t &R)
    {
        Q_ = Q;
        R_ = R;
    }

    void getFinalPreValue(state_vector_t &final)
    {
        final = final_x_;
    }

    void getProcessPreValue(state_vector_t &stage)
    {
        stage = estimate_x_;
    }

private:
    LTISystem<STATE_DIM, CONTROL_DIM, OUTPUT_DIM, SCALAR> lti_sys;
    state_vector_t estimate_x_, final_x_;
    state_matrix_t P_;
    process_noise_t Q_;
    measure_noise_t R_, S_;
    output_vector_t estimate_y_, delta_y_;
    Eigen::Matrix<SCALAR, STATE_DIM, OUTPUT_DIM> K_;
};