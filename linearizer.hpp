#ifndef _LINEARIZER_H_
#define _LINEARIZER_H_

#include <cmath>
#include <iostream>
#include <functional>
#include "cppTypes.hpp"

//Computes the linearization of a system dynamics function through numerical finite differencing
//http://en.wikipedia.org/wiki/Numerical_differentiation#Practical_considerations_using_floating_point_arithmetic

template <size_t STATE_DIM, size_t CONTROL_DIM, typename SCALAR>
class LinerizerNumDiff
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    typedef Eigen::Matrix<SCALAR, STATE_DIM, 1> stateVector;
    typedef Eigen::Matrix<SCALAR, CONTROL_DIM, 1> controlVector;
    typedef Eigen::Matrix<SCALAR, STATE_DIM, STATE_DIM> jacobianStateMatrix;
    typedef Eigen::Matrix<SCALAR, STATE_DIM, CONTROL_DIM> jacobianControlMatrix;

    //void(current_state, current_control, current_time, result)
    typedef std::function<void(const stateVector&, const controlVector&, const SCALAR&, stateVector&)> systemDynamic;

    //dyn can be assignmented by std::bind()
    //https://blog.csdn.net/nanjiye/article/details/52164279
    LinerizerNumDiff(systemDynamic dyn, bool isSymmetricDerivative = true)
        : dynamic_func_(dyn), isSymmetricDerivative_(isSymmetricDerivative)
    {
        dfdx_.setZero();
        dfdu_.setZero();
        eps_ = sqrt(Eigen::NumTraits<SCALAR>::epsilon());
    }

    jacobianStateMatrix getDerivativeState(const stateVector& x, const controlVector& u, const SCALAR t = SCALAR(0))
    {
        if (!isSymmetricDerivative_)
            dynamic_func_(x, u, t, f_ref);

        for (size_t i = 0; i < STATE_DIM; i++)
        {
            SCALAR h = eps_ * std::max(std::abs<SCALAR>(x(i)), SCALAR(1.0));
            SCALAR x_plusDelta = x(i) + h;
            SCALAR delta = h;

            stateVector x_perturbed = x;
            x_perturbed(i) = x_plusDelta;

            // evaluate dynamics at perturbed state
            stateVector f_plusDelta;
            dynamic_func_(x_perturbed, u, t, f_plusDelta);

            if (isSymmetricDerivative_)
            {
                SCALAR x_subDelta = x(i) - h;

                x_perturbed = x;
                x_perturbed(i) = x_subDelta;

                stateVector f_subDelta;
                dynamic_func_(x_perturbed, u, t, f_subDelta);

                dfdx_.col(i) = (f_plusDelta - f_subDelta) / (2 * delta);
            }
            else
            {
                dfdx_.col(i) = (f_plusDelta - f_ref) / delta;
            }
        }

        return dfdx_;
    }

    jacobianControlMatrix getDerivativeControl(const stateVector& x, const controlVector& u, const SCALAR t = SCALAR(0)){
        if (!isSymmetricDerivative_)
            dynamic_func_(x, u, t, f_ref);

        for (size_t i = 0; i < CONTROL_DIM; i++)
        {
            SCALAR h = eps_ * std::max(std::abs<SCALAR>(x(i)), SCALAR(1.0));
            SCALAR u_plusDelta = u(i) + h;
            SCALAR delta = h;

            controlVector u_perturbed = u;
            u_perturbed(i) = u_plusDelta;

            stateVector f_plusDelta;
            dynamic_func_(x, u_perturbed, t, f_plusDelta);

            if (isSymmetricDerivative_)
            {
                SCALAR u_subDelta = u(i) - h;

                u_perturbed = u;
                u_perturbed(i) = u_subDelta;

                stateVector f_subDelta;
                dynamic_func_(x, u_perturbed, t, f_subDelta);

                dfdu_.col(i) = (f_plusDelta - f_subDelta) / (2 * delta);
            }
            else
            {
                dfdu_.col(i) = (f_plusDelta - f_ref) / delta;
            }
        }

        return dfdu_;
    }

protected:
    systemDynamic dynamic_func_;
    bool isSymmetricDerivative_;
    SCALAR  eps_;
    jacobianStateMatrix dfdx_;
    jacobianControlMatrix dfdu_;
    stateVector f_ref;
};


#endif