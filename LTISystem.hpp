#ifndef _LTISYSTEM_H_
#define _LTISYSTEM_H_

#include "cppTypes.hpp"
#include <unsupported/Eigen/MatrixFunctions>

enum DiscreteMethod{FORWARD_EULER = 0, BACKWARD_EULER, TUSTIN, EIGEN_EXP};

template<size_t STATE_DIM, size_t CONTROL_DIM, size_t OUTPUT_DIM, typename SCALAR>
class LTISystem{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    using stateVector = typename Eigen::Matrix<SCALAR, STATE_DIM, 1>;
    using controlVector = typename Eigen::Matrix<SCALAR, CONTROL_DIM, 1>;
    using outputVector = typename Eigen::Matrix<SCALAR, OUTPUT_DIM, 1>;
    using stateMatrix = typename Eigen::Matrix<SCALAR, STATE_DIM, STATE_DIM>;
    using controlMatrix = typename Eigen::Matrix<SCALAR, STATE_DIM, CONTROL_DIM>;
    using outputMatrix = typename Eigen::Matrix<SCALAR, OUTPUT_DIM, STATE_DIM>;

    LTISystem(const stateMatrix& statematrix, const controlMatrix& controlmatrix, const outputMatrix& outputmatrix):
    A_(statematrix), B_(controlmatrix), C_(outputmatrix)
    {
        constTermC_.setZero();
        disA_.setZero();
        disB_.setZero();
        disConstTermC_.setZero();
    }
    ~LTISystem(){}

    void computeConDynamics(const stateVector& curr_state, const controlVector& curr_control, 
        const SCALAR& curr_time, stateVector& derived_state)
    {
        derived_state = A_ * curr_state + B_ * curr_control;
    }

    void computeDisDynamics(const stateVector& curr_state, const controlVector& curr_control, 
        const SCALAR& curr_time, stateVector& derived_state)
    {
        derived_state = disA_ * curr_state + disB_ * curr_control;
    }

    void computeOutput(const stateVector& curr_state, const controlVector& curr_control, 
        const SCALAR& curr_time, outputVector& derived_out)
    {
        stateVector new_state;
        computeDisDynamics(curr_state, curr_control, curr_time, new_state);
        derived_out = C_ * new_state;
    }
    /*
    For marginally stable systems, starting from any initial conditions, the solution will neither
    converge to zero nor diverge to infinity. This also calls stabilty in the sense of Lyapunov.
    */
    bool isStable(const stateMatrix& statematrix, bool isConSystem){
        bool flag;
        if(isConSystem){
            Eigen::EigenSolver<Eigen::Matrix<double, STATE_DIM, STATE_DIM>> es(statematrix.template cast<double>());
            Eigen::MatrixXcd evals = es.eigenvalues();
            Eigen::MatrixXd evalsReal;
            evalsReal = evals.real();
            bool allSmallerZero = (evalsReal.array() < 0).all();
            bool existEqualZero = (evalsReal.array() == 0).any();
            bool existGreaterZero = (evalsReal.array() > 0).any();
            
            if(allSmallerZero) flag = true;
            else if(existEqualZero && !existGreaterZero){
                throw std::runtime_error("exist zero eigenvalues and no positive values, system may be marginal stable(or\
                refer as Lyapunov stable). If all the Jordan blocks associated with zero eigenvalues have size one, system\
                is marginal stable. please use MATLAB to check: jordan(A)");
                }
            else if(existGreaterZero) flag = false;
        }
        if(!isConSystem){
            Eigen::EigenSolver<Eigen::Matrix<double, STATE_DIM, STATE_DIM>> es(statematrix.template cast<double>());
            Eigen::MatrixXcd evals = es.eigenvalues();
            Eigen::MatrixXd evalsAbs;
            evalsAbs = evals.array().abs();
            bool allSmallerOne = (evalsAbs.array() < 0).all();
            bool existEqualOne = (evalsAbs.array() == 0).any();
            bool existGreaterOne = (evalsAbs.array() > 0).any();
            
            if(allSmallerOne) flag = true;
            else if(existEqualOne && !existGreaterOne){
                throw std::runtime_error("exist eigenvalues whose abs equal one and no greater one, system may be marginal stable(or\
                refer as Lyapunov stable). If all the Jordan blocks associated with such eigenvalues have size one, system\
                is marginal stable. please use MATLAB to check: jordan(A)");
                }
            else if(existGreaterOne) flag = false;
        }

        return flag;
    }

    bool isControllable(const stateMatrix& statematrix, const controlMatrix& controlmatrix){
        Eigen::Matrix<double, STATE_DIM, STATE_DIM * CONTROL_DIM> Con;
        Con.block(0, 0, STATE_DIM, CONTROL_DIM) = controlmatrix.template cast<double>();
        for(size_t i = 1; i < STATE_DIM; ++i){
            Con.block(0, i * CONTROL_DIM, STATE_DIM, CONTROL_DIM) = 
                statematrix.template cast<double>() * Con.block(0, (i-1) * CONTROL_DIM, STATE_DIM, CONTROL_DIM);
        }
        Eigen::FullPivLU<Eigen::Matrix<double, STATE_DIM, STATE_DIM * CONTROL_DIM>> LUdecomposition(Con);

        return LUdecomposition.rank() == STATE_DIM;
    }

    bool isOberservable(const stateMatrix& statematrix, const outputMatrix& outputmatrix){
        Eigen::Matrix<double, OUTPUT_DIM * STATE_DIM, STATE_DIM> Ob;
        Ob.block(0, 0, OUTPUT_DIM, STATE_DIM) = outputmatrix.template cast<double>();
        for (size_t i = 1; i < STATE_DIM; i++)
        {
            Ob.block(i * OUTPUT_DIM, 0, OUTPUT_DIM, STATE_DIM) =
                Ob.block((i - 1) * OUTPUT_DIM, 0, OUTPUT_DIM, STATE_DIM) * statematrix.template cast<double>();
        }
        Eigen::FullPivLU<Eigen::Matrix<double, OUTPUT_DIM * STATE_DIM, STATE_DIM>> LUdecomposition(Ob);
        return LUdecomposition.rank() == STATE_DIM;
    }
    /*
    LTI system is stabilizable if there exists matrix K such that system xdot = (A-BK)x is stable
    for LTI system, controllabe equals stabilizable
    so when LTI system is controllable, LQR have an unique solution, which can be solved using backwards method
    */
    bool isStabilizable(const stateMatrix& statematrix, const controlMatrix& controlmatrix){
        return isControllable(statematrix, controlmatrix);
    }

    void updateSystem(const stateMatrix& statematrix, const controlMatrix& controlmatrix, const outputMatrix& outputmatrix, 
        const stateVector& constterm)
    {
        A_  = statematrix;
        B_ = controlmatrix;
        C_ =  outputmatrix;
        constTermC_ = constterm;
    }

    stateMatrix A() const {return A_;}
    controlMatrix B() const {return B_;}
    stateVector constTermC() const {return constTermC_;}
    stateMatrix disA() const {return disA_;}
    controlMatrix disB() const {return disB_;}
    stateVector disConstTermC() const {return disConstTermC_;}
    outputMatrix C() const{return C_;}

    void discreteSystem(const SCALAR& time_step, DiscreteMethod method, bool haveConstTerm){
        if(!haveConstTerm){
            switch(method){
                case DiscreteMethod::FORWARD_EULER:
                {
                    disA_ = stateMatrix::Identity() + time_step * A_;
                    disB_ = time_step * B_;
                    break;
                }
                case DiscreteMethod::EIGEN_EXP:
                {
                    Eigen::Matrix<SCALAR, STATE_DIM + CONTROL_DIM, STATE_DIM + CONTROL_DIM> AB_, expm;
                    AB_.setZero();
                    AB_.block(0, 0, STATE_DIM, STATE_DIM) = A_;
                    AB_.block(0, STATE_DIM, STATE_DIM, CONTROL_DIM) = B_;
                    AB_ *= time_step;
                    expm = AB_.exp(); 
                    disA_ = expm.block(0, 0, STATE_DIM, STATE_DIM);
                    disB_ = expm.block(0, STATE_DIM, STATE_DIM, CONTROL_DIM) ;
                    break;
                }
                case DiscreteMethod::BACKWARD_EULER:
                {
                    disA_ = (stateMatrix::Identity() - time_step * A_).inverse();
                    disB_ = time_step * disA_ * B_;
                    break;
                }
                case DiscreteMethod::TUSTIN:
                {
                    disA_ = (stateMatrix::Identity() + time_step * A_ / 2) * (stateMatrix::Identity() - time_step * A_ / 2).inverse();
                    disB_ = (stateMatrix::Identity() - time_step * A_ / 2).inverse() * B_ * std::sqrt(time_step);
                    break;
                }
                default:
                    throw std::runtime_error("unknown discrete method");
            }
        }
        else if(haveConstTerm){
            std::cout << "In this case, will use  Eigen's exp() to calculate discrete matrix" << std::endl;
            Eigen::Matrix<SCALAR, STATE_DIM + CONTROL_DIM + 1, STATE_DIM + CONTROL_DIM + 1> ABC_, expm;
            ABC_.setZero();
            ABC_.block(0, 0, STATE_DIM, STATE_DIM) = A_;
            ABC_.block(0, STATE_DIM, STATE_DIM, CONTROL_DIM) = B_;
            ABC_.block(0, STATE_DIM + CONTROL_DIM, STATE_DIM, 1) = constTermC_;
            ABC_ *= time_step;
            expm = ABC_.exp();
            
            disA_ = expm.block(0, 0, STATE_DIM, STATE_DIM) ;
            disB_ = expm.block(0, STATE_DIM, STATE_DIM, CONTROL_DIM) ;
            disConstTermC_ = expm.block(0, STATE_DIM + CONTROL_DIM, STATE_DIM, 1);
        }
    }

private:
    stateMatrix A_;
    controlMatrix B_;
    outputMatrix C_;//y = Cx
    stateVector constTermC_;//xdot = A * x + B * u + constTermC_ 
    stateMatrix disA_;//discrete
    controlMatrix disB_;
    stateVector disConstTermC_;

};

#endif