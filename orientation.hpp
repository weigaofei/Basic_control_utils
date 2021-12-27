#ifndef _ORIENTATION_H_
#define _ORIENTATION_H_

#include <cmath>
#include <iostream>
#include "cppTypes.hpp"

//mainly reference:
//https://github1s.com/machines-in-motion/kino_dynamic_opt/blob/HEAD/momentumopt/src/momentumopt/utilities/OrientationUtils.cpp
//https://github1s.com/mit-biomimetics/Cheetah-Software/blob/HEAD/common/include/Math/orientation_tools.h

enum CoordinateAxis {X, Y, Z};

template <typename T>
T square(T a) {
  return a * a;
}

template <typename T>
templateMat3<typename T::Scalar> crossMatrix(const Eigen::MatrixBase<T>& v);

template <typename T>
templateMat3<T> coordinateRotation(CoordinateAxis axis, const T& theta);

template <typename T>
templateVec3<typename T::Scalar> quatToRPY(const Eigen::MatrixBase<T>& q);

template <typename T>
templateQuat<typename T::Scalar> rpyToQuat(const Eigen::MatrixBase<T>& rpyAngle);

template <typename T>
templateMat3<typename T::Scalar> quatToRotationMatrix(const Eigen::MatrixBase<T>& q);

template <typename T>
templateVec3<typename T::Scalar> rotationMatrixToRPY(const Eigen::MatrixBase<T>& R);

template <typename T>
templateMat3<typename T::Scalar> bodyToWorldMatrix(const Eigen::MatrixBase<T>& rpyAngle);

template <typename T, typename T2>
templateVec3<typename T::Scalar> eulerAngleRateToAngVelWorld(const Eigen::MatrixBase<T>& rpyAngle,const Eigen::MatrixBase<T2>& rpyAngleRate);

template <typename T>
templateVec3<typename T::Scalar> quaternionToso3(const Eigen::MatrixBase<T>& q);

template <typename T>
templateVec3<typename T::Scalar> so3ToQuaternion(const Eigen::MatrixBase<T>& rotation_vector);

template <typename T>
templateMat4<typename T::Scalar> quatToProductMatrix(const Eigen::MatrixBase<T>& quat, bool isRightProduct);

template <typename T>
templateMat4c3<typename T::Scalar> jacobianQuatRateWrtEulerAngleRate(const Eigen::MatrixBase<T>& rotation_vector);

template <typename T>
templateMat3c4<typename T::Scalar> jacobianAngVelWorldWrtQuatRate(const Eigen::MatrixBase<T>& quat);

template <typename T, typename T2>
templateVec3<typename T::Scalar> eulerAngleRateToAngVelWorldUsingQuat(const Eigen::MatrixBase<T>& rpyAngle,const Eigen::MatrixBase<T2>& rpyAngleRate);

template <typename T, typename T2>
templateQuat<typename T::Scalar> orienError_quat(const Eigen::MatrixBase<T>& desired_quat, const Eigen::MatrixBase<T2>& current_quat);

template <typename T, typename T2>
templateVec3<typename T::Scalar> orienError_so3(const Eigen::MatrixBase<T>& desired_quat, const Eigen::MatrixBase<T2>& current_quat);

template <typename T, typename T2, typename T3>
templateVec3<typename T::Scalar> angVelWorldSwingToGoal(const Eigen::MatrixBase<T>& desired_quat,
                                          const Eigen::MatrixBase<T2>& current_quat,
                                          const T3& time_step);

#include "orientation.inl"


#endif