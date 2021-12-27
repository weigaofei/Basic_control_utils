#include "orientation.hpp"

template <typename T>
templateMat3<typename T::Scalar> crossMatrix(const Eigen::MatrixBase<T>& v) {
  static_assert(T::ColsAtCompileTime == 1 && T::RowsAtCompileTime == 3,
                "Must have 3x1 matrix");
  templateMat3<typename T::Scalar> m;
  m << 0, -v[2], v[1], 
      v[2], 0, -v[0], 
      -v[1], v[0], 0;
  return m;
}

template <typename T>
templateMat3<T> coordinateRotation(CoordinateAxis axis, const T& theta){
  T s = std::sin(theta);
  T c = std::cos(theta);
  templateMat3<T> R;
  if(axis == CoordinateAxis::X){
    R << 1, 0, 0, 0, c, -s, 0, s, c;
  }
  else if(axis == CoordinateAxis::Y){
    R << c, 0, s, 0, 1, 0, -s, 0, c;
  }
  else if(axis == CoordinateAxis::Z){
    R << c, -s, 0, s, c, 0, 0, 0, 1;
  }
  
  return R;
}

template <typename T>
templateVec3<typename T::Scalar> quatToRPY(const Eigen::MatrixBase<T>& q){
  static_assert(T::ColsAtCompileTime == 1 && T::RowsAtCompileTime == 4,
                "Must have 4x1 quat");
  templateVec3<typename T::Scalar> rpy;
  typename T::Scalar as = std::min(-2. * (q[1] * q[3] - q[0] * q[2]), .99999);
  rpy(2) = std::atan2(2 * (q[1] * q[2] + q[0] * q[3]),
                 square(q[0]) + square(q[1]) - square(q[2]) - square(q[3]));
  rpy(1) = std::asin(as);
  rpy(0) = std::atan2(2 * (q[2] * q[3] + q[0] * q[1]),
                 square(q[0]) - square(q[1]) - square(q[2]) + square(q[3]));
  return rpy;
}

template <typename T>
templateQuat<typename T::Scalar> rpyToQuat(const Eigen::MatrixBase<T>& rpyAngle){
  static_assert(T::ColsAtCompileTime == 1 && T::RowsAtCompileTime == 3,
                "Must have 3x1 vec");
  templateMat3<typename T::Scalar> body2world = templateMat3<typename T::Scalar>::Zero();
  body2world = bodyToWorldMatrix(rpyAngle);
  templateQuat<typename T::Scalar> quat(body2world);

  return quat;
}

template <typename T>
templateMat3<typename T::Scalar> quatToRotationMatrix(const Eigen::MatrixBase<T>& q){
  static_assert(T::ColsAtCompileTime == 1 && T::RowsAtCompileTime == 4, "Must have 4x1 quat");
  typename T::Scalar e0 = q(0);
  typename T::Scalar e1 = q(1);
  typename T::Scalar e2 = q(2);
  typename T::Scalar e3 = q(3);

  templateMat3<typename T::Scalar> R;

  R << 1 - 2 * (e2 * e2 + e3 * e3), 2 * (e1 * e2 + e0 * e3), 2 * (e1 * e3 - e0 * e2),
        2 * (e1 * e2 - e0 * e3), 1 - 2 * (e1 * e1 + e3 * e3), 2 * (e2 * e3 + e0 * e1),
        2 * (e1 * e3 + e0 * e2),  2 * (e2 * e3 - e0 * e1), 1 - 2 * (e1 * e1 + e2 * e2);

  return R;
}

template <typename T>
templateVec3<typename T::Scalar> rotationMatrixToRPY(const Eigen::MatrixBase<T>& R){
  templateQuat<typename T::Scalar> quat(R);
  templateVec4<typename T::Scalar> quaVec;
  quaVec << quat.w(), quat.vec();

  return quatToRPY(quaVec);
}

template <typename T>
templateMat3<typename T::Scalar> bodyToWorldMatrix(const Eigen::MatrixBase<T>& rpyAngle){
  static_assert(T::ColsAtCompileTime == 1 && T::RowsAtCompileTime == 3,
                "must have 3x1 vector");
  templateMat3<typename T::Scalar> m = coordinateRotation(CoordinateAxis::Z, rpyAngle[2]) *
                               coordinateRotation(CoordinateAxis::Y, rpyAngle[1]) *
                               coordinateRotation(CoordinateAxis::X, rpyAngle[0]);
  return m;
}

template <typename T, typename T2>
templateVec3<typename T::Scalar> eulerAngleRateToAngVelWorld(const Eigen::MatrixBase<T>& rpyAngle,const Eigen::MatrixBase<T2>& rpyAngleRate){
  static_assert(T::ColsAtCompileTime == 1 && T::RowsAtCompileTime == 3,
                "Must have 3x1 rpyAngle");
  static_assert(T2::ColsAtCompileTime == 1 && T2::RowsAtCompileTime == 3,
                "Must have 3x1 rpyAngleRate");
  templateMat3<typename T::Scalar> R;
  templateVec3<typename T::Scalar> angvel;
  T spitch = std::sin(rpyAngle[1]);
  T cpitch = std::cos(rpyAngle[1]);
  T syaw = std::sin(rpyAngle[2]);
  T cyaw = std::cos(rpyAngle[2]);

  R << cpitch * cyaw, -syaw, 0, 
      cpitch * syaw, cyaw, 0,
      -spitch, 0, 1;
  angvel = R * rpyAngleRate;

  return angvel;
}

template <typename T>
templateVec3<typename T::Scalar> quaternionToso3(const Eigen::MatrixBase<T>& q){
  static_assert(T::ColsAtCompileTime == 1 && T::RowsAtCompileTime == 4,
                "Must have 4x1 quat");
  templateQuat<typename T::Scalar> normalized_quaternion(q[0], q[1], q[2], q[3]);
  normalized_quaternion.normalize();
  if (normalized_quaternion.w() < 0.0) {
    normalized_quaternion.w() *= -1.0;
    normalized_quaternion.vec() *= -1.0;
  }
  if (std::abs(1.0-normalized_quaternion.w()) < 0.01) {
    return 2.0*normalized_quaternion.vec();
  } else {
    typename T::Scalar theta = 2.0*std::acos(normalized_quaternion.w());
    return theta/std::sin(0.5*theta) * normalized_quaternion.vec();
  }
}

template <typename T>
templateVec3<typename T::Scalar> so3ToQuaternion(const Eigen::MatrixBase<T>& rotation_vector){
    static_assert(T::ColsAtCompileTime == 1 && T::RowsAtCompileTime == 3,
                "Must have 3x1 rotation_vector");
    templateQuat<typename T::Scalar> quaternion;
    if (rotation_vector.norm() < 0.01) {
      quaternion = templateQuat<typename T::Scalar>(1.0, 0.5*rotation_vector.x(), 0.5*rotation_vector.y(), 0.5*rotation_vector.z());
    } else {
      typename T::Scalar norm = std::cos(rotation_vector.norm()/2.0);
      quaternion.w() = norm;
      quaternion.vec() = rotation_vector/rotation_vector.norm() * std::sin(0.5*rotation_vector.norm());
    }
    quaternion.normalize();
    if (quaternion.w() < 0.0) {
      quaternion.w() *= -1.0;
      quaternion.vec() *= -1.0;
    }
    return quaternion;
}

template <typename T>
templateMat4<typename T::Scalar> quatToProductMatrix(const Eigen::MatrixBase<T>& quat, bool isRightProduct){
    static_assert(T::ColsAtCompileTime == 1 && T::RowsAtCompileTime == 4,
                  "Must have 4x1 quat");
    templateMat4<typename T::Scalar> q_matrix = quat.w() * templateMat4<typename T::Scalar>::Identity();
    q_matrix.template bottomLeftCorner<3,1>() = quat.vec();
    q_matrix.template topRightCorner<1,3>() = - quat.vec();
    if (!isRightProduct) { q_matrix.template bottomRightCorner<3,3>() += crossMatrix(quat.vec()); }
    else               { q_matrix.template bottomRightCorner<3,3>() -= crossMatrix(quat.vec()); }

    return q_matrix;
}

template <typename T>
templateMat4c3<typename T::Scalar> jacobianQuatRateWrtEulerAngleRate(const Eigen::MatrixBase<T>& rotation_vector){
  static_assert(T::ColsAtCompileTime == 1 && T::RowsAtCompileTime == 3,
                  "Must have 3x1 rotation_vector");
  templateMat4c3<typename T::Scalar> jacobian = templateMat4c3<typename T::Scalar>::Zero();
  if (rotation_vector.norm() < 0.01) {
    jacobian.template topLeftCorner<1,3>() = -0.25*rotation_vector;
    jacobian.template bottomLeftCorner<3,3>() = 0.5*templateMat3<typename T::Scalar>::Identity();
  } 
  else {
    typename T::Scalar factor1 = 0.5*std::sin(0.5*rotation_vector.norm())/rotation_vector.norm();
    typename T::Scalar factor2 = ( std::cos(0.5*rotation_vector.norm())*rotation_vector.norm() - 2.0*std::sin(0.5*rotation_vector.norm()) )
                      / ( 2.0*std::pow(rotation_vector.norm(), 3.0) );
    jacobian.template topLeftCorner<1,3>() = -factor1*rotation_vector.transpose();
    jacobian.template bottomLeftCorner<3,3>() = 2.0*factor1*templateMat3<typename T::Scalar>::Identity() + factor2*(rotation_vector*rotation_vector.transpose());
  }

  return jacobian;
}

template <typename T>
templateMat3c4<typename T::Scalar> jacobianAngVelWorldWrtQuatRate(const Eigen::MatrixBase<T>& quat){
  static_assert(T::ColsAtCompileTime == 1 && T::RowsAtCompileTime == 4,
                  "Must have 4x1 quatRate");
  templateQuat<typename T::Scalar> quaternion(quat[0], quat[1], quat[2], quat[3]);
  templateMat3c4<typename T::Scalar> jacobianAngularVelocityWrtQuaternion = quatToProductMatrix(quaternion, false).template bottomRows<3>();
  jacobianAngularVelocityWrtQuaternion.template bottomLeftCorner<3,1>() *= -1.0;

  jacobianAngularVelocityWrtQuaternion *= 2;

	return jacobianAngularVelocityWrtQuaternion;
}

template <typename T, typename T2>
templateVec3<typename T::Scalar> eulerAngleRateToAngVelWorldUsingQuat(const Eigen::MatrixBase<T>& rpyAngle,const Eigen::MatrixBase<T2>& rpyAngleRate){
  static_assert(T::ColsAtCompileTime == 1 && T::RowsAtCompileTime == 3,
                "Must have 3x1 rpyAngle");
  static_assert(T2::ColsAtCompileTime == 1 && T2::RowsAtCompileTime == 3,
                "Must have 3x1 rpyAngleRate");
  templateVec3<typename T::Scalar> angvel;
  angvel = jacobianAngVelWorldWrtQuatRate(rpyToQuat(rpyAngle)) * jacobianQuatRateWrtEulerAngleRate(rpyAngle) * rpyAngleRate;

  return angvel;
}

template <typename T, typename T2>
templateQuat<typename T::Scalar> orienError_quat(const Eigen::MatrixBase<T>& desired_quat, const Eigen::MatrixBase<T2>& current_quat){
  static_assert(T::ColsAtCompileTime == 1 && T::RowsAtCompileTime == 4,
                "Must have 4x1 quat");
  static_assert(T2::ColsAtCompileTime == 1 && T2::RowsAtCompileTime == 4,
                "Must have 4x1 quat");
  templateQuat<typename T::Scalar> q_desired(desired_quat[0], desired_quat[1], desired_quat[2], desired_quat[3]);
  templateQuat<typename T::Scalar> q_current(current_quat[0], current_quat[1], current_quat[2], current_quat[3]);

  if(q_current.dot(q_desired) < 0) {q_current.w() *= -1.0; q_current.vec() *= -1.0;}
  templateQuat<typename T::Scalar> deltaQuat = q_desired * q_current.inverse();

  return deltaQuat;
}

template <typename T, typename T2>
templateVec3<typename T::Scalar> orienError_so3(const Eigen::MatrixBase<T>& desired_quat, const Eigen::MatrixBase<T2>& current_quat){
  static_assert(T::ColsAtCompileTime == 1 && T::RowsAtCompileTime == 4,
                "Must have 4x1 quat");
  static_assert(T2::ColsAtCompileTime == 1 && T2::RowsAtCompileTime == 4,
                "Must have 4x1 quat");
  templateQuat<typename T::Scalar> quat_error = orienError_quat(desired_quat, current_quat);
  templateVec3<typename T::Scalar> ori_error_so3 = quaternionToso3(quat_error);

  return ori_error_so3;
}

template <typename T, typename T2, typename T3>
templateVec3<typename T::Scalar> angVelWorldSwingToGoal(const Eigen::MatrixBase<T>& desired_quat,
                                          const Eigen::MatrixBase<T2>& current_quat,
                                          const T3& time_step){
  static_assert(T::ColsAtCompileTime == 1 && T::RowsAtCompileTime == 4,
                "Must have 4x1 quat");
  static_assert(T2::ColsAtCompileTime == 1 && T2::RowsAtCompileTime == 4,
                "Must have 4x1 quat");
  templateVec3<typename T::Scalar> delta_ori_so3 = orienError_so3(desired_quat, current_quat);
  templateVec3<typename T::Scalar> current_ori_so3 = quaternionToso3(current_quat);

  templateVec3<typename T::Scalar> angvel = eulerAngleRateToAngVelWorldUsingQuat(current_ori_so3, delta_ori_so3 / T::Scalar);

  return angvel;
}