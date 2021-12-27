#ifndef _CPPTYPES_H_
#define _CPPTYPES_H_

#include <Eigen/Dense>
#include <stdexcept>
#include <iostream>

template <typename T>
using templateVec3 = typename Eigen::Matrix<T, 3, 1>;

template <typename T>
using templateVec4 = typename Eigen::Matrix<T, 4, 1>;

template <typename T>
using templateMat3 = typename Eigen::Matrix<T, 3, 3>;

template <typename T>
using templateMat4 = typename Eigen::Matrix<T, 4, 4>;

template <typename T>
using templateMat4c3 = typename Eigen::Matrix<T, 4, 3>;

template <typename T>
using templateMat3c4 = typename Eigen::Matrix<T, 3, 4>;

template <typename T>
using templateQuat = typename Eigen::Quaternion<T>;

template <typename T>
using templateDVec = typename Eigen::Matrix<T, Eigen::Dynamic, 1>;

template <typename T>
using templateDMat = typename Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

#endif