//
// Created by ep1 on 11/26/20.
//

#ifndef RLML_UTILITY_H
#define RLML_UTILITY_H

#include "math.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/Geometry>

namespace RLML
{
    namespace utility
    {
        float toRadian(const float theta);
        float normAngle(const float angle);
        Eigen::Vector3f R2ypr(const Eigen::Matrix3f &R);
        Eigen::Matrix3f ypr2R(const Eigen::Vector3f &ypr);
//        template <typename Derived>
//        void reduceVector(std::vector<Derived> &v, std::vector<int> status);
        template <typename Derived>
        void reduceVector(std::vector<Derived> &v, std::vector<int> status)
        {
            int j = 0;
            for (int i = 0; i < int(v.size()); i++)
                if (status[i])
                    v[j++] = v[i];
            v.resize(j);
        }
    }
}


#endif //RLML_UTILITY_H
