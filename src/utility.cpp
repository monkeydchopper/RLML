//
// Created by ep1 on 11/26/20.
//

#include "utility.h"


namespace RLML
{
    namespace utility
    {

        float toRadian(const float theta)
        {
            return theta * (float)M_PI / 180.f;
        }

        float normAngle(const float angle)
        {
            float norm = angle;
            if(angle > (float)M_PI) norm = (angle - 2*(float)M_PI);
            if(angle < -(float)M_PI) norm = (angle + 2*(float)M_PI);
            return norm;
        }

        Eigen::Vector3f R2ypr(const Eigen::Matrix3f &R)
        {
            Eigen::Vector3f n = R.col(0);
            Eigen::Vector3f o = R.col(1);
            Eigen::Vector3f a = R.col(2);

            Eigen::Vector3f ypr(3);

            float y = std::atan2(n(1), n(0));
            float r = std::atan2(o(2),a(2));
            float p = std::atan2(-n(2), o(2) * std::sin(r) + a(2) * std::cos(r));
            ypr(0) = y;
            ypr(1) = p;
            ypr(2) = r;

            return ypr / (float)M_PI * 180.f;
        }


        Eigen::Matrix3f ypr2R(const Eigen::Vector3f &ypr)
        {

            float y = ypr(0) / 180.f * (float)M_PI;
            float p = ypr(1) / 180.f * (float)M_PI;
            float r = ypr(2) / 180.f * (float)M_PI;

            Eigen::Matrix<float, 3, 3> Rz;
            Rz << std::cos(y), -std::sin(y), 0,
                    std::sin(y), std::cos(y), 0,
                    0, 0, 1;

            Eigen::Matrix<float, 3, 3> Ry;
            Ry << std::cos(p), 0., std::sin(p),
                    0., 1., 0.,
                    -std::sin(p), 0., std::cos(p);

            Eigen::Matrix<float, 3, 3> Rx;
            Rx << 1., 0., 0.,
                    0., std::cos(r), -std::sin(r),
                    0., std::sin(r), std::cos(r);

            return Rz * Ry * Rx;
        }



    }
}