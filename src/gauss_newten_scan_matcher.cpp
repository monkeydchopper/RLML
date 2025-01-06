#include "gauss_newten_scan_matcher.h"

namespace RLML
{

GaussNewtenScanMatcher::GaussNewtenScanMatcher() : max_iterations_(25)
{

}

void GaussNewtenScanMatcher::setCorrelativeGrid(float resolution, float standard_deviation)
{
    correlative_grid_ = std::shared_ptr<CorrelativeGrid>(new CorrelativeGrid(resolution, standard_deviation));
    prev_correlative_grid_ = std::shared_ptr<CorrelativeGrid>(new CorrelativeGrid(resolution, standard_deviation));

}

void GaussNewtenScanMatcher::bilinearInterpolation(const Eigen::Vector2f& coords, std::shared_ptr<CorrelativeGrid> map,
                                                   float& value, Eigen::Vector2f& derivative)
{
    if(map->isOutOfMap(coords)) {
        value = 0.0;
        derivative = Eigen::Vector2f::Zero();
        return;
    }

    Eigen::Vector2i coords_floor(coords.cast<int>());
    Eigen::Vector2f factors(coords - coords_floor.cast<float>());
    Eigen::Vector2f factors_inv(1.0f - factors[0], 1.0f - factors[1]);

    int size_x = map->getSizeX();
    int index = coords_floor[1] * size_x + coords_floor[0];
    float intensities[4];

    intensities[0] = static_cast<float>(map->getGridValue(index)) / 255.0f;
    intensities[1] = static_cast<float>(map->getGridValue(index + 1)) / 255.0f;
    intensities[2] = static_cast<float>(map->getGridValue(index + size_x)) / 255.0f;
    intensities[3] = static_cast<float>(map->getGridValue(index + size_x + 1)) / 255.0f;

    float dx1 = intensities[0] - intensities[1];
    float dx2 = intensities[2] - intensities[3];

    float dy1 = intensities[0] - intensities[2];
    float dy2 = intensities[1] - intensities[3];

    value = ((intensities[0] * factors_inv[0] + intensities[1] * factors[0]) * factors_inv[1]) +
            ((intensities[2] * factors_inv[0] + intensities[3] * factors[0]) * (factors[1]));

//    if(intensities[0] > 0.99) value = intensities[0];
//
    derivative = Eigen::Vector2f(-((dx1 * factors_inv[0]) + (dx2 * factors[0])),
                                 -((dy1 * factors_inv[1]) + (dy2 * factors[1])));
//    derivative = Eigen::Vector2f(-((dx1 * factors[0]) + (dx2 * factors[0])),
//                                 -((dy1 * factors[1]) + (dy2 * factors[1])));

//    std::cout << coords_floor.x() << " " << coords_floor.y() << " " << intensities[0] << " "
//              << intensities[1] << " " << intensities[2] << " " << intensities[3] << std::endl;
}



Eigen::Vector3f GaussNewtenScanMatcher::matchScanToMap(const Eigen::Vector3f& prior_pose,
                               std::shared_ptr<CorrelativeGrid> map,
                               std::shared_ptr<CorrelativeGrid> running_map,
                               const PointCloud& scan_data,
                               const std::vector<Eigen::Vector2f>& reflective_markers,
                               const std::map<int, std::shared_ptr<LandMark>>& associations,
                               const Eigen::Matrix3d& sigma_odom,
                               const Eigen::Matrix2d& sigma_reflector,
                               const double& sigma_scan,
                               Eigen::Matrix3d& sigma_state,
                               double& final_cost)
{

    Eigen::Vector3f pose(map->getMapPose(prior_pose));
    Eigen::Vector3f running_pose(running_map->getMapPose(prior_pose));
    Eigen::Vector3f pose_delta = running_pose - pose;
//    std::cout << "prior pose " << std::endl << prior_pose << std::endl;
//    std::cout << "pose " << std::endl << pose << std::endl;
//    std::cout << "running_pose " << std::endl << running_pose << std::endl;
    Eigen::Vector3f last_pose = pose;
    Eigen::Vector3f constraint_pose = pose;

    float resolution = map->getResolution();

//    std::cout << "original sigma" << std::endl << sigma_odom << std::endl;
    Eigen::Matrix3d scale_mat;
    scale_mat <<
              1./resolution, 0, 0,
            0, 1./resolution,0,
            0,0,1;
    Eigen::Matrix3d sigma_odom_in_map = scale_mat * sigma_odom * scale_mat.transpose();
    Eigen::Matrix3d odom_info_mat;
    odom_info_mat = sigma_odom_in_map.inverse();
    Eigen::Matrix2d reflective_info_mat = sigma_reflector.inverse() * resolution * resolution;

    double scan_info = 1. / sigma_scan;


    double cost = 0, last_cost = 0;
    float sum_value = 0;
    std::stringstream cost_ss;

    for(int i = 0; i < max_iterations_; ++i) {
        Eigen::Affine2f transform(Eigen::Translation2f(pose[0], pose[1]) * Eigen::Rotation2Df(pose[2]));

        float sin_rot = std::sin(pose[2]);
        float cos_rot = std::cos(pose[2]);

        Eigen::Matrix3f H = Eigen::Matrix3f::Zero();
        Eigen::Vector3f JTr = Eigen::Vector3f::Zero();
        cost = 0;

        float value;
        sum_value = 0;
        Eigen::Vector2f derivative;


        for(size_t j = 0; j < scan_data.size(); j = j + 1) {
            Eigen::Vector2f scan_point = scan_data[j];
            Eigen::Vector2f point = scan_point / resolution;
            Eigen::Vector2f transform_point = transform * point;
            double adaptive_scan_info = scan_info;


            bilinearInterpolation(transform_point, map, value, derivative);
//            std::cout << scan_point.x() << " " << scan_point.y() << " " << transform_point.x() << " "
//                      << transform_point.y() << " " << value << " " << derivative.x() << " " <<
//                      derivative.y() << std::endl;

            sum_value += value;
            float r = 1.0f - value;
            cost += r * r * adaptive_scan_info;
            JTr[0] += derivative[0] * r * adaptive_scan_info;
            JTr[1] += derivative[1] * r * adaptive_scan_info;

            float rot_deriv = ((-sin_rot * point.x() - cos_rot * point.y()) * derivative[0] +
                               (cos_rot * point.x() - sin_rot * point.y()) * derivative[1]);

            JTr[2] += rot_deriv * r * adaptive_scan_info;

            H(0, 0) += derivative[0] * derivative[0] * adaptive_scan_info;
            H(1, 1) += derivative[1] * derivative[1] * adaptive_scan_info;
            H(2, 2) += rot_deriv * rot_deriv * adaptive_scan_info;

            H(0, 1) += derivative[0] * derivative[1] * adaptive_scan_info;
            H(0, 2) += derivative[0] * rot_deriv * adaptive_scan_info;
            H(1, 2) += derivative[1] * rot_deriv * adaptive_scan_info;

//            std::cout << "point " << std::endl << point << std::endl;
//            std::cout << "transform point " << std::endl << transform * point << std::endl;
//            std::cout << "derivative " << std::endl << derivative << std::endl;
//            std::cout << "cost " << std::endl << cost << std::endl;
//            std::cout << "JTr " << std::endl << JTr << std::endl;
//            std::cout << "H " << std::endl << H << std::endl;

        }
        cost_ss << "old map cost "  << cost << " ";



//        Eigen::Affine2f running_transform = Eigen::Translation2f(running_pose[0], running_pose[1]) * Eigen::Rotation2Df(running_pose[2]);
//
//
//
//        for(int j = 0; j < scan_data.size(); j = j + 1) {
//            Eigen::Vector2f scan_point = scan_data[j];
//            Eigen::Vector2f point = scan_point / resolution;
//            bilinearInterpolation(running_transform * point, running_map, value, derivative);
//
//            float r = 1.0f - value;
//            cost += r * r * scan_info;
//            JTr[0] += derivative[0] * r * scan_info;
//            JTr[1] += derivative[1] * r * scan_info;
//
//            float rot_deriv = ((-sin_rot * point.x() - cos_rot * point.y()) * derivative[0] +
//                               (cos_rot * point.x() - sin_rot * point.y()) * derivative[1]);
//
//            JTr[2] += rot_deriv * r * scan_info;
//
//            H(0, 0) += derivative[0] * derivative[0] * scan_info;
//            H(1, 1) += derivative[1] * derivative[1] * scan_info;
//            H(2, 2) += rot_deriv * rot_deriv * scan_info;
//
//            H(0, 1) += derivative[0] * derivative[1] * scan_info;
//            H(0, 2) += derivative[0] * rot_deriv * scan_info;
//            H(1, 2) += derivative[1] * rot_deriv * scan_info;
//
////            std::cout << "point " << std::endl << point << std::endl;
////            std::cout << "transform point " << std::endl << transform * point << std::endl;
////            std::cout << "derivative " << std::endl << derivative << std::endl;
////            std::cout << "cost " << std::endl << cost << std::endl;
////            std::cout << "JTr " << std::endl << JTr << std::endl;
////            std::cout << "H " << std::endl << H << std::endl;
//
//        }
//        cost_ss << "running map cost   " << cost << "   ";
        H(1, 0) = H(0, 1);
        H(2, 0) = H(0, 2);
        H(2, 1) = H(1, 2);


//         add reflector
        for(auto const& asso:associations)
        {
            Eigen::Vector2f point = reflective_markers[asso.first] / resolution;
            Eigen::Vector2f marker_pos = map->getMapCoords(asso.second->getMapPos());

//            std::cout << "point " << std::endl << point << std::endl;
//            std::cout << "transform point " << std::endl << transform * point << std::endl;
//            std::cout << "marker_pos " << std::endl << marker_pos << std::endl;

            Eigen::Vector2f err = marker_pos - transform * point;

            //test
//            if(std::sqrt(point.x() * point.x() + point.y() * point.y()) < 200)
//            {
//                reflective_info_mat  = reflective_info_mat * 50;
//            }
            //test

            cost += err.transpose() * reflective_info_mat.cast<float>() * err;

//            std::cout << "err " << std::endl << err << std::endl;
//            std::cout << "cost " << std::endl << cost << std::endl;

            Eigen::Matrix<float, 2, 3> jacobian;
            jacobian <<
                     -1, 0, sin_rot * point.x() + cos_rot * point.y(),
                    0, -1, -cos_rot * point.x() + sin_rot * point.y();

            JTr += - jacobian.transpose() * reflective_info_mat.cast<float>()  * err;

            H += jacobian.transpose() * reflective_info_mat.cast<float>()  * jacobian;

//            std::cout << "jacobian " << std::endl << jacobian << std::endl;
//            std::cout << "JTr " << std::endl << JTr << std::endl;
//            std::cout << "H " << std::endl << H << std::endl;
        }
        cost_ss << "add reflective cost " << cost << " ";


        // add odom factor

        Eigen::Vector3f err = constraint_pose - pose;
        cost += err.transpose() * odom_info_mat.cast<float>()  * err;
        cost_ss << "add odom cost " << cost << " ";

        Eigen::Matrix<float, 3, 3> jacobian;
        jacobian <<
                 -1, 0, 0,
                0, -1, 0,
                0, 0, -1;
        JTr += - jacobian.transpose() * odom_info_mat.cast<float>()  * err;

        H += jacobian.transpose() * odom_info_mat.cast<float>()  * jacobian;



        if(i > 0 && cost > last_cost) {
            pose = last_pose;
            cost_ss << " cost becomes larger" << std::endl;
            break;
        }

        if(i > 0 && (std::abs(last_cost - cost) / cost) < 0.001f ) {
            cost_ss << " cost change small" << std::endl;
            break;

        }

        last_cost = cost;



        if(H.determinant() != 0.0f) {
            last_pose = pose;
            Eigen::Vector3f delta = H.inverse() * JTr;
            pose += delta;
            running_pose = pose + pose_delta;
//            std::cout << "pose " << std::endl << pose << std::endl;
//            std::cout << "running_pose " << std::endl << running_pose << std::endl;
        }
    }


    std::cout << cost_ss.str();

    final_cost = sum_value;
    pose[2] = normalizeAngle(pose[2]);

    // calculate covariance
    Eigen::Matrix3f cumulative_info = Eigen::Matrix3f::Zero();
    Eigen::Affine2f transform(Eigen::Translation2f(pose[0], pose[1]) * Eigen::Rotation2Df(pose[2]));

    float sin_rot = std::sin(pose[2]);
    float cos_rot = std::cos(pose[2]);
    float value;
    Eigen::Vector3f derivative_scan = Eigen::Vector3f::Zero();
    Eigen::Vector2f derivative = Eigen::Vector2f::Zero();
    for(size_t j = 0; j < scan_data.size(); j = j + 2) {
        Eigen::Vector2f scan_point = scan_data[j];
        Eigen::Vector2f point = scan_point / resolution;

        bilinearInterpolation(transform * point, map, value, derivative);

        derivative_scan.x() = derivative[0];
        derivative_scan.y() = derivative[1];
        float rot_deriv = ((-sin_rot * point.x() - cos_rot * point.y()) * derivative[0] +
                           (cos_rot * point.x() - sin_rot * point.y()) * derivative[1]);
        derivative_scan.z() = rot_deriv;

        cumulative_info += derivative_scan * scan_info * derivative_scan.transpose();
    }

    for(auto const& asso:associations)
    {
        Eigen::Vector2f point = reflective_markers[asso.first] / resolution;
        Eigen::Matrix<float, 2, 3> derivative_reflector;
        derivative_reflector <<
                             -1, 0, sin_rot * point.x() + cos_rot * point.y(),
                0, -1, -cos_rot * point.x() + sin_rot * point.y();

        cumulative_info += derivative_reflector.transpose() * reflective_info_mat.cast<float>()  * derivative_reflector;

    }


    Eigen::Matrix<float, 3, 3> derivative_odom;
    derivative_odom <<
                    -1, 0, 0,
            0, -1, 0,
            0, 0, -1;
    cumulative_info += derivative_odom.transpose() * odom_info_mat.cast<float>()  * derivative_odom;


    Eigen::Matrix3d cumulative_info_inv = cumulative_info.inverse().cast<double>();
    Eigen::Matrix3d scale_mat_inv = scale_mat.inverse();
    sigma_state = scale_mat_inv * cumulative_info_inv * scale_mat_inv.transpose();

//    std::cout << "updated sigma" << std::endl << sigma_state << std::endl;
    return map->getWorldPose(pose);
}

Eigen::Vector3f GaussNewtenScanMatcher::matchScanToMap(const Eigen::Vector3f& prior_pose,
                                                        std::shared_ptr<CorrelativeGrid> map,
                                                        const PointCloud& scan_data,
                                                       const std::vector<Eigen::Vector2f>& reflective_markers,
                                                       const std::map<int, std::shared_ptr<LandMark>>& associations,
                                                       const Eigen::Matrix3d& sigma_odom,
                                                       const Eigen::Matrix2d& sigma_reflector,
                                                       const double& sigma_scan,
                                                       Eigen::Matrix3d& sigma_state)
{
    Eigen::Vector3f pose(map->getMapPose(prior_pose));
//    std::cout << "prior pose " << std::endl << prior_pose << std::endl;
//    std::cout << "pose " << std::endl << pose << std::endl;
    Eigen::Vector3f last_pose = pose;
    Eigen::Vector3f constraint_pose = pose;

    float resolution = map->getResolution();

//    std::cout << "original sigma" << std::endl << sigma_odom << std::endl;
    Eigen::Matrix3d scale_mat;
    scale_mat <<
    1./resolution, 0, 0,
    0, 1./resolution,0,
    0,0,1;
    Eigen::Matrix3d sigma_odom_in_map = scale_mat * sigma_odom * scale_mat.transpose();
    Eigen::Matrix3d odom_info_mat;
    odom_info_mat = sigma_odom_in_map.inverse();

    Eigen::Matrix2d reflective_info_mat = sigma_reflector.inverse() * resolution * resolution;
    double scan_info = 1. / sigma_scan;


    double cost = 0, last_cost = 0;

    for(int i = 0; i < max_iterations_; ++i) {
        Eigen::Affine2f transform(Eigen::Translation2f(pose[0], pose[1]) * Eigen::Rotation2Df(pose[2]));

        float sin_rot = std::sin(pose[2]);
        float cos_rot = std::cos(pose[2]);

        Eigen::Matrix3f H = Eigen::Matrix3f::Zero();
        Eigen::Vector3f JTr = Eigen::Vector3f::Zero();
        cost = 0;

        float value;
        Eigen::Vector2f derivative;
        for(const Eigen::Vector2f& scan_point : scan_data) {
            Eigen::Vector2f point = scan_point / resolution;
            bilinearInterpolation(transform * point, map, value, derivative);

            float r = 1.0f - value;
            cost += r * r * scan_info;
            JTr[0] += derivative[0] * r * scan_info;
            JTr[1] += derivative[1] * r * scan_info;

            float rot_deriv = ((-sin_rot * point.x() - cos_rot * point.y()) * derivative[0] +
                  (cos_rot * point.x() - sin_rot * point.y()) * derivative[1]);

            JTr[2] += rot_deriv * r * scan_info;

            H(0, 0) += derivative[0] * derivative[0] * scan_info;
            H(1, 1) += derivative[1] * derivative[1] * scan_info;
            H(2, 2) += rot_deriv * rot_deriv * scan_info;

            H(0, 1) += derivative[0] * derivative[1] * scan_info;
            H(0, 2) += derivative[0] * rot_deriv * scan_info;
            H(1, 2) += derivative[1] * rot_deriv * scan_info;

//            std::cout << "point " << std::endl << point << std::endl;
//            std::cout << "transform point " << std::endl << transform * point << std::endl;
//            std::cout << "derivative " << std::endl << derivative << std::endl;
//            std::cout << "cost " << std::endl << cost << std::endl;
//            std::cout << "JTr " << std::endl << JTr << std::endl;
//            std::cout << "H " << std::endl << H << std::endl;

        }
//        std::cout << "map cost " << std::endl << cost << std::endl;
        H(1, 0) = H(0, 1);
        H(2, 0) = H(0, 2);
        H(2, 1) = H(1, 2);


        for(auto const& asso:associations)
        {
            Eigen::Vector2f point = reflective_markers[asso.first] / resolution;
            Eigen::Vector2f marker_pos = map->getMapCoords(asso.second->getMapPos());

//            std::cout << "point " << std::endl << point << std::endl;
//            std::cout << "transform point " << std::endl << transform * point << std::endl;
//            std::cout << "marker_pos " << std::endl << marker_pos << std::endl;

            Eigen::Vector2f err = marker_pos - transform * point;

            cost += err.transpose() * reflective_info_mat.cast<float>() * err;

//            std::cout << "err " << std::endl << err << std::endl;
//            std::cout << "cost " << std::endl << cost << std::endl;

            Eigen::Matrix<float, 2, 3> jacobian;
            jacobian <<
                     -1, 0, sin_rot * point.x() + cos_rot * point.y(),
                    0, -1, -cos_rot * point.x() + sin_rot * point.y();

            JTr += - jacobian.transpose() * reflective_info_mat.cast<float>()  * err;

            H += jacobian.transpose() * reflective_info_mat.cast<float>()  * jacobian;

//            std::cout << "jacobian " << std::endl << jacobian << std::endl;
//            std::cout << "JTr " << std::endl << JTr << std::endl;
//            std::cout << "H " << std::endl << H << std::endl;
        }
//        std::cout << "add reflective cost " << std::endl << cost << std::endl;




        // add odom factor

        Eigen::Vector3f err = constraint_pose - pose;
        cost += err.transpose() * odom_info_mat.cast<float>()  * err;
//        std::cout << "add odom cost " << std::endl << cost << std::endl;

        Eigen::Matrix<float, 3, 3> jacobian;
        jacobian <<
        -1, 0, 0,
        0, -1, 0,
        0, 0, -1;
        JTr += - jacobian.transpose() * odom_info_mat.cast<float>()  * err;

        H += jacobian.transpose() * odom_info_mat.cast<float>()  * jacobian;



        if(i > 0 && cost > last_cost) {
            pose = last_pose;
            break;
        }

        if(i > 0 && (std::abs(last_cost - cost) / cost) < 0.001f ) {
             break;
        }

        last_cost = cost;



        if(H.determinant() != 0.0f) {
            last_pose = pose;
            Eigen::Vector3f delta = H.inverse() * JTr;
            pose += delta;

        }
    }
    pose[2] = normalizeAngle(pose[2]);

    // calculate covariance
    Eigen::Matrix3f cumulative_info = Eigen::Matrix3f::Zero();
    Eigen::Affine2f transform(Eigen::Translation2f(pose[0], pose[1]) * Eigen::Rotation2Df(pose[2]));

    float sin_rot = std::sin(pose[2]);
    float cos_rot = std::cos(pose[2]);
    float value;
    Eigen::Vector3f derivative_scan = Eigen::Vector3f::Zero();
    Eigen::Vector2f derivative = Eigen::Vector2f::Zero();
    for(const Eigen::Vector2f& scan_point : scan_data) {
        Eigen::Vector2f point = scan_point / resolution;

        bilinearInterpolation(transform * point, map, value, derivative);

        derivative_scan.x() = derivative[0];
        derivative_scan.y() = derivative[1];
        float rot_deriv = ((-sin_rot * point.x() - cos_rot * point.y()) * derivative[0] +
                           (cos_rot * point.x() - sin_rot * point.y()) * derivative[1]);
        derivative_scan.z() = rot_deriv;

        cumulative_info += derivative_scan * scan_info * derivative_scan.transpose();
    }

    for(auto const& asso:associations)
    {
        Eigen::Vector2f point = reflective_markers[asso.first] / resolution;
        Eigen::Matrix<float, 2, 3> derivative_reflector;
        derivative_reflector <<
                 -1, 0, sin_rot * point.x() + cos_rot * point.y(),
                0, -1, -cos_rot * point.x() + sin_rot * point.y();

        cumulative_info += derivative_reflector.transpose() * reflective_info_mat.cast<float>()  * derivative_reflector;
    }


    Eigen::Matrix<float, 3, 3> derivative_odom;
    derivative_odom <<
             -1, 0, 0,
            0, -1, 0,
            0, 0, -1;
    cumulative_info += derivative_odom.transpose() * odom_info_mat.cast<float>()  * derivative_odom;


    Eigen::Matrix3d cumulative_info_inv = cumulative_info.inverse().cast<double>();
    Eigen::Matrix3d scale_mat_inv = scale_mat.inverse();
    sigma_state = scale_mat_inv * cumulative_info_inv * scale_mat_inv.transpose();

//    std::cout << "updated sigma" << std::endl << sigma_state << std::endl;
    return map->getWorldPose(pose);
}


Eigen::Vector3f GaussNewtenScanMatcher::matchScanToMap(const Eigen::Vector3f& prior_pose, std::shared_ptr<CorrelativeGrid> map,
                               const PointCloud& scan_data)
{
    Eigen::Vector3f pose(map->getMapPose(prior_pose));
    Eigen::Vector3f last_pose = pose;

    double cost = 0, last_cost = 0;

    for(int i = 0; i < max_iterations_; ++i) {
        Eigen::Affine2f transform(Eigen::Translation2f(pose[0], pose[1]) * Eigen::Rotation2Df(pose[2]));

        float sin_rot = sin(pose[2]);
        float cos_rot = cos(pose[2]);

        Eigen::Matrix3f H = Eigen::Matrix3f::Zero();
        Eigen::Vector3f JTr = Eigen::Vector3f::Zero();
        cost = 0;

        float value;
        float resolution = map->getResolution();
        Eigen::Vector2f derivative;
        for(const Eigen::Vector2f& scan_point : scan_data) {
            Eigen::Vector2f point = scan_point / resolution;
            bilinearInterpolation(transform * point, map, value, derivative);

            float r = 1.0f - value;
            cost += r * r;
            JTr[0] += derivative[0] * r;
            JTr[1] += derivative[1] * r;

            float rot_deriv = ((-sin_rot * point.x() - cos_rot * point.y()) * derivative[0] +
                               (cos_rot * point.x() - sin_rot * point.y()) * derivative[1]);

            JTr[2] += rot_deriv * r;

            H(0, 0) += derivative[0] * derivative[0];
            H(1, 1) += derivative[1] * derivative[1];
            H(2, 2) += rot_deriv * rot_deriv;

            H(0, 1) += derivative[0] * derivative[1];
            H(0, 2) += derivative[0] * rot_deriv;
            H(1, 2) += derivative[1] * rot_deriv;

//            std::cout << "point " << std::endl << point << std::endl;
//            std::cout << "transform point " << std::endl << transform * point << std::endl;
//            std::cout << "derivative " << std::endl << derivative << std::endl;
//            std::cout << "cost " << std::endl << cost << std::endl;
//            std::cout << "JTr " << std::endl << JTr << std::endl;
//            std::cout << "H " << std::endl << H << std::endl;

        }

        H(1, 0) = H(0, 1);
        H(2, 0) = H(0, 2);
        H(2, 1) = H(1, 2);


        if(i > 0 && cost > last_cost) {
            pose = last_pose;
            break;
        }
        last_cost = cost;



        if(H.determinant() != 0.0f) {
            last_pose = pose;
            Eigen::Vector3f delta = H.inverse() * JTr;
            pose += delta;
            pose[2] = normalizeAngle(pose[2]);
        }
    }

    return map->getWorldPose(pose);
}


Eigen::Vector3f GaussNewtenScanMatcher::matchScan(const std::shared_ptr<LaserScan>& scan,
                          const std::vector<std::shared_ptr<LaserScan>>& base_scan,
                          const std::vector<std::shared_ptr<LaserScan>>& prev_scan,
                          const std::vector<Eigen::Vector2f>& reflective_markers,
                          const std::map<int, std::shared_ptr<LandMark>>& associations,
                          const Eigen::Matrix3d& sigma_odom,
                          const Eigen::Matrix2d& sigma_reflector,
                          const double& sigma_scan,
                          Eigen::Matrix3d& sigma_state,
                          double& final_cost)
{
    correlative_grid_->addScans(base_scan);
    prev_correlative_grid_->addScans(prev_scan);
    return matchScanToMap( scan->getPose(), prev_correlative_grid_, correlative_grid_, scan->getRawPointCloud(),
            reflective_markers, associations, sigma_odom, sigma_reflector, sigma_scan, sigma_state, final_cost);
}



Eigen::Vector3f GaussNewtenScanMatcher::matchScan(const std::shared_ptr<LaserScan>& scan,
                                                  const std::vector<std::shared_ptr<LaserScan>>& base_scan,
                                                  const std::vector<Eigen::Vector2f>& reflective_markers,
                                                  const std::map<int, std::shared_ptr<LandMark>>& associations,
                                                  const Eigen::Matrix3d& sigma_odom,
                                                  const Eigen::Matrix2d& sigma_reflector,
                                                  const double& sigma_scan,
                                                  Eigen::Matrix3d& sigma_state)
{
    correlative_grid_->addScans(base_scan);

    return matchScanToMap( scan->getPose(), correlative_grid_, scan->getRawPointCloud(), reflective_markers, associations,
                          sigma_odom, sigma_reflector, sigma_scan, sigma_state);
}


Eigen::Vector3f GaussNewtenScanMatcher::matchScan(const std::shared_ptr<LaserScan>& scan,
                          const std::vector<std::shared_ptr<LaserScan>>& base_scan)
{
    correlative_grid_->addScans(base_scan);

    return matchScanToMap(scan->getPose(), correlative_grid_, scan->getRawPointCloud());

}

} // namespace RLML
