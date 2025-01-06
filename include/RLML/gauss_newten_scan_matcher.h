#ifndef GAUSS_NEWTEN_SCAN_MATCHER_H
#define GAUSS_NEWTEN_SCAN_MATCHER_H

#include <iostream>
#include <memory>
#include "laser_scan.h"
#include "correlative_grid.h"
#include "math_func.h"

namespace RLML
{

class GaussNewtenScanMatcher
{
public:
    GaussNewtenScanMatcher();
    ~GaussNewtenScanMatcher() {}

    void setCorrelativeGrid(float resolution, float standard_deviation);
    Eigen::Vector3f matchScan(const std::shared_ptr<LaserScan>& scan,
                              const std::vector<std::shared_ptr<LaserScan>>& base_scan,
                              const std::vector<std::shared_ptr<LaserScan>>& prev_scan,
                              const std::vector<Eigen::Vector2f>& reflective_markers,
                              const std::map<int, std::shared_ptr<LandMark>>& associations,
                              const Eigen::Matrix3d& sigma_odom,
                              const Eigen::Matrix2d& sigma_reflector,
                              const double& sigma_scan,
                              Eigen::Matrix3d& sigma_state,
                              double& final_cost);
    Eigen::Vector3f matchScan(const std::shared_ptr<LaserScan>& scan,
                              const std::vector<std::shared_ptr<LaserScan>>& base_scan,
                              const std::vector<Eigen::Vector2f>& reflective_markers,
                              const std::map<int, std::shared_ptr<LandMark>>& associations,
                              const Eigen::Matrix3d& sigma_odom,
                              const Eigen::Matrix2d& sigma_reflector,
                              const double& sigma_scan,
                              Eigen::Matrix3d& sigma_state);
    Eigen::Vector3f matchScan(const std::shared_ptr<LaserScan>& scan,
                              const std::vector<std::shared_ptr<LaserScan>>& base_scan);
    Eigen::Vector3f matchScanToMap(const Eigen::Vector3f& prior_pose, std::shared_ptr<CorrelativeGrid> map,
                                   const PointCloud& scan_data,
                                   const std::vector<Eigen::Vector2f>& reflective_markers,
                                   const std::map<int, std::shared_ptr<LandMark>>& associations,
                                   const Eigen::Matrix3d& sigma_odom,
                                   const Eigen::Matrix2d& sigma_reflector,
                                   const double& sigma_scan,
                                   Eigen::Matrix3d& sigma_state);
    Eigen::Vector3f matchScanToMap(const Eigen::Vector3f& prior_pose, std::shared_ptr<CorrelativeGrid> map,
                                   std::shared_ptr<CorrelativeGrid> running_map,
                                   const PointCloud& scan_data,
                                   const std::vector<Eigen::Vector2f>& reflective_markers,
                                   const std::map<int, std::shared_ptr<LandMark>>& associations,
                                   const Eigen::Matrix3d& sigma_odom,
                                   const Eigen::Matrix2d& sigma_reflector,
                                   const double& sigma_scan,
                                   Eigen::Matrix3d& sigma_state,
                                   double& final_cost);
    Eigen::Vector3f matchScanToMap(const Eigen::Vector3f& prior_pose, std::shared_ptr<CorrelativeGrid> map,
                                   const PointCloud& scan_data);
    std::shared_ptr<CorrelativeGrid> getCorrelativeGrid() {  return prev_correlative_grid_; }

private:
    void bilinearInterpolation(const Eigen::Vector2f& coords, std::shared_ptr<CorrelativeGrid> map,
                               float& value, Eigen::Vector2f& derivative);

private:
    int max_iterations_;
    std::shared_ptr<CorrelativeGrid> correlative_grid_;
    std::shared_ptr<CorrelativeGrid> prev_correlative_grid_;
};

}

#endif // GAUSS_NEWTEN_SCAN_MATCHING_H
