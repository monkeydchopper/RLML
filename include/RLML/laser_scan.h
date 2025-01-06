#ifndef LASER_SCAN_H
#define LASER_SCAN_H

#include <iostream>
#include <Eigen/Geometry>
#include <vector>
#include <memory>
#include "map"
#include "landmark.h"
#include "utility.h"

typedef std::vector<Eigen::Vector2f> PointCloud;


namespace RLML
{

    class LandMark;

class Range
{
public:
    Range() : min_(std::numeric_limits<float>::max(), std::numeric_limits<float>::max()),
              max_(std::numeric_limits<float>::min(), std::numeric_limits<float>::min()) {}
    void addPoint(const Eigen::Vector2f& point)
    {
        if(point[0] < min_[0]) {
            min_[0] = point[0];
        }
        else if(point[0] > max_[0]) {
            max_[0] = point[0];
        }

        if(point[1] < min_[1]) {
            min_[1] = point[1];
        }
        else if(point[1] > max_[1]) {
            max_[1] = point[1];
        }
    }

    void addRange(const Range& box)
    {
        addPoint(box.min_);
        addPoint(box.max_);
    }

    const Eigen::Vector2f& getMin() { return min_; }
    const Eigen::Vector2f& getMax() { return max_; }

private:
    Eigen::Vector2f min_;
    Eigen::Vector2f max_;
};

class LaserScan
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    LaserScan() {
        lidar_to_odom_ = Eigen::Vector3f(0,0,0);
    }
    LaserScan(const PointCloud& points)
    {
        raw_points_ = points;
        lidar_to_odom_ = Eigen::Vector3f(0,0,0);
    }
    ~LaserScan() {}

    const PointCloud& getRawPointCloud() { return raw_points_; }
    const PointCloud& getTransformedPointCloud() { return transformed_points_; }
    void updatePointCloud(std::vector<int> status)
    {
        utility::reduceVector(raw_points_, status);
        utility::reduceVector(transformed_points_, status);
    }

    void setCalibration(const Eigen::Vector3f& lidar_to_odom){lidar_to_odom_ = lidar_to_odom;}
    Eigen::Vector3f getCalibration(){return lidar_to_odom_;}
    void setId(int id) { id_ = id; }
    int getId() { return id_; }
    void setPose(const Eigen::Vector3f& pose) {
        pose_ = pose;
        Eigen::AngleAxisf rotation(pose_[2], Eigen::Vector3f(0, 0, 1));
        scan_pose_ = pose_ + rotation * lidar_to_odom_;
    }
    Eigen::Vector3f getPose() { return pose_; }
    Eigen::Vector3f getScanPose() {return scan_pose_;}
    Range getRange() { return range_; }

    void addLandmark(const Eigen::Vector2f& landmark_in_fm_pos, const std::shared_ptr<LandMark>& landmark_ptr)
    {
        landmarks_associations_[landmark_ptr] = landmark_in_fm_pos;
    }

    void deleteLandmark(const std::shared_ptr<LandMark>& landmark_ptr)
    {
        landmarks_associations_.erase(landmark_ptr);
    }

    std::map<std::shared_ptr<LandMark>, Eigen::Vector2f> getLandmarks()
    {
        return landmarks_associations_;
    }

    void replaceLandmark(const std::shared_ptr<LandMark>& landmark_discard, const std::shared_ptr<LandMark>& landmark_replace)
    {
        landmarks_associations_[landmark_replace] = landmarks_associations_[landmark_discard];
//        landmarks_associations_.erase(landmark_discard);
    }

    Eigen::Vector2f transformLandmark(const Eigen::Vector2f& landmark_in_fm_pos)
    {
        Eigen::Affine2f transform(Eigen::Translation2f(pose_[0], pose_[1]) * Eigen::Rotation2Df(pose_[2]));
        return transform * landmark_in_fm_pos;
    }

    void transformPointCloud()
    {
        Range range;
        transformed_points_.clear();
        Eigen::Affine2f transform(Eigen::Translation2f(pose_[0], pose_[1]) * Eigen::Rotation2Df(pose_[2]));

        for(const Eigen::Vector2f& point : raw_points_) {
            Eigen::Vector2f transformed_point = transform * point;
            transformed_points_.push_back(transformed_point);
            range.addPoint(transformed_point);
        }

        range.addPoint(pose_.head<2>());

        range_ = range;
    }




private:
    std::map<std::shared_ptr<LandMark>, Eigen::Vector2f> landmarks_associations_;
    PointCloud raw_points_;
    PointCloud transformed_points_;
    Eigen::Vector3f pose_;
    Eigen::Vector3f scan_pose_;
    Eigen::Vector3f lidar_to_odom_;
    Range range_;
    int id_;
};

} // namespace RLML

#endif // LASER_SCAN_H
