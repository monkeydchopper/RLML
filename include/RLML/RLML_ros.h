#ifndef RLML_ROS_H
#define RLML_ROS_H

#include <ros/ros.h>
#include <message_filters/subscriber.h>
#include <message_filters/time_synchronizer.h>
#include <message_filters/sync_policies/approximate_time.h>
#include <nav_msgs/OccupancyGrid.h>
#include <nav_msgs/Path.h>
#include <visualization_msgs/MarkerArray.h>
#include <sensor_msgs/LaserScan.h>
#include <nav_msgs/Odometry.h>
#include <sensor_msgs/Imu.h>
#include <std_srvs/Empty.h>
#include <tf/transform_listener.h>
#include <tf/transform_broadcaster.h>
#include <tf_conversions/tf_eigen.h>
#include <tf/transform_datatypes.h>
#include <tf2_ros/static_transform_broadcaster.h>
#include <geometry_msgs/TransformStamped.h>

#include <thread>
#include <mutex>
#include <fstream>
#include <functional>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl_conversions/pcl_conversions.h>

#include <xmlrpcpp/XmlRpcException.h>
#include "map_builder.h"
#include "utility.h"
#include "landmark.h"


struct sycData
{
    sensor_msgs::LaserScanConstPtr scan_msg;
    nav_msgs::Odometry::ConstPtr odom_msg;
};


struct robotState
{
    float v, omega;
    Eigen::Vector3f pose;
    ros::Time stamp;
    robotState(){}
    robotState(const float v_, const float omega_, const Eigen::Vector3f& pose_, const ros::Time stamp_):
    v(v_), omega(omega_), pose(pose_), stamp(stamp_)
    {}
};

class fusionSlamRos
{
public:
    fusionSlamRos();
    ~fusionSlamRos();

private:
    void scanCallback(const sensor_msgs::LaserScanConstPtr& scan_msg, const nav_msgs::Odometry::ConstPtr& odom_msg);
    void odomCallback(const nav_msgs::Odometry::ConstPtr& odom_msg);
    void processData();
    void cmdProcess();
    void publishPoseUpToDate();
    void publishPose(const Eigen::Vector3f& pose, const ros::Time& stamp);
    void publishOccupancyGridMap();
    void publishProbabilityGridMap();
    void publishCorrelativeGrid();
    void publishLoop();
    void publishPath();
    void publishConstraintList();
    void publishCurLandmarks(const ros::Time& stamp);
    void publishAllLandmarks();
    void setCorrelativeTranslationTable();
    bool optimizationCallback(std_srvs::Empty::Request& req, std_srvs::Empty::Response& res);

    void saveMap();
    void updateMap();
    void saveOccupancyGridMap();
    template <typename Derived>
    void loadCovariance(const ros::NodeHandle& nh,const std::string& cov_name, Eigen::DenseBase<Derived>& covariance);

private:
    // frame id
    std::string base_frame_, map_frame_, lidar_frame_;
    std::string save_path_, load_path_;

    bool use_publish_pose_thread_;
    double publish_freq_, pose_freq_;

    // flags
    bool got_lidar_tf_;
    bool initialize_;

    tf::TransformListener tf_listener_;
    tf::TransformBroadcaster tf_broadcaster_;

    ros::Publisher map_pub_;
    ros::Publisher correlative_grid_pub_;
    ros::Publisher pose_pub_;
    ros::Publisher path_pub_;
    ros::Publisher constraint_list_pub_;
    ros::Publisher landmarks_pub_, cur_landmarks_pub_;
    ros::Subscriber odom_sub_;
    ros::ServiceServer optimization_srv_, configuration_srv_;

    message_filters::Subscriber<sensor_msgs::LaserScan>  lidar_subscriber_;
    message_filters::Subscriber<nav_msgs::Odometry>  odom_subscriber_;
    typedef message_filters::sync_policies::ApproximateTime<
            sensor_msgs::LaserScan,
            nav_msgs::Odometry> syncPolicy;
    typedef message_filters::Synchronizer<syncPolicy> Sync;
    std::shared_ptr<Sync> sync_;

    tf::StampedTransform lidar_tf_;
    Eigen::Affine3d lidar_trans_;
    Eigen::Vector3f lidar_to_odom_;

    robotState cur_state_;
    double previous_t_;
    bool mapping_;
    bool use_reflector_, use_landmark_area_, use_reflector_relocalize_, preset_lidar_tf_;

    std::shared_ptr<std::thread> proc_th_ptr_;
    std::deque<sycData> dataques_;


    RLML::MapBuilder map_builder_;
    std::shared_ptr<std::thread> publish_thread_;
    std::shared_ptr<std::thread> cmdThPtr_;
    std::shared_ptr<std::thread> poseThPtr_;
    std::mutex mBuf_, mPub_;

    char* correlative_translation_table_;

    double reflector_max_length_, reflector_min_length_, reflector_mean_length_;
    double intensity_min_;

    float lidar_trans_x_, lidar_trans_y_, lidar_yaw_;

    Eigen::Vector3f initial_pose_;

    bool relocalize_, tracking_succeed_;


    bool receive_first_lidar_, receive_first_odom_;
    double publish_time_, publish_gap_, check_lidar_time_, check_odom_time_;
    std::chrono::time_point<std::chrono::steady_clock> check_data_processing_time_;

    float scan_min_dis_;

    pcl::PointCloud<pcl::PointXYZI> laser_cloud_;

    std::ofstream ofs_;
    std::string log_path_, err_path_;

};

#endif // RLML_ROS_H
