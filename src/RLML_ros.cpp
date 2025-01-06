
#include "RLML_ros.h"

fusionSlamRos::fusionSlamRos() : got_lidar_tf_(false), initialize_(false),
relocalize_(false),  tracking_succeed_(true), receive_first_lidar_(false),
receive_first_odom_(false), publish_time_(0), check_lidar_time_(0), check_odom_time_(0)
{
    ros::NodeHandle nh;
    ros::NodeHandle private_nh("~");

    std::string scan_topic, map_topic, odom_topic;
    float relocalize_min_response;
    private_nh.param("scan_topic", scan_topic, std::string("scan"));
    private_nh.param("map_topic", map_topic, std::string("map"));
    private_nh.param("odom_topic", odom_topic, std::string("odom"));

    private_nh.param("base_frame", base_frame_, std::string("base_footprint"));
    private_nh.param("map_frame", map_frame_, std::string("map"));
    private_nh.param("lidar_frame", lidar_frame_, std::string("sensor_Link"));
    private_nh.param("use_publish_pose_thread", use_publish_pose_thread_, false);
    private_nh.param("publish_map_freq", publish_freq_, 0.25);
    private_nh.param("publish_pose_freq", pose_freq_, 60.);
    private_nh.param("mapping", mapping_, false);
    private_nh.param("use_reflector", use_reflector_, true);
    private_nh.param("use_reflector_relocalize", use_reflector_relocalize_, true);
    private_nh.param("preset_lidar_tf", preset_lidar_tf_, true);
    private_nh.param("lidar_trans_x", lidar_trans_x_, 0.857f);
    private_nh.param("lidar_trans_y", lidar_trans_y_, 0.0f);
    private_nh.param("lidar_yaw", lidar_yaw_, 0.0f);
    private_nh.param("reflector_max_length", reflector_max_length_, 0.1);
    private_nh.param("reflector_min_length", reflector_min_length_, 0.03);
    private_nh.param("intensity_min", intensity_min_, 500.);
    private_nh.param("relocalize_min_response", relocalize_min_response, 0.5f);
    private_nh.param("scan_min_dis", scan_min_dis_, 0.4f);

    std::cout << "fusion slam 3.2.0" <<std::endl;

    if (preset_lidar_tf_)
    {
        Eigen::Vector3f lidar_to_odom(lidar_trans_x_, lidar_trans_y_, lidar_yaw_);
        std::cout << "lidar_to_odom " << std::endl << lidar_to_odom << std::endl;
        lidar_to_odom_ = lidar_to_odom;

        tf::Transform lidar_tf(tf::createQuaternionFromYaw(lidar_to_odom_[2]), tf::Vector3(lidar_to_odom_[0], lidar_to_odom_[1], 0.0f));
        lidar_tf_ = tf::StampedTransform(lidar_tf, ros::Time::now(), base_frame_, lidar_frame_);
        tf::transformTFToEigen(lidar_tf_, lidar_trans_);

        static tf2_ros::StaticTransformBroadcaster static_tf_broadcaster;
        geometry_msgs::TransformStamped static_transform_stamped;

        static_transform_stamped.header.stamp = ros::Time::now();
        static_transform_stamped.header.frame_id = base_frame_;
        static_transform_stamped.child_frame_id = lidar_frame_;
        static_transform_stamped.transform.translation.x = lidar_to_odom_[0];
        static_transform_stamped.transform.translation.y = lidar_to_odom_[1];
        static_transform_stamped.transform.translation.z = 0.0;
        tf::Quaternion quat = tf::createQuaternionFromYaw(lidar_to_odom_[2]);
        static_transform_stamped.transform.rotation.x = quat.x();
        static_transform_stamped.transform.rotation.y = quat.y();
        static_transform_stamped.transform.rotation.z = quat.z();
        static_transform_stamped.transform.rotation.w = quat.w();

        static_tf_broadcaster.sendTransform(static_transform_stamped);
    }

    private_nh.param("log_path", log_path_, std::string("log/log.txt"));
    private_nh.param("err_path", err_path_, std::string("log/err.txt"));

    // uncomment to log into log files
    if (freopen(log_path_.c_str(), "w", stdout) == nullptr) {
        std::cerr << "Failed to redirect stdout to log file: " << log_path_ << std::endl;
        std::exit(EXIT_FAILURE);
    }

    if (freopen(err_path_.c_str(), "w", stderr) == nullptr) {
        std::cerr << "Failed to redirect stderr to error file: " << err_path_ << std::endl;
        std::exit(EXIT_FAILURE);
    }


    std::cout << "fusion slam 3.2.0" <<std::endl;
    std::cerr << "fusion slam 3.2.0" <<std::endl;
    publish_gap_ = 1. / pose_freq_;

    std::cout << "mapping mode: " << mapping_ << std::endl;
    std::cout << "use reflector: " << use_reflector_ << std::endl;
    std::cout << "use reflector relocalize: " << use_reflector_relocalize_ << std::endl;
    std::cout << "reflector_max_length: " << reflector_max_length_ << std::endl;
    std::cout << "reflector_min_length: " << reflector_min_length_ << std::endl;
    std::cout << "intensity_min: " << intensity_min_ << std::endl;
    std::cout << "relocalize_min_response: " << relocalize_min_response << std::endl;
    std::cout << "localization frequency: " << pose_freq_ << std::endl;
    reflector_mean_length_ = (reflector_min_length_ + reflector_max_length_)/2.;


    private_nh.param("save_path", save_path_, std::string("log/"));
    private_nh.param("load_path", load_path_, std::string("log/"));
    std::cout << "save path: " << save_path_ << std::endl;
    std::cout << "load path: " << load_path_ << std::endl;

    double v_covariance_param = 0, omega_covariance_param = 0, scan_probability_covariance_param = 0;
    float err_cost = 0.15, alert_cost = 0.25;
    Eigen::Matrix3d initial_pose_covariance = Eigen::Matrix3d::Zero();
    initial_pose_covariance <<
            0.05 * 0.05, 0, 0,
            0, 0.05 * 0.05, 0,
            0, 0, 0.03 * 0.03;
    Eigen::Matrix3d odom_model_covariance = Eigen::Matrix3d::Zero();
    odom_model_covariance <<
            0.05 * 0.05, 0, 0,
            0, 0.05 * 0.05, 0,
            0, 0, 0.01 * 0.01;
    Eigen::Matrix2d reflector_covariance = Eigen::Matrix2d::Zero();
    reflector_covariance <<
            0.025 * 0.025, 0,
            0, 0.025 * 0.025;
    Eigen::Matrix3d relocalize_covariance = Eigen::Matrix3d::Zero();
    relocalize_covariance <<
            0.03 * 0.03, 0, 0,
            0, 0.03 * 0.03, 0,
            0, 0, 0.01 * 0.01;
    Eigen::Matrix<float, 1, 3> initial_pose_trans = Eigen::Matrix<float, 1, 3>::Zero();

    private_nh.param("err_cost", err_cost, 0.15f);
    private_nh.param("alert_cost", alert_cost, 0.25f);

    private_nh.param("v_covariance_param", v_covariance_param, 0.1);
    private_nh.param("omega_covariance_param", omega_covariance_param, 0.05);
    private_nh.param("scan_probability_covariance_param", scan_probability_covariance_param, 0.5);

    // load initial_pose_covariance
    loadCovariance(private_nh, "initial_pose_covariance", initial_pose_covariance);
    // load odom_model_covariance
    loadCovariance(private_nh, "odom_model_covariance", odom_model_covariance);
    // load reflector_covariance
    loadCovariance(private_nh, "reflector_covariance", reflector_covariance);
    // load relocalize_covariance
    loadCovariance(private_nh, "relocalize_covariance", relocalize_covariance);
    // load initial pose
    loadCovariance(private_nh, "initial_pose", initial_pose_trans);
    initial_pose_ = initial_pose_trans.transpose();



    map_builder_.setCostThreshold(err_cost, alert_cost);
    map_builder_.setInitialCovariance(initial_pose_covariance);
    map_builder_.setOdomCovariance(odom_model_covariance);
    map_builder_.setReflectorCovariance(reflector_covariance);
    map_builder_.setRelocalizeCovariance(relocalize_covariance);
    map_builder_.setVelocityCovarianceParam(v_covariance_param);
    map_builder_.setOmegaCovarianceParam(omega_covariance_param);
    map_builder_.setScanCovariance(scan_probability_covariance_param);
    map_builder_.setRelocalizeMinResponse(relocalize_min_response);

    double occupancy_grid_map_resolution;
    double min_update_distance, min_update_orientation;
    int scan_buffer_size;
    double loop_scan_search_distance;
    int loop_match_min_chain_size;
    double loop_closure_min_response;
    double loop_closure_xy_variance_threshold, loop_closure_angle_variance_threshold;
    int optimize_every_n_constrains;

    double loop_closure_xy_search_range, loop_closure_angle_search_range;
    double loop_closure_grid_resolution;
    double loop_closure_coarse_xy_search_resolution, loop_closure_fine_xy_search_range;
    double loop_closure_coarse_angle_search_resolution;
    double loop_closure_fine_angle_search_range, loop_closure_fine_angle_search_resolution;

    bool use_correlative_scan_matcher;


    private_nh.param("occupancy_grid_map_resolution", occupancy_grid_map_resolution, 0.05);
    private_nh.param("min_update_distance", min_update_distance, 0.15);
    private_nh.param("min_update_orientation", min_update_orientation, 5.0);
    private_nh.param("scan_buffer_size", scan_buffer_size, 30);
    private_nh.param("loop_scan_search_distance", loop_scan_search_distance, 10.0);
    private_nh.param("loop_match_min_chain_size", loop_match_min_chain_size, 5);
    private_nh.param("loop_closure_min_response", loop_closure_min_response, 0.65);
    private_nh.param("loop_closure_xy_variance_threshold", loop_closure_xy_variance_threshold, 0.01);
    private_nh.param("loop_closure_angle_variance_threshold", loop_closure_angle_variance_threshold, 0.05);
    private_nh.param("optimize_every_n_constrains", optimize_every_n_constrains, 20);

    private_nh.param("loop_closure_xy_search_range", loop_closure_xy_search_range, 5.0);
    private_nh.param("loop_closure_angle_search_range", loop_closure_angle_search_range, 20.0);
    private_nh.param("loop_closure_grid_resolution", loop_closure_grid_resolution, 0.05);
    private_nh.param("loop_closure_coarse_xy_search_resolution", loop_closure_coarse_xy_search_resolution, 0.1);
    private_nh.param("loop_closure_fine_xy_search_range", loop_closure_fine_xy_search_range, 0.1);
    private_nh.param("loop_closure_coarse_angle_search_resolution", loop_closure_coarse_angle_search_resolution, 2.0);
    private_nh.param("loop_closure_fine_angle_search_range", loop_closure_fine_angle_search_range, 2.0);
    private_nh.param("loop_closure_fine_angle_search_resolution", loop_closure_fine_angle_search_resolution, 0.2);

    map_builder_.setOccupancyGridMapResolution(occupancy_grid_map_resolution);
    map_builder_.setMinUpdateDistance(min_update_distance);
    map_builder_.setMinUpdateOrientation(RLML::degToRad(min_update_orientation));
    map_builder_.setScanBufferSize(scan_buffer_size);
    map_builder_.setLoopScanSearchDistance(loop_scan_search_distance);
    map_builder_.setLoopMatchMinChainSize(loop_match_min_chain_size);
    map_builder_.setLoopClosureMinResponse(loop_closure_min_response);
    map_builder_.setLoopClosureXYVarianceThreshold(loop_closure_xy_variance_threshold);
    map_builder_.setLoopClosureAngleVarianceThreshold(loop_closure_angle_variance_threshold);
    map_builder_.setOptimizeEveryNConstraints(optimize_every_n_constrains);

    map_builder_.setLoopClosureXYSearchRange(loop_closure_xy_search_range);
    map_builder_.setLoopClosureAngleSearchRange(RLML::degToRad(loop_closure_angle_search_range));
    map_builder_.setLoopClosureGridResolution(loop_closure_grid_resolution);
    map_builder_.setLoopClosureCoarseXYSearchResolution(loop_closure_coarse_xy_search_resolution);
    map_builder_.setLoopClosureFineXYSearchRange(loop_closure_fine_xy_search_range);
    map_builder_.setLoopClosureCoarseAngleSearchResolution(RLML::degToRad(loop_closure_coarse_angle_search_resolution));
    map_builder_.setLoopClosureFineAngleSearchRange(RLML::degToRad(loop_closure_fine_angle_search_range));
    map_builder_.setLoopClosureFineAngleSearchResolution(RLML::degToRad(loop_closure_fine_angle_search_resolution));

    private_nh.param("use_correlative_scan_matcher", use_correlative_scan_matcher, false);
    if(use_correlative_scan_matcher) {
        double csm_xy_search_range;
        double csm_angle_search_range;
        double csm_grid_resolution;
        double csm_xy_search_resolution;
        double csm_angle_search_resolution;

        private_nh.param("csm_xy_search_range", csm_xy_search_range, 0.2);
        private_nh.param("csm_angle_search_range", csm_angle_search_range, 20.0);
        private_nh.param("csm_grid_resolution", csm_grid_resolution, 0.01);
        private_nh.param("csm_xy_search_resolution", csm_xy_search_resolution, 0.02);
        private_nh.param("csm_angle_search_resolution", csm_angle_search_resolution, 0.5);

        map_builder_.useCorrelativeScanMatcher(true);
        map_builder_.setCSMXYSearchRange(csm_xy_search_range);
        map_builder_.setCSMAngleSearchRange(RLML::degToRad(csm_angle_search_range));
        map_builder_.setCSMGridResolution(csm_grid_resolution);
        map_builder_.setCSMXYSearchResolution(csm_xy_search_resolution);
        map_builder_.setCSMAngleSearchResolution(RLML::degToRad(csm_angle_search_resolution));
    }





    map_builder_.initialize();
    if(!mapping_)
    {
        if (!map_builder_.loadScans(load_path_)) {
            std::cerr << "Failed to load scans. Exiting..." << std::endl;
            std::exit(EXIT_FAILURE);
        }
        std::cout << "initial pose " << std::endl << initial_pose_ << std::endl;
    }

    map_pub_ = nh.advertise<nav_msgs::OccupancyGrid>("map", 1, true);
    correlative_grid_pub_ = nh.advertise<nav_msgs::OccupancyGrid>("correlative_grid", 1, true);
    pose_pub_ = nh.advertise<geometry_msgs::PoseStamped>("robot_pose", 1, true);
    path_pub_ = nh.advertise<nav_msgs::Path>("path", 1, true);
    constraint_list_pub_ = nh.advertise<visualization_msgs::MarkerArray>("constraint_list", 1, true);
    landmarks_pub_ = nh.advertise<visualization_msgs::MarkerArray>( "landmarks", 1);
    cur_landmarks_pub_ = nh.advertise<visualization_msgs::MarkerArray>( "cur_landmarks", 1);
    odom_sub_ = nh.subscribe(odom_topic, 10, &fusionSlamRos::odomCallback, this);


    lidar_subscriber_.subscribe(nh, scan_topic, 10);
    odom_subscriber_.subscribe(nh, odom_topic, 10);
    sync_.reset(new Sync(syncPolicy(20), lidar_subscriber_, odom_subscriber_));
    sync_->registerCallback(std::bind(&fusionSlamRos::scanCallback, this, std::placeholders::_1, std::placeholders::_2));
    if(!mapping_)
    {
        std::cout << "publishing map, wait..." << std::endl;
        publishConstraintList();
        publishProbabilityGridMap();
        publishPath();
        std::cout << "publishing map done..." << std::endl;
    }
    publish_thread_ = std::make_shared<std::thread>(&fusionSlamRos::publishLoop, this);


    check_data_processing_time_ = std::chrono::steady_clock::now();
    proc_th_ptr_ = std::make_shared<std::thread>(&fusionSlamRos::processData, this);
    cmdThPtr_ = std::make_shared<std::thread>(&fusionSlamRos::cmdProcess, this);

    if(use_publish_pose_thread_)
    {
        poseThPtr_ = std::make_shared<std::thread>(&fusionSlamRos::publishPoseUpToDate, this);
    }

    optimization_srv_ = nh.advertiseService("optimization", &fusionSlamRos::optimizationCallback, this);

    setCorrelativeTranslationTable();
}

fusionSlamRos::~fusionSlamRos()
{
    fclose(stdout);
    cmdThPtr_->join();
    proc_th_ptr_->join();
    publish_thread_->join();
    if(use_publish_pose_thread_)
    {
        poseThPtr_->join();
    }
    if(correlative_translation_table_) {
        delete[] correlative_translation_table_;
    }
}

void fusionSlamRos::setCorrelativeTranslationTable()
{
    correlative_translation_table_ = new char[256];

    for (int i = 0; i < 256; i++) {
        correlative_translation_table_[i] = char(i * 100 / 255);
    }
}

bool fusionSlamRos::optimizationCallback(std_srvs::Empty::Request& req, std_srvs::Empty::Response& res)
{
    map_builder_.doPoseAdjustment();
    return true;
}





void fusionSlamRos::scanCallback(const sensor_msgs::LaserScanConstPtr& scan_msg,
                               const nav_msgs::Odometry::ConstPtr& odom_msg)
{
    if(receive_first_lidar_)
    {
        if(scan_msg->header.stamp.toSec() - check_lidar_time_ > 0.2)
        {
            std::cerr << "lidar and odom collect process abnormal" << std::endl;
            std::cout << "lidar and odom collect process abnormal" << std::endl;

        }
    } else
    {
        std::cout << "receive first lidar and odom" << std::endl;
        receive_first_lidar_ = true;
    }

    check_lidar_time_ = scan_msg->header.stamp.toSec();

    if(!preset_lidar_tf_ && !got_lidar_tf_) {
        try {
            tf_listener_.waitForTransform(base_frame_, lidar_frame_, ros::Time(0), ros::Duration(0.5));
            tf_listener_.lookupTransform(base_frame_, lidar_frame_, ros::Time(0), lidar_tf_);
            tf::transformTFToEigen(lidar_tf_, lidar_trans_);
            Eigen::Matrix4f lidar_trans_mat = lidar_trans_.matrix().cast<float>();
            Eigen::Matrix3f lidar_to_odom_rot = lidar_trans_mat.block(0,0,3,3);
            Eigen::Vector3f euler = lidar_to_odom_rot.eulerAngles(2,1,0);
            float yaw = euler(0,0);
            Eigen::Vector3f lidar_to_odom(lidar_trans_mat(0,3), lidar_trans_mat(1,3), yaw);
            std::cout << "lidar_to_odom " << std::endl << lidar_to_odom << std::endl;
            lidar_to_odom_ = lidar_to_odom;
            got_lidar_tf_ = true;
        }
        catch (const tf::TransformException& e) {
            ROS_ERROR("%s", e.what());
            return;
        }
    }


    sycData datafm;
    datafm.scan_msg = scan_msg;
    datafm.odom_msg = odom_msg;

    mBuf_.lock();
    dataques_.push_back(datafm);
    std::cout << "get a new frame, frame size: " << dataques_.size() << std::endl;
    mBuf_.unlock();
}


void fusionSlamRos::odomCallback(const nav_msgs::Odometry::ConstPtr& odom_msg)
{
    if(receive_first_odom_)
    {
        if(odom_msg->header.stamp.toSec() - check_odom_time_ > 0.2)
        {
            std::cerr << "odom datastream not stable, stuck for a long time, timestamp: " << std::setprecision(18) << odom_msg->header.stamp.toSec() << std::endl;
            std::cout << "odom datastream not stable, stuck for a long time, timestamp: " << std::setprecision(18) << odom_msg->header.stamp.toSec() << std::endl;

        }
    } else
    {
        std::cout << "receive first odom" << std::endl;
        receive_first_odom_ = true;
    }

    check_odom_time_ = odom_msg->header.stamp.toSec();
}



void fusionSlamRos::processData()
{
    while (ros::ok())
    {
        auto t1 = std::chrono::steady_clock::now();
        auto process_gap = std::chrono::duration_cast<std::chrono::duration<double>>(t1 - check_data_processing_time_).count() * 1000.;
        if( process_gap > 200)
        {
            std::stringstream ss;
            ss << "data processing thread abnormal gap " << process_gap << std::endl;
            std::cout << ss.str() << std::endl;
        }
        check_data_processing_time_ = t1;
        mBuf_.lock();
        if(!dataques_.empty())
        {
            sycData datafm = dataques_.front();
            dataques_.pop_front();
            std::cout << "process a frame, frame size: " << dataques_.size() << std::endl;

            if(!mapping_ && !relocalize_ && !dataques_.empty())
            {
                while(!dataques_.empty()) dataques_.pop_front();
                std::cout << "empty the queue when relocalizing" << dataques_.size() << std::endl;
            }
            mBuf_.unlock();

            // process data

            sensor_msgs::LaserScanConstPtr scan_msg = datafm.scan_msg;
            nav_msgs::Odometry::ConstPtr odom_msg = datafm.odom_msg;



            // All reflector points
            std::vector<std::pair<float,Eigen::Vector2f>> reflector_points;
            // All reflector centers
            std::vector<Eigen::Vector2f> centers;
            // Now reflector points, ids and center
            std::vector<Eigen::Vector2f> reflector;
            std::vector<int> reflector_id;
            Eigen::Vector2f now_center(0.f, 0.f);

            PointCloud scan_points;


            laser_cloud_.clear();
            for(size_t i = 0; i < scan_msg->ranges.size(); ++i) {
                float angle = scan_msg->angle_min + i * scan_msg->angle_increment;
                if(scan_msg->ranges[i] > scan_min_dis_ && scan_msg->ranges[i] < scan_msg->range_max) {
                    Eigen::Vector3f scan_point(
                            scan_msg->ranges[i] * std::cos(angle), scan_msg->ranges[i] * std::sin(angle), 0.0);

                    scan_point = static_cast<Eigen::Affine3f>(lidar_trans_) * scan_point;
                    Eigen::Vector2f laser_point(scan_point.x(), scan_point.y());


                    scan_points.emplace_back(laser_point);


                    pcl::PointXYZI point;
                    point.x = scan_point.x();
                    point.y = scan_point.y();
                    point.z = 0;
                    point.intensity = scan_msg->intensities[i];
                    laser_cloud_.push_back(point);

                    if(scan_msg->intensities[i] > intensity_min_)
                    {
                        // Add the first point
                        if (reflector.empty())
                        {
                            reflector.push_back(laser_point);
                            reflector_id.push_back(i);
                            now_center += laser_point;
                        }
                        else
                        {
                            const int last_id = reflector_id.back();
                            // Add connected points
                            if (i - last_id == 1)
                            {
                                reflector.push_back(laser_point);
                                reflector_id.push_back(i);
                                now_center += laser_point;
                            }
                            else
                            {
                                // Calculate reflector length
                                const float reflector_length = std::hypotf(reflector.front().x() - reflector.back().x(),
                                                                           reflector.front().y() - reflector.back().y());

//                        std::cout << "reflector_length "  << reflector_length << std::endl;

                                // Add good reflector and its center
                                if ((reflector_length > reflector_min_length_) && (reflector_length < reflector_max_length_))
                                {
                                    Eigen::Vector2f center = now_center / reflector.size();
//                                    centers.push_back(center);
                                    reflector_points.push_back(std::make_pair(std::abs(reflector_length - reflector_mean_length_), center));
                                }
                                // Update now reflector
                                reflector.clear();
                                reflector_id.clear();
                                now_center.setZero(2);
                                reflector_id.push_back(i);
                                reflector.push_back(laser_point);
                                now_center += laser_point;
                            }
                        }
                    }
                }
            }


            if (!reflector.empty())
            {
                const float reflector_length = std::hypotf(reflector.front().x() - reflector.back().x(),
                                                           reflector.front().y() - reflector.back().y());
                if ((reflector_length > reflector_min_length_) && (reflector_length < reflector_max_length_))
                {
                    Eigen::Vector2f center = now_center / reflector.size();
//                    centers.push_back(center);
                    reflector_points.push_back(std::make_pair(std::abs(reflector_length - reflector_mean_length_), center));
                }
            }
            std::shared_ptr<RLML::LaserScan> laser_scan(new RLML::LaserScan(scan_points));
            laser_scan->setCalibration(lidar_to_odom_);

            if(!use_reflector_) reflector_points.clear();

            std::sort(reflector_points.begin(), reflector_points.end(), [](auto &left, auto &right) {
                return left.first < right.first;
            });

            for(const auto& reflector:reflector_points)
            {
                centers.push_back(reflector.second);
            }
            std::stringstream ss;
            ss << "detect reflectors " << centers.size() << " ";

            if(!mapping_ && !relocalize_)
            {
                laser_scan->setPose(initial_pose_);
                bool automatic_relocalize = use_reflector_ & use_reflector_relocalize_;
                if(map_builder_.relocalize(laser_scan, centers, automatic_relocalize))
                {
                    relocalize_ = true;
                    initial_pose_ = laser_scan->getPose();
                } else
                {
                    if(automatic_relocalize)
                    {
                        std::cout << "auto relocalize fail" << std::endl;
                    } else
                    {
                        std::cout << "manual relocalize fail" << std::endl;
                    }

                    continue;
                }
            }

            if(!mapping_ && !tracking_succeed_)
            {
                std::cout << "tracking fail " << std::endl;

                continue;
            }


            float v = 0, omega = 0;
            omega = (float)odom_msg->twist.twist.angular.z;
            v = (float)odom_msg->twist.twist.linear.x;
            ss << "v " << v << " omega " << omega << std::endl;
            std::cout << ss.str();
            map_builder_.addOdom(v, omega, (float)(scan_msg->header.stamp.toSec() - previous_t_));
            previous_t_ = scan_msg->header.stamp.toSec();


            mPub_.lock();
            if(!mapping_)
            {
                tracking_succeed_ = map_builder_.localizeLaserScan(initial_pose_, laser_scan, centers);
            } else
            {
                map_builder_.addLaserScan(laser_scan, centers);
            }
            Eigen::Vector3f pose = laser_scan->getPose();
            cur_state_ = robotState(v, omega, pose, scan_msg->header.stamp);
            mPub_.unlock();


            if(!initialize_) {
                initialize_ = true;
            }

            if(!use_publish_pose_thread_)
            {
                publishPose(pose, scan_msg->header.stamp);
                tf::Transform base_to_map(tf::createQuaternionFromYaw(pose[2]), tf::Vector3(pose[0], pose[1], 0.0f));
                tf_broadcaster_.sendTransform(tf::StampedTransform(base_to_map, scan_msg->header.stamp, map_frame_, base_frame_));

            }

            if(!mapping_)
            {
                publishCurLandmarks(scan_msg->header.stamp);

            }
        }
        else
        {
            mBuf_.unlock();
        }

        std::chrono::milliseconds dura(1);
        std::this_thread::sleep_for(dura);
    }
}




void fusionSlamRos::cmdProcess()
{

    ros::Rate rate(100);
    while (ros::ok())
    {

        char c = (char)getchar();
        if(c == 's' && mapping_)
        {
            mPub_.lock();
            printf("save map...\n");
            saveMap();
            mPub_.unlock();
            printf("program shutting down...\n");
            ros::shutdown();
        }
        if(c == 'u')
        {
            if(mapping_)
            {
                mPub_.lock();
                printf("update map...\n");
                updateMap();
                mPub_.unlock();
            } else
            {
                mPub_.lock();
                printf("update map...\n");
                map_builder_.discardScan();
                publishProbabilityGridMap();
                mPub_.unlock();
            }

        }
        if(c == 'q')
        {
            printf("program shutting down...\n");

            ros::shutdown();
        }


        rate.sleep();
    }
}



void fusionSlamRos::publishPoseUpToDate()
{
//    ros::Rate rate(pose_freq_);

    while (ros::ok()) {
        if(initialize_ && use_publish_pose_thread_) {
            double publish_time;

            publish_time = ros::Time::now().toSec();


            if(publish_time_ < 1.)
            {
                publish_time_ = publish_time;
            }
            else if(publish_time - publish_time_ > publish_gap_)
            {
                publish_time_ = publish_time;
                mBuf_.lock();
                robotState cur_state = cur_state_;
                Eigen::Vector3f cur_pose = cur_state.pose;
                std::deque<sycData> dataques = dataques_;
                mBuf_.unlock();
                std::stringstream ss;
                ss << "current stamp: " << std::setprecision(18) << cur_state.stamp.toSec() << " update pose stamp: " << std::setprecision(14) << publish_time - cur_state.stamp.toSec() << std::endl;
                std::cout << ss.str();
                for (size_t i = 0; i < dataques.size(); ++i) {
                    float dt = (float)(dataques[i].scan_msg->header.stamp.toSec() - cur_state.stamp.toSec());
                    if(dt < 0.001)
                    {
                        continue;
                    }
                    cur_state.stamp = dataques[i].scan_msg->header.stamp;
//                std::cout << "update pose stamp: " << std::setprecision(14) << cur_state.stamp.toSec() << std::endl;
                    cur_state.v = (float)dataques[i].odom_msg->twist.twist.linear.x;
                    cur_state.omega = (float)dataques[i].odom_msg->twist.twist.angular.z;
                    float angular_half_delta = cur_state.pose[2] + cur_state.omega * dt / 2.f;
                    Eigen::Vector3f increment;
                    if(cur_state.omega > 0.01f)
                    {
                        increment[0] = (cur_state.v / cur_state.omega) * (std::sin(cur_state.pose[2] + cur_state.omega * dt) - std::sin(cur_state.pose[2]));
                        increment[1] = (cur_state.v / cur_state.omega) * std::cos(cur_state.pose[2])  - (cur_state.v / cur_state.omega) * std::cos(cur_state.pose[2] + cur_state.omega * dt);
                        increment[2] =  cur_state.omega * dt;
                    } else
                    {
                        increment[0] = cur_state.v * dt * std::cos(angular_half_delta);
                        increment[1] = cur_state.v * dt * std::sin(angular_half_delta);
                        increment[2] =  cur_state.omega * dt;
                    }
                    cur_state.pose = cur_state.pose + increment;
                }
                ros::Time cur_stamp = ros::Time(publish_time);
                float dt = (float)(cur_stamp.toSec() - cur_state.stamp.toSec());
//            std::cout << "current stamp: " << std::setprecision(14) << cur_stamp.toSec() << std::endl;

                if(dt < 0.0)
                {
                    cur_stamp = cur_state.stamp;
                    cur_pose = cur_state.pose;
                    publishPose(cur_pose, cur_stamp);
                    tf::Transform map_to_base(tf::createQuaternionFromYaw(cur_pose[2]), tf::Vector3(cur_pose[0], cur_pose[1], 0.0f));
                    tf_broadcaster_.sendTransform(tf::StampedTransform(map_to_base, cur_stamp, map_frame_, base_frame_));
                }
                else if(dt > 0.2)
                {
                    std::cout << "data missing, cannot update pose" << std::endl;
                }
                // 0 < dt < 0.2
                else if(!tracking_succeed_)
                {
                    mBuf_.lock();
                    mPub_.lock();
                    std::cerr << "Cannot update pose. Might be sensor data stream not stable or tracking fail" << std::setprecision(18) << ros::Time::now().toSec() << std::endl;
                    std::cout << "Cannot update pose. Might be sensor data stream not stable or tracking fail" <<  std::setprecision(18) << ros::Time::now().toSec() << std::endl;
                    ros::shutdown();
                    mPub_.unlock();
                    mBuf_.unlock();

                }
                else
                {
                    float angular_half_delta = cur_state.pose[2] + cur_state.omega * dt / 2.f;
                    Eigen::Vector3f increment;
                    if(cur_state.omega > 0.01f)
                    {
                        increment[0] = (cur_state.v / cur_state.omega) * (std::sin(cur_state.pose[2] + cur_state.omega * dt) - std::sin(cur_state.pose[2]));
                        increment[1] = (cur_state.v / cur_state.omega) * std::cos(cur_state.pose[2])  - (cur_state.v / cur_state.omega) * std::cos(cur_state.pose[2] + cur_state.omega * dt);
                        increment[2] =  cur_state.omega * dt;
                    } else
                    {
                        increment[0] = cur_state.v * dt * std::cos(angular_half_delta);
                        increment[1] = cur_state.v * dt * std::sin(angular_half_delta);
                        increment[2] =  cur_state.omega * dt;
                    }
                    cur_pose = cur_state.pose + increment;
                    publishPose(cur_pose, cur_stamp);
                    tf::Transform map_to_base(tf::createQuaternionFromYaw(cur_pose[2]), tf::Vector3(cur_pose[0], cur_pose[1], 0.0f));
                    tf_broadcaster_.sendTransform(tf::StampedTransform(map_to_base, cur_stamp, map_frame_, base_frame_));

                }

            }

        }
//        rate.sleep();
        std::chrono::milliseconds dura(1);
        std::this_thread::sleep_for(dura);
    }
}

void fusionSlamRos::publishPath()
{
    std::vector<Eigen::Vector3f> path = map_builder_.getPath();
    nav_msgs::Path path_msg;

    path_msg.header.frame_id = map_frame_;
    path_msg.header.stamp = ros::Time::now();

    for(const Eigen::Vector3f& pose : path) {
        geometry_msgs::PoseStamped pose_msg;
        pose_msg.header.stamp = ros::Time::now();
        pose_msg.header.frame_id = map_frame_;
        pose_msg.pose.position.x = pose[0];
        pose_msg.pose.position.y = pose[1];
        pose_msg.pose.position.z = 0;

        tf::Quaternion q;
        q.setRPY(0.0, 0.0, pose[2]);
        pose_msg.pose.orientation.x = q.x();
        pose_msg.pose.orientation.y = q.y();
        pose_msg.pose.orientation.z = q.z();
        pose_msg.pose.orientation.w = q.w();
        path_msg.poses.push_back(pose_msg);
    }

    path_pub_.publish(path_msg);
}

void fusionSlamRos::publishPose(const Eigen::Vector3f& pose, const ros::Time& stamp)
{
    geometry_msgs::PoseStamped pose_msg;
    pose_msg.header.stamp = stamp;
    pose_msg.header.frame_id = map_frame_;
    pose_msg.pose.position.x = pose[0];
    pose_msg.pose.position.y = pose[1];
    pose_msg.pose.position.z = 0;

    tf::Quaternion q;
    q.setRPY(0.0, 0.0, pose[2]);
    pose_msg.pose.orientation.x = q.x();
    pose_msg.pose.orientation.y = q.y();
    pose_msg.pose.orientation.z = q.z();
    pose_msg.pose.orientation.w = q.w();
    pose_pub_.publish(pose_msg);
}

void fusionSlamRos::publishOccupancyGridMap()
{
    std::shared_ptr<RLML::OccupancyGridMap> map = map_builder_.getOccupancyGridMap();

    nav_msgs::OccupancyGrid map_msg;
    Eigen::Vector2f origin = map->getOrigin();
    map_msg.header.stamp = ros::Time::now();
    map_msg.info.origin.position.x = origin.x();
    map_msg.info.origin.position.y = origin.y();
    map_msg.info.origin.orientation.x = 0;
    map_msg.info.origin.orientation.y = 0;
    map_msg.info.origin.orientation.z = 0;
    map_msg.info.origin.orientation.w = 1;
    map_msg.info.resolution = map->getResolution();
    map_msg.info.width = map->getSizeX();
    map_msg.info.height = map->getSizeY();
    map_msg.data.resize(map_msg.info.width * map_msg.info.height, -1);

    for(size_t i = 0; i < map_msg.data.size(); ++i) {
        if(map->isFree(i)) {
            map_msg.data[i] = 0;
        }
        else if(map->isOccupied(i)) {
            map_msg.data[i] = 100;
        }
    }

    map_pub_.publish(map_msg);
}

void fusionSlamRos::publishProbabilityGridMap()
{
    std::shared_ptr<RLML::ProbabilityGridMap> map = map_builder_.getProbabilityGridMap();

    nav_msgs::OccupancyGrid map_msg;
    Eigen::Vector2f origin = map->getOrigin();
    map_msg.header.stamp = ros::Time::now();
    map_msg.info.origin.position.x = origin.x();
    map_msg.info.origin.position.y = origin.y();
    map_msg.info.origin.orientation.x = 0;
    map_msg.info.origin.orientation.y = 0;
    map_msg.info.origin.orientation.z = 0;
    map_msg.info.origin.orientation.w = 1;
    map_msg.info.resolution = map->getResolution();
    map_msg.info.width = map->getSizeX();
    map_msg.info.height = map->getSizeY();
    map_msg.data.resize(map_msg.info.width * map_msg.info.height, -1);

    for(size_t i = 0; i < map_msg.data.size(); ++i) {
        int value = map->getGridValue(i);
        if(value == RLML::LogOdds_Unknown) {
            map_msg.data[i] = -1;
        }
        else {
            map_msg.data[i] = map->getGridValue(i);
        }
    }

    map_pub_.publish(map_msg);
}

void fusionSlamRos::publishCorrelativeGrid()
{  
    std::shared_ptr<RLML::CorrelativeGrid> correlative_grid = map_builder_.getCorrelativeGrid();

    if(correlative_grid->getSize() == 0) {
        return;
    }

    nav_msgs::OccupancyGrid map_msg;
    Eigen::Vector2f origin = correlative_grid->getOrigin();
    map_msg.header.stamp = ros::Time::now();
    map_msg.info.origin.position.x = origin.x();
    map_msg.info.origin.position.y = origin.y();
    map_msg.info.origin.orientation.x = 0;
    map_msg.info.origin.orientation.y = 0;
    map_msg.info.origin.orientation.z = 0;
    map_msg.info.origin.orientation.w = 1;
    map_msg.info.resolution = correlative_grid->getResolution();
    map_msg.info.width = correlative_grid->getSizeX();
    map_msg.info.height = correlative_grid->getSizeY();
    map_msg.data.resize(map_msg.info.width * map_msg.info.height, 0);

    for(size_t i = 0; i < map_msg.data.size(); ++i) {
        map_msg.data[i] = correlative_translation_table_[correlative_grid->getGridValue(i)];
    }

    correlative_grid_pub_.publish(map_msg);
}

void fusionSlamRos::publishConstraintList()
{
    std::vector<Eigen::Vector2f> graph_nodes;
    std::vector<std::pair<Eigen::Vector2f, Eigen::Vector2f>> graph_edges;
    map_builder_.getGraph(graph_nodes, graph_edges);

    visualization_msgs::MarkerArray marker_array;

    visualization_msgs::Marker marker;
    marker.header.frame_id = map_frame_;
    marker.header.stamp = ros::Time::now();
    marker.action = visualization_msgs::Marker::ADD;
    marker.id = 0;
    marker.type = visualization_msgs::Marker::SPHERE;
    marker.pose.position.x = 0.0;
    marker.pose.position.y = 0.0;
    marker.pose.position.z = 0.0;
    marker.pose.orientation.w = 1.0;
    marker.pose.orientation.x = 0.0;
    marker.pose.orientation.y = 0.0;
    marker.pose.orientation.z = 0.0;
    marker.scale.x = 0.01;
    marker.scale.y = 0.01;
    marker.scale.z = 0.01;
    marker.color.r = 0.0;
    marker.color.g = 0.0;
    marker.color.b = 0.5;
    marker.color.a = 1.0;
    marker.lifetime = ros::Duration(0);

    visualization_msgs::Marker edge;
    edge.header.frame_id = map_frame_;
    edge.header.stamp = ros::Time::now();
    edge.action = visualization_msgs::Marker::ADD;
    edge.id = 0;
    edge.type = visualization_msgs::Marker::LINE_STRIP;
    edge.scale.x = 0.01;
    edge.scale.y = 0.01;
    edge.scale.z = 0.01;
    edge.color.a = 1.0;
    edge.color.r = 0.0;
    edge.color.g = 1.0;
    edge.color.b = 0.0;

    int id = 0;
    for (size_t i = 0; i < graph_nodes.size(); ++i) {
        marker.id = id;
        marker.pose.position.x = graph_nodes[i](0);
        marker.pose.position.y = graph_nodes[i](1);
        marker_array.markers.push_back(visualization_msgs::Marker(marker));
        id++;
    }

    for (size_t i = 0; i < graph_edges.size(); ++i) {
        edge.points.clear();
        geometry_msgs::Point p;
        p.x = graph_edges[i].first(0);
        p.y = graph_edges[i].first(1);
        edge.points.push_back(p);
        p.x = graph_edges[i].second(0);
        p.y = graph_edges[i].second(1);
        edge.points.push_back(p);
        edge.id = id;
        marker_array.markers.push_back(visualization_msgs::Marker(edge));
        id++;
    }

    constraint_list_pub_.publish(marker_array);
}




void fusionSlamRos::publishCurLandmarks(const ros::Time& stamp)
{
    visualization_msgs::MarkerArray markers;
    std::vector<Eigen::Vector2f> landmarks = map_builder_.getCurLandmarks();
    for (size_t i = 0; i < landmarks.size(); ++i) {
        visualization_msgs::Marker marker;
        marker.header.frame_id = base_frame_;
        marker.header.stamp = stamp;
        marker.ns = "cur_reflective_marker";
        marker.id = i;
        marker.type = visualization_msgs::Marker::SPHERE;
        marker.action = visualization_msgs::Marker::ADD;
        marker.pose.position.x = landmarks[i].x();
        marker.pose.position.y = landmarks[i].y();
        marker.pose.position.z = 1;
        marker.pose.orientation.x = 0.;
        marker.pose.orientation.y = 0.;
        marker.pose.orientation.z = 0;
        marker.pose.orientation.w = 1;
        marker.scale.x = 1.0;
        marker.scale.y = 1.0;
        marker.scale.z = 0.1;
        marker.color.a = 0.4; // Don't forget to set the alpha!
        marker.color.r = 0.0;
        marker.color.g = 1.0;
        marker.color.b = 0.0;
        marker.lifetime = ros::Duration(0.1);


        markers.markers.push_back(marker);
    }
    cur_landmarks_pub_.publish(markers);

}


void fusionSlamRos::publishAllLandmarks()
{
    visualization_msgs::MarkerArray markers;
    std::vector<std::shared_ptr<RLML::LandMark>> landmarks = map_builder_.getLandmarks();

    for (size_t i = 0; i < landmarks.size(); ++i) {
        visualization_msgs::Marker marker;
        marker.header.frame_id = "map";
        marker.header.stamp = ros::Time::now();
        marker.ns = "reflective_marker";
        marker.id = i;
        marker.type = visualization_msgs::Marker::CUBE;
        marker.action = visualization_msgs::Marker::ADD;
        marker.pose.position.x = landmarks[i]->getMapPos().x();
        marker.pose.position.y = landmarks[i]->getMapPos().y();
        marker.pose.position.z = 0;
        marker.pose.orientation.x = 0.;
        marker.pose.orientation.y = 0.;
        marker.pose.orientation.z = 0;
        marker.pose.orientation.w = 1;
        marker.scale.x = 1.0;
        marker.scale.y = 1.0;
        marker.scale.z = 0.1;
        marker.color.a = 0.5; // Don't forget to set the alpha!
        marker.color.r = 0.0;
        marker.color.g = 0.0;
        marker.color.b = 1.0;
        marker.lifetime = ros::Duration(1./publish_freq_ + 2.);


        markers.markers.push_back(marker);
    }


    landmarks_pub_.publish(markers);

}

void fusionSlamRos::publishLoop()
{
    ros::Rate rate(publish_freq_);

    while (ros::ok()) {
        if(initialize_) {
            if(mapping_)
            {
                publishConstraintList();
                publishProbabilityGridMap();
                publishPath();

            }

            publishCorrelativeGrid();


        }
        publishAllLandmarks();
        rate.sleep();
    }
}


void fusionSlamRos::saveMap()
{
    map_builder_.updateLandmarks();
    map_builder_.doPoseAdjustment();
    map_builder_.saveScans(save_path_, lidar_to_odom_);
    saveOccupancyGridMap();
}



void fusionSlamRos::updateMap()
{
    map_builder_.updateLandmarks();
    map_builder_.doPoseAdjustment();
}


void fusionSlamRos::saveOccupancyGridMap()
{
    std::shared_ptr<RLML::OccupancyGridMap> map = map_builder_.getOccupancyGridMap();

    nav_msgs::OccupancyGrid map_msg;
    map_msg.info.resolution = map->getResolution();
    map_msg.info.width = map->getSizeX();
    map_msg.info.height = map->getSizeY();
    map_msg.data.resize(map_msg.info.width * map_msg.info.height, -1);

    for(size_t i = 0; i < map_msg.data.size(); ++i) {
        if(map->isFree(i)) {
            map_msg.data[i] = 0;
        }
        else if(map->isOccupied(i)) {
            map_msg.data[i] = 100;
        }
    }

    std::string mapdatafile = save_path_ + "map.pgm";
    ROS_INFO("Writing map occupancy data to %s", mapdatafile.c_str());
    FILE* out = fopen(mapdatafile.c_str(), "w");
    if (!out)
    {
        ROS_ERROR("Couldn't save map file to %s", mapdatafile.c_str());
        return;
    }

    fprintf(out, "P5\n# CREATOR: Map_generator.cpp %.3f m/pix\n%d %d\n255\n",
            map_msg.info.resolution, map_msg.info.width, map_msg.info.height);
    for(unsigned int y = 0; y < map_msg.info.height; y++) {
        for(unsigned int x = 0; x < map_msg.info.width; x++) {
            unsigned int i = x + (map_msg.info.height - y - 1) * map_msg.info.width;
            if (map_msg.data[i] == 0) { //occ [0,0.1)
                fputc(254, out);
            } else if (map_msg.data[i] == +100) { //occ (0.65,1]
                fputc(000, out);
            } else { //occ [0.1,0.65]
                fputc(205, out);
            }
        }
    }

    fclose(out);


    std::string mapmetadatafile = save_path_ + "map.yaml";
    ROS_INFO("Writing map occupancy data to %s", mapmetadatafile.c_str());
    FILE* yaml = fopen(mapmetadatafile.c_str(), "w");


    /*
resolution: 0.100000
origin: [0.000000, 0.000000, 0.000000]
#
negate: 0
occupied_thresh: 0.65
free_thresh: 0.196

     */
    Eigen::Vector2f origin = map->getOrigin();

    fprintf(yaml, "image: %s\nresolution: %f\norigin: [%f, %f, %f]\nnegate: 0\noccupied_thresh: 0.65\nfree_thresh: 0.196\n\n",
            mapdatafile.c_str(), map->getResolution(), origin.x(), origin.y(), 0.0  );

    fclose(yaml);

    ROS_INFO("Done\n");

}


template <typename Derived>
void fusionSlamRos::loadCovariance(const ros::NodeHandle& nh, const std::string& cov_name, Eigen::DenseBase<Derived>& covariance)
{
    XmlRpc::XmlRpcValue covariance_config;
    nh.getParam(cov_name, covariance_config);
//    ROS_ASSERT(covariance_config.getType() == XmlRpc::XmlRpcValue::TypeArray);
    int mat_rows = covariance.rows();
    int mat_cols = covariance.cols();
    for (int i = 0; i < mat_rows; i++)
    {
        for (int j = 0; j < mat_cols; j++)
        {
            try
            {
                // These matrices can cause problems if all the types
                // aren't specified with decimal points. Handle that
                // using string streams.
                std::ostringstream ostr;
                ostr << covariance_config[mat_rows * i + j];
                std::istringstream istr(ostr.str());
                istr >> covariance(i, j);
            }
            catch(XmlRpc::XmlRpcException &e)
            {
                throw e;
            }
        }
    }
}

