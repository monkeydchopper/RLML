#include "map_builder.h"




namespace RLML
{
typedef g2o::BlockSolver<g2o::BlockSolverTraits<-1, -1>> SlamBlockSolver;
typedef g2o::LinearSolverCholmod<SlamBlockSolver::PoseMatrixType> SlamLinearSolver;

MapBuilder::MapBuilder() :
    occupancy_grid_map_resolution_(0.05f), min_update_distance_(0.15f),
    min_update_orientation_(degToRad(5.0f)), scan_buffer_size_(30),
    loop_scan_search_distance_(10.0f), loop_closure_min_response_(0.65f), loop_closure_xy_variance_threshold_(0.01f),
    loop_closure_angle_variance_threshold_(0.1f) , loop_match_min_chain_size_(5), loop_closure_xy_search_range_(5.0f),
    loop_closure_angle_search_range_(degToRad(20.f)), loop_closure_grid_resolution_(0.05f),
    loop_closure_coarse_xy_search_resolution_(0.1f), loop_closure_fine_xy_search_range_(0.2f),
    loop_closure_coarse_angle_search_resolution_(degToRad(2.f)),
    loop_closure_fine_angle_search_range_(degToRad(2.f)),
    loop_closure_fine_angle_search_resolution_(degToRad(0.2f)),
    optimize_every_n_constraint_(20), gn_scan_matcher_grid_resolution_(0.05f), use_correlative_scan_matcher_(false),
    csm_xy_search_range_(0.1f), csm_angle_search_range_(degToRad(10.f)),
    csm_grid_resolution_(0.02f), csm_xy_search_resolution_(0.02f),
    csm_angle_search_resolution_(degToRad(0.2f)), search_space_standard_deviation_(0.03f),
    new_loop_constraint_count_(0), loop_constraint_count_(0), got_first_scan_(false), got_first_odom_(false), edge_count_(0), keyframe_select_cnt_(0),
    prev_cost_(0)
{

    sigma_t_minus_1_ <<
    0.05 * 0.05, 0, 0,
    0, 0.05 * 0.05, 0,
    0, 0, 0.03 * 0.03;
    sigma_t_ = sigma_t_minus_1_;
    sigma_u_ = Eigen::Matrix2d::Zero();
    sigma_model_ <<
     0.05 * 0.05, 0, 0,
    0, 0.05 * 0.05, 0,
    0, 0, 0.01 * 0.01;

    sigma_reflector_ <<
    0.025 * 0.025, 0,
    0, 0.025 * 0.025;

    sigma_scan_ = 0.5 * 0.5;

    sigma_k_minus_1_  = sigma_t_minus_1_;
    sigma_relocalize_ <<
    0.03 * 0.03, 0, 0,
    0, 0.03 * 0.03, 0,
    0, 0, 0.01 * 0.01;



}

void MapBuilder::initialize()
{
    std::cout << "occupancy_grid_map_resolution : " << occupancy_grid_map_resolution_ << std::endl;
    std::cout << "min_update_distance : " << min_update_distance_ << std::endl;
    std::cout << "min_update_orientation : " << min_update_orientation_ << std::endl;
    std::cout << "scan_buffer_size : " << scan_buffer_size_ << std::endl;
    std::cout << "loop_scan_search_distance : " << loop_scan_search_distance_ << std::endl;
    std::cout << "loop_match_min_chain_size : " << loop_match_min_chain_size_ << std::endl;
    std::cout << "loop_closure_min_response : " << loop_closure_min_response_ << std::endl;
    std::cout << "loop_closure_xy_variance_threshold : " << loop_closure_xy_variance_threshold_ << std::endl;
    std::cout << "loop_closure_angle_variance_threshold : " << loop_closure_angle_variance_threshold_ << std::endl;
    std::cout << "optimize_every_n_constraint : " << optimize_every_n_constraint_ << std::endl;

    std::cout << "loop_closure_xy_search_range : " << loop_closure_xy_search_range_ << std::endl;
    std::cout << "loop_closure_angle_search_range : " << loop_closure_angle_search_range_ << std::endl;
    std::cout << "loop_closure_grid_resolution : " << loop_closure_grid_resolution_ << std::endl;
    std::cout << "loop_closure_coarse_xy_search_resolution : " << loop_closure_coarse_xy_search_resolution_ << std::endl;
    std::cout << "loop_closure_fine_xy_search_range : " << loop_closure_fine_xy_search_range_ << std::endl;
    std::cout << "loop_closure_coarse_angle_search_resolution : " << loop_closure_coarse_angle_search_resolution_ << std::endl;
    std::cout << "loop_closure_fine_angle_search_range : " << loop_closure_fine_angle_search_range_ << std::endl;
    std::cout << "loop_closure_fine_angle_search_resolution : " << loop_closure_fine_angle_search_resolution_ << std::endl;

    gauss_newten_scan_matcher_.setCorrelativeGrid(gn_scan_matcher_grid_resolution_, search_space_standard_deviation_);

    loop_closure_scan_matcher_.setCoarseXYSearchRange(loop_closure_xy_search_range_);
    loop_closure_scan_matcher_.setCoarseXYSearchResolution(loop_closure_coarse_xy_search_resolution_);
    loop_closure_scan_matcher_.setFineXYSearchRange(loop_closure_fine_xy_search_range_);
    loop_closure_scan_matcher_.setFineXYSearchResolution(loop_closure_grid_resolution_);
    loop_closure_scan_matcher_.setCoarseAngleSearchRange(loop_closure_angle_search_range_);
    loop_closure_scan_matcher_.setCoarseAngleSearchResolution(loop_closure_coarse_angle_search_resolution_);
    loop_closure_scan_matcher_.setFineAngleSearchRange(loop_closure_fine_angle_search_range_);
    loop_closure_scan_matcher_.setFineAngleSearchResolution(loop_closure_fine_angle_search_resolution_);
    loop_closure_scan_matcher_.setCorrelativeGrid(loop_closure_grid_resolution_, search_space_standard_deviation_);

    if(use_correlative_scan_matcher_) {
        std::cout << "csm_xy_search_range : " << csm_xy_search_range_ << std::endl;
        std::cout << "csm_angle_search_range : " << csm_angle_search_range_ << std::endl;
        std::cout << "csm_grid_resolution : " << csm_grid_resolution_ << std::endl;
        std::cout << "csm_xy_search_resolution : " << csm_xy_search_resolution_ << std::endl;
        std::cout << "csm_angle_search_resolution : " << csm_angle_search_resolution_ << std::endl;

        correlative_scan_matcher_.setFineXYSearchRange(csm_xy_search_range_);
        correlative_scan_matcher_.setFineXYSearchResolution(csm_grid_resolution_);
        correlative_scan_matcher_.setFineAngleSearchRange(csm_angle_search_range_);
        correlative_scan_matcher_.setFineAngleSearchResolution(csm_angle_search_resolution_);
        correlative_scan_matcher_.setCorrelativeGrid(csm_grid_resolution_, search_space_standard_deviation_);
    }
}


bool MapBuilder::getClosestScans(const std::shared_ptr<LaserScan>& laser_scan, std::vector<std::shared_ptr<LaserScan>>& map_scans)
{
    // Vector to store pairs of distance and corresponding LaserScan
    std::vector<std::pair<float, std::shared_ptr<LaserScan>>> dis_scan_chain;

    int n = scans_.size();
    for (int i = 0; i < n; ++i) {
        // Compute the Euclidean distance between the input scan and the current scan in the map
        float distance = (laser_scan->getPose().head<2>() - scans_[i]->getPose().head<2>()).norm();

        // Add the scan to the list if the distance is less than 30
        if (distance < 30) {
            dis_scan_chain.emplace_back(distance, scans_[i]);
        }
    }

    // If no nearby scans are found, print a message and return false
    if (dis_scan_chain.empty()) {
        std::cout << "No adjacent map to relocalize" << std::endl;
        return false;
    }

    // Sort the list of scans by distance in ascending order
    std::sort(dis_scan_chain.begin(), dis_scan_chain.end(),
              [](const std::pair<float, std::shared_ptr<LaserScan>>& left,
                 const std::pair<float, std::shared_ptr<LaserScan>>& right) {
                  return left.first < right.first;
              });

    // Determine the number of closest scans to return (limited by scan_buffer_size_)
    int buffer_size = std::min(static_cast<int>(dis_scan_chain.size()), static_cast<int>(scan_buffer_size_));

    // Collect the closest scans into the output vector
    for (int i = 0; i < buffer_size; ++i) {
        map_scans.emplace_back(dis_scan_chain[i].second);
    }

    return true;
}



bool MapBuilder::checkPose(const std::shared_ptr<LaserScan>& laser_scan)
{
    std::shared_ptr<LaserScan>& last_laser_scan = scans_.back();

    Eigen::Vector3f delta_pose = last_laser_scan->getPose() - laser_scan->getPose();

    float update_distance = delta_pose[0] * delta_pose[0] + delta_pose[1] * delta_pose[1];
    float update_angle = fabs(normalizeAngle(delta_pose[2]));

    if(update_distance >= min_update_distance_ * min_update_distance_ || update_angle >= min_update_orientation_) {
        return true;
    }
    else {
        return false;
    }
}

void MapBuilder::addRunningScan(std::shared_ptr<LaserScan> laser_scan)
{
    running_scan_.push_back(laser_scan);

    while (running_scan_.size() > scan_buffer_size_) {
        running_scan_.erase(running_scan_.begin());
    }
}

void MapBuilder::addVertex(g2o::SparseOptimizer& optimizer, const std::shared_ptr<LaserScan>& scan)
{
    Eigen::Vector3f pose = scan->getPose();
    g2o::VertexSE2* vertex = new g2o::VertexSE2;
    vertex->setEstimate(g2o::SE2(pose[0], pose[1], pose[2]));
    vertex->setId(scan->getId());
    optimizer.addVertex(vertex);
}


void MapBuilder::addVertex(g2o::SparseOptimizer& optimizer, const std::shared_ptr<LandMark>& landmark, const int offset_id)
{
    g2o::VertexPointXY* vertex = new g2o::VertexPointXY;
    vertex->setId(offset_id + landmark->getId());
    vertex->setEstimate(landmark->getMapPos().cast<double>());
    optimizer.addVertex(vertex);
}


void MapBuilder::addEdge(g2o::SparseOptimizer& optimizer,
        std::shared_ptr<LaserScan> source_scan, const Eigen::Vector3f& source_pose,
        std::shared_ptr<LaserScan> target_scan, const Eigen::Vector3f& target_pose,
        const Eigen::Matrix3d& information)
{
    Eigen::AngleAxisf rotation(-source_pose[2], Eigen::Vector3f(0, 0, 1));
    Eigen::Vector3f delta_pose = rotation * (target_pose - source_pose);

    g2o::EdgeSE2* edge = new g2o::EdgeSE2;

    int source_id = source_scan->getId();
    int target_id = target_scan->getId();

    edge->vertices()[0] = optimizer.vertex(source_id);
    edge->vertices()[1] = optimizer.vertex(target_id);

    g2o::SE2 measurement(delta_pose[0], delta_pose[1], delta_pose[2]);
//    std::cout << "add pose edge " << edge_count_ << std::endl;
    edge->setId(edge_count_);
    edge->setMeasurement(measurement);
    edge->setInformation(information);
    edge_count_++;

    optimizer.addEdge(edge);
}


void MapBuilder::addEdge(g2o::SparseOptimizer& optimizer,
        const poseEdge& pose_edge)
{
    g2o::EdgeSE2* edge = new g2o::EdgeSE2;

    int source_id = pose_edge.from_scan_->getId();
    int target_id = pose_edge.to_scan_->getId();

    edge->vertices()[0] = optimizer.vertex(source_id);
    edge->vertices()[1] = optimizer.vertex(target_id);

    g2o::SE2 measurement(pose_edge.mes_[0], pose_edge.mes_[1], pose_edge.mes_[2]);
//    std::cout << "add pose edge " << edge_count_ << std::endl;
    edge->setId(edge_count_);
    edge->setMeasurement(measurement);
    edge->setInformation(pose_edge.info_);
    edge_count_++;

    optimizer.addEdge(edge);
}



void MapBuilder::addEdge(g2o::SparseOptimizer& optimizer,
                         const landmarkEdge& landmark_edge,
                         const int offset_id)
{
    g2o::EdgeSE2PointXY* edge =  new g2o::EdgeSE2PointXY;
    edge->vertices()[0] = optimizer.vertex(landmark_edge.laser_->getId());
    edge->vertices()[1] = optimizer.vertex(offset_id + landmark_edge.landmark_->getId());
//    std::cout << "add landmark edge " << edge_count_ << std::endl;

    edge->setId(edge_count_);
    edge->setMeasurement(landmark_edge.mes_.cast<double>());
    edge->setInformation(landmark_edge.info_);
    edge_count_++;
    optimizer.addEdge(edge);
}

std::shared_ptr<LaserScan> MapBuilder::getClosestScan(const std::shared_ptr<LaserScan>& base_scan,
                                                     const std::vector<std::shared_ptr<LaserScan>>& chain)
{
    float min_distance = std::numeric_limits<float>::max();
    std::shared_ptr<LaserScan> closest_scan;

    for(const std::shared_ptr<LaserScan>& scan : chain) {
        float distance = (scan->getPose().head<2>() - base_scan->getPose().head<2>()).norm();
        if(min_distance > distance) {
            min_distance = distance;
            closest_scan = scan;
        }
    }

    return closest_scan;
}

void MapBuilder::detectLoopClosure(const std::shared_ptr<LaserScan>& scan, const std::vector<Eigen::Vector2f>& reflective_markers)
{
    std::vector<std::shared_ptr<LaserScan>> scan_chain;
    std::vector<std::vector<std::shared_ptr<LaserScan>>> scan_chains;

    int n = scans_.size();
    for(int i = 0; i < n; ++i) {
        float distance = (scan->getPose().head<2>() - scans_[i]->getPose().head<2>()).norm();

        if(distance < loop_scan_search_distance_) {
            if(std::find(running_scan_.begin(), running_scan_.end(), scans_[i]) == running_scan_.end()) {
                scan_chain.push_back(scans_[i]);
            }
            else {
                scan_chain.clear();
            }
        }
        else {
            if(scan_chain.size() > loop_match_min_chain_size_) {
                scan_chains.push_back(scan_chain);
                scan_chain.clear();
            }
            else {
                scan_chain.clear();
            }
        }
    }

    if(scan_chains.empty()) {
        return;
    }

    // construct previous landmarks
    std::vector<std::shared_ptr<LandMark>> previous_landmarks = landmarks_;
    for (size_t i = 0; i < running_scan_.size(); ++i) {
        std::map<std::shared_ptr<LandMark>, Eigen::Vector2f> landmarks_associations = running_scan_[i]->getLandmarks();
        for(auto const& asso:landmarks_associations)
        {
            auto it = std::find(previous_landmarks.begin(), previous_landmarks.end(), asso.first);
            if(it != previous_landmarks.end())
            {
                previous_landmarks.erase(it);
            }
        }
    }

    std::vector<Eigen::Vector3f> pose_candidates;
    for(const std::vector<std::shared_ptr<LaserScan>>& chain : scan_chains) {
        Eigen::Vector3f csm_pose;
        Eigen::Matrix3f covariance;
        auto t1 = std::chrono::steady_clock::now();
        float response = loop_closure_scan_matcher_.multiResMatchScan(scan, chain, csm_pose, covariance);
        auto t2 = std::chrono::steady_clock::now();
        auto delta_t = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);

        if(response < loop_closure_min_response_) {
            continue;
        }

        if(covariance(0, 0) > loop_closure_xy_variance_threshold_ ||
           covariance(1, 1) > loop_closure_xy_variance_threshold_) {
            std::cout << "Reject loop closure with large covariance : " << std::endl;
            std::cout << covariance << std::endl;
            continue;
        }

        Eigen::Vector3f gnsm_pose = gauss_newten_scan_matcher_.matchScanToMap(csm_pose, loop_closure_scan_matcher_.getCorrelativeGrid(),
                                                                  scan->getRawPointCloud());

        float delta_trans = (gnsm_pose.head<2>() - csm_pose.head<2>()).norm();
        float delta_angle = normalizeAngle(gnsm_pose[2] - csm_pose[2]);
        if(fabs(delta_trans) > 0.15 && fabs(delta_angle) > 0.1) {
            std::cout <<  "Loop closure gauss newten scan matcher result got large jump with delta_trans = " << delta_trans
                      << " delta_angle = " << delta_angle << std::endl;
        }



        new_loop_constraint_count_++;
        loop_constraint_count_++;
        std::cout << "Find loop closure with response " << response << ". Cost time " << delta_t.count() * 1000.0 << "ms." << std::endl;
        std::shared_ptr<LaserScan> closest_scan = getClosestScan(scan, chain);

        // push loop constraint
        poseEdge pe;
        pe.from_scan_ = scan;
        pe.to_scan_ = closest_scan;
//        Eigen::AngleAxisf rotation(-csm_pose[2], Eigen::Vector3f(0, 0, 1));
//        Eigen::Vector3f delta_pose = rotation * (closest_scan->getPose() - csm_pose);
        Eigen::AngleAxisf rotation(-gnsm_pose[2], Eigen::Vector3f(0, 0, 1));
        Eigen::Vector3f delta_pose = rotation * (closest_scan->getPose() - gnsm_pose);
        pe.mes_ = delta_pose;
        Eigen::Matrix3d info_mat;

        info_mat = sigma_relocalize_.inverse();
        pe.info_ = info_mat * response;
        std::cout << "loop info " << std::endl << pe.info_ << std::endl;
        poseEdges_.push_back(pe);
        pose_candidates.push_back(csm_pose);
    }

    float avg_x = 0.f, avg_y = 0.f, sin_theta = 0.f, cos_theta = 0.f;
    for (size_t j = 0; j < pose_candidates.size(); ++j) {
        avg_x += pose_candidates[j][0];
        avg_y += pose_candidates[j][1];
        sin_theta += std::sin(pose_candidates[j][2]);
        cos_theta += std::cos(pose_candidates[j][2]);
    }

    // find associations in the map
    Eigen::Vector3f relocalize_pose(avg_x/pose_candidates.size(), avg_y/pose_candidates.size(), std::atan2(sin_theta, cos_theta));
    std::map<int, std::shared_ptr<LandMark>> associations_all, associations_ransac;
    associateLandmarksRANSAC(relocalize_pose, reflective_markers, previous_landmarks, associations_all, associations_ransac);

    if(associations_ransac.size() > 0)
    {
        std::cout << "loop find landmarks size" << associations_ransac.size() << std::endl;
    }

    // add associations into scan and landmark
    for(auto const& asso:associations_ransac)
    {
        asso.second->addObservations(scan, reflective_markers[asso.first]);
        scan->addLandmark(reflective_markers[asso.first],asso.second);
    }

}



bool MapBuilder::relocalize(const std::shared_ptr<LaserScan>& scan, const std::vector<Eigen::Vector2f>& reflective_markers, const bool use_reflector_relocalize)
{
    if(use_reflector_relocalize)
    {
        // use three distances
        std::vector<float> distances;
        if(reflective_markers.size() < 3)
        {
            std::cout << "reflective markers not enough" << std::endl;
            std::cerr << "reflective markers not enough" << std::endl;
            return false;
        }
        std::vector<Eigen::Vector2f> selected_markers;
        std::sample(reflective_markers.begin(), reflective_markers.end(), std::back_inserter(selected_markers),
                    3, std::mt19937{std::random_device{}()});
        // dis between markers
        distances.push_back((selected_markers[0] - selected_markers[1]).norm());
        distances.push_back((selected_markers[0] - selected_markers[2]).norm());
        distances.push_back((selected_markers[1] - selected_markers[2]).norm());

        // check for inputs
        if(std::abs(distances[0] - distances[1]) < 0.2f ||
                std::abs(distances[0] - distances[2]) < 0.2f ||
                std::abs(distances[1] - distances[2]) < 0.2f)
        {
            std::cout << "there is ambiguity in landmark relocalization" << std::endl;
            std::cerr << "there is ambiguity in landmark relocalization" << std::endl;
            return false;
        }

        // search for marker0
        std::vector<int> landmark_index(3,-1);
        bool found_relocalize = false;
        for (size_t i = 0; i < landmark_dis_tab_.size(); ++i) {
            bool found_marker1 = false, found_marker2 = false;
            for(const auto& landmark_dis:landmark_dis_tab_[i])
            {
                if(std::abs(landmark_dis.second - distances[0]) < 0.15f)
                {
                    found_marker1 = true;
                    landmark_index[1] = landmark_dis.first;
                }
                if(std::abs(landmark_dis.second - distances[1]) < 0.15f)
                {
                    found_marker2 = true;
                    landmark_index[2] = landmark_dis.first;
                }
            }
            if(found_marker1 && found_marker2 &&
            ((landmarks_[landmark_index[1]]->getMapPos() - landmarks_[landmark_index[2]]->getMapPos()).norm() - distances[2]) < 0.15f)
            {
                landmark_index[0] = i;
                found_relocalize = true;
                break;
            }
        }

        if(!found_relocalize)
        {
            std::cout << "cannot use reflector to find relocalization in this frame" << std::endl;
            std::cerr << "cannot use reflector to find relocalization in this frame" << std::endl;
            return false;
        }

        std::vector<int> flags;
        std::vector<Eigen::Vector2f> mes, est;
        for (size_t i = 0; i < landmark_index.size(); ++i) {
            mes.push_back(landmarks_[landmark_index[i]]->getMapPos());
            est.push_back(selected_markers[i]);
        }
        Eigen::Matrix<float,6,4> A;
        Eigen::Matrix<float,6,1> b;
        A <<
          est[0].x(), -est[0].y(), 1, 0,
                est[0].y(),  est[0].x(), 0, 1,
                est[1].x(), -est[1].y(), 1, 0,
                est[1].y(),  est[1].x(), 0, 1,
                est[2].x(), -est[2].y(), 1, 0,
                est[2].y(),  est[2].x(), 0, 1;
        if((A.transpose() * A).determinant() < 1e-2)
        {
            std::cout << "reflector estimate pose error" << std::endl;
            return false;
        }
        b << mes[0].x(), mes[0].y(),mes[1].x(),mes[1].y(),mes[2].x(),mes[2].y();
//        std::cout << "A = " << A << std::endl;
//        std::cout << "A inverse = " << (A.transpose() * A).inverse() * A.transpose() << std::endl;
//            std::cout << "b = " << b << std::endl;

        Eigen::Vector4f X = (A.transpose() * A).inverse() * A.transpose() * b;
        std::cout << "initial pose x = " << X(2) << " y = " << X(3) << " theta = " << atan2(X(1), X(0)) << std::endl;

        Eigen::Vector3f initial_pose(X(2),X(3),atan2(X(1), X(0)));
        Eigen::Vector3f pose_diff = initial_pose - scan->getPose();
        scan->setPose(initial_pose);

        std::cout << "pose jump from intial pose setting x = " << pose_diff(0) << " y = " << pose_diff(1) << " theta = " << pose_diff(2) << std::endl;
    }


    std::vector<std::shared_ptr<LaserScan>> map_scans;
    if(!getClosestScans(scan, map_scans))
    {
        std::cout << "cannot get closest scans" << std::endl;
        return false;
    }
    Eigen::Vector3f csm_pose;
    Eigen::Matrix3f covariance;
    auto t1 = std::chrono::steady_clock::now();
    float response = loop_closure_scan_matcher_.multiResMatchScan(scan, map_scans, csm_pose, covariance);
    auto t2 = std::chrono::steady_clock::now();
    auto delta_t = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);

    if(response < relocalize_min_response_) {
        std::cout << "reject relocalize with small response " << std::endl << response << std::endl;
        std::cerr << "reject relocalize with small response " << std::endl << response << std::endl;

        return false;
    }
    else
    {
        std::cout << "relocalize with response " << std::endl << response << std::endl;
    }

//    if(covariance(0, 0) > loop_closure_xy_variance_threshold_ ||
//       covariance(1, 1) > loop_closure_xy_variance_threshold_) {
//            std::cout << "reject relocalize with large covariance" << std::endl;
//            std::cout << covariance << std::endl;
//        return false;
//    } else
//    {
        std::cout << "relocalize with covariance" << std::endl;
        std::cout << covariance << std::endl;
//    }

    Eigen::Vector3f gnsm_pose = gauss_newten_scan_matcher_.matchScanToMap(csm_pose, loop_closure_scan_matcher_.getCorrelativeGrid(),
                                                                          scan->getRawPointCloud());

    float delta_trans = (gnsm_pose.head<2>() - csm_pose.head<2>()).norm();
    float delta_angle = normalizeAngle(gnsm_pose[2] - csm_pose[2]);
    if(fabs(delta_trans) > 0.15 && fabs(delta_angle) > 0.1) {
        std::cout <<  "relocalize gauss newten scan matcher result got large jump with delta_trans = " << delta_trans
                  << " delta_angle = " << delta_angle << std::endl;
    }

    std::cout <<   "reflector relocalize pose: " << scan->getPose().x() << " " << scan->getPose().y() << " "
    << scan->getPose().z() << std::endl;
    std::cout <<   "matching relocalize pose: " << gnsm_pose.x() << " " << gnsm_pose.y() << " "
              << gnsm_pose.z() << std::endl;

    // find associations in the map
    std::map<int, std::shared_ptr<LandMark>> associations_all, associations_ransac;
    std::vector<std::shared_ptr<LandMark>> previous_landmarks = landmarks_;
    associateLandmarksRANSAC(gnsm_pose, reflective_markers, previous_landmarks, associations_all, associations_ransac);

    if(associations_ransac.empty())
    {
        std::cout << "cannot find reflector association with this initial pose, relocalize fail" << std::endl;
        std::cerr << "cannot find reflector association with this initial pose, relocalize fail" << std::endl;
        return false;
    }

    Eigen::Matrix3d sigma_update = sigma_t_;
    gnsm_pose = gauss_newten_scan_matcher_.matchScanToMap( gnsm_pose,  loop_closure_scan_matcher_.getCorrelativeGrid(), scan->getRawPointCloud(), reflective_markers, associations_ransac,
                                                                     sigma_t_, sigma_reflector_, sigma_scan_,sigma_update);

    scan->setPose(gnsm_pose);
    return true;
}


bool MapBuilder::localRelocalize(const std::shared_ptr<LaserScan>& scan)
{

    if(running_scan_.empty())
    {
        std::cout << "cannot local relocalize, not enough running scan" << std::endl;
        return false;
    }


    scan->setPose(last_scan_->getPose());
    std::cout << "pose before relocalize: " << std::endl << last_scan_->getPose() << std::endl;

    std::vector<std::shared_ptr<LaserScan>> map_scans;
    if(!getClosestScans(scan, map_scans))
    {
        std::cout << "cannot get closest scans" << std::endl;
        return false;
    }    map_scans.insert(map_scans.end(), running_scan_.begin(),running_scan_.end());
    Eigen::Vector3f csm_pose;
    Eigen::Matrix3f covariance;
    auto t1 = std::chrono::steady_clock::now();
    float response = loop_closure_scan_matcher_.multiResMatchScan(scan, map_scans, csm_pose, covariance);
    auto t2 = std::chrono::steady_clock::now();
    auto delta_t = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);

    if(response < relocalize_min_response_) {
        std::cout << "reject relocalize with small response " << std::endl << response << std::endl;
        return false;
    }
    else
    {
        std::cout << "relocalize with response " << std::endl << response << std::endl;
    }

//    if(covariance(0, 0) > loop_closure_xy_variance_threshold_ ||
//       covariance(1, 1) > loop_closure_xy_variance_threshold_) {
//        std::cout << "reject relocalize with large covariance" << std::endl;
////            std::cout << covariance << std::endl;
//        return false;
//    }

    Eigen::Vector3f gnsm_pose = gauss_newten_scan_matcher_.matchScanToMap(csm_pose, loop_closure_scan_matcher_.getCorrelativeGrid(),
                                                                          scan->getRawPointCloud());

    float delta_trans = (gnsm_pose.head<2>() - csm_pose.head<2>()).norm();
    float delta_angle = normalizeAngle(gnsm_pose[2] - csm_pose[2]);
    if(fabs(delta_trans) > 0.15 && fabs(delta_angle) > 0.1) {
        std::cout <<  "Loop closure gauss newten scan matcher result got large jump with delta_trans = " << delta_trans
                  << " delta_angle = " << delta_angle << std::endl;
    }

//    std::cout << "local relocalize succeed" << std::endl;

    scan->setPose(gnsm_pose);
    last_scan_ = scan;
    std::cout << "pose after relocalize: " << std::endl << gnsm_pose << std::endl;

    return true;
}

bool MapBuilder::doPoseAdjustment()
{
    edge_count_ = 0;
    if(poseEdges_.empty())
    {
        std::cout << "nothing to optimize" << std::endl;
        return true;
    }
    g2o::SparseOptimizer optimizer;
//    SlamBlockSolver::LinearSolverType* linear_solver = new SlamLinearSolver;
//    SlamBlockSolver* solver_ptr = new SlamBlockSolver(linear_solver);
//    g2o::OptimizationAlgorithmLevenberg* solver = new g2o::OptimizationAlgorithmLevenberg(solver_ptr); // L-M
//    g2o::SparseOptimizer* optimizer;
    auto linearSolver = g2o::make_unique<SlamLinearSolver>();
    linearSolver->setBlockOrdering(false);
    g2o::OptimizationAlgorithmLevenberg* solver = new g2o::OptimizationAlgorithmLevenberg(g2o::make_unique<SlamBlockSolver>(std::move(linearSolver)));
    optimizer.setAlgorithm(solver);
    optimizer.setVerbose(false);


    std::vector<landmarkEdge> landmarkEdges;

    // add all the pose vertices
    for (size_t i = 0; i < scans_.size(); ++i) {
        addVertex(optimizer, scans_[i]);
        std::map<std::shared_ptr<LandMark>, Eigen::Vector2f> landmarks = scans_[i]->getLandmarks();
        for(auto const& ld:landmarks)
        {
            landmarkEdge le;
            le.laser_ = scans_[i];
            le.landmark_ = ld.first;
            le.mes_ = ld.second;
            le.info_ = sigma_reflector_.inverse();
            landmarkEdges.push_back(le);
        }
    }

    // add constraints between poses
    for (size_t i = 0; i < poseEdges_.size(); ++i) {
        addEdge(optimizer, poseEdges_[i]);
    }

    // add landmarks vertices
    for (size_t i = 0; i < landmarks_.size(); ++i) {
        addVertex(optimizer, landmarks_[i], scans_.size());
    }


    // add constraints between pose and landmark
    for (size_t i = 0; i < landmarkEdges.size(); ++i) {
        addEdge(optimizer, landmarkEdges[i], scans_.size());
    }

    g2o::OptimizableGraph::Vertex* v = optimizer.vertex(0);
    v->setFixed(true);

    optimizer.initializeOptimization();

    auto t1 = std::chrono::steady_clock::now();
    int iter = optimizer.optimize(100);
    auto t2 = std::chrono::steady_clock::now();
    auto delta_t = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);

    if (iter > 0) {
        std::cout << "global optimization finished after " << iter << " iterations. Cost time " <<
                     delta_t.count() * 1000.0 << "ms." << std::endl;
    }
    else {
        std::cout << "global optimization failed, result might be invalid!" << std::endl;
        return false;
    }

    g2o::SparseOptimizer::VertexContainer nodes = optimizer.activeVertices();


    assert(nodes.size() == (scans_.size() + landmarks_.size()));
    for (size_t i = 0; i < nodes.size(); ++i) {
        if(i < scans_.size())
        {
            // update poses
            double estimate[3];
            if(nodes[i]->getEstimateData(estimate))
            {
                Eigen::Vector3f pose(static_cast<float>(estimate[0]),
                                     static_cast<float>(estimate[1]),
                                     static_cast<float>(estimate[2]));
                scans_[i]->setPose(pose);
                scans_[i]->transformPointCloud();
            }
            else {
                std::cout <<"Could not get estimated pose from Optimizer!" << std::endl;
            }
        } else
        {
            // update landmarks
            double estimate[2];
            std::shared_ptr<LandMark> landmark;
//            bool got_data = nodes[i]->getEstimateData(estimate);
//            bool got_landmark = getLandmarkByID(nodes[i]->id() - scans_.size(), landmark);
            if(nodes[i]->getEstimateData(estimate) && getLandmarkByID(nodes[i]->id() - scans_.size(), landmark))
//            if(got_data && got_landmark)
            {
                Eigen::Vector2f pos(estimate[0], estimate[1]);
                landmark->setMapPos(pos);
            } else
            {
                std::cout <<"Could not get estimated landmark from Optimizer!" << std::endl;
            }
        }
    }


    return true;


}


bool MapBuilder::doLocalPoseAdjustment(const std::vector<std::shared_ptr<LandMark>>& running_landmarks)
{
    edge_count_ = 0;
    g2o::SparseOptimizer optimizer;
//    SlamBlockSolver::LinearSolverType* linear_solver = new SlamLinearSolver;
//    SlamBlockSolver* solver_ptr = new SlamBlockSolver(linear_solver);
//    g2o::OptimizationAlgorithmLevenberg* solver = new g2o::OptimizationAlgorithmLevenberg(solver_ptr); // L-M
    auto linearSolver = g2o::make_unique<SlamLinearSolver>();
    linearSolver->setBlockOrdering(false);
    g2o::OptimizationAlgorithmLevenberg* solver = new g2o::OptimizationAlgorithmLevenberg(g2o::make_unique<SlamBlockSolver>(std::move(linearSolver)));
    optimizer.setAlgorithm(solver);
    optimizer.setVerbose(false);


    std::vector<landmarkEdge> landmarkEdges;
    std::vector<std::shared_ptr<LaserScan>> associate_scans = running_scan_;

    // add all the running scans to vertex
    for (size_t i = 0; i < associate_scans.size(); ++i) {
//        std::cout << "add scan vertex " << associate_scans[i]->getId() << std::endl;
        addVertex(optimizer, associate_scans[i]);
        if(associate_scans[i]->getId() == 0)
        {
            g2o::OptimizableGraph::Vertex* v = optimizer.vertex(0);
            v->setFixed(true);
        }
    }

    // find all the scans associate with running landmarks
    // set scan not in running scans to be fixed
    // add scans to vertices, push landmarkEdges
    for (size_t i = 0; i < running_landmarks.size(); ++i) {

//        std::cout << "i " << i << "add landmark vertex " << running_landmarks[i]->getId() << std::endl;
//        std::cout << "add landmark vertex id " << running_landmarks[i]->getId() + scans_.size() << std::endl;
        addVertex(optimizer, running_landmarks[i], scans_.size());
        std::map<std::shared_ptr<LaserScan>, Eigen::Vector2f> observations = running_landmarks[i]->getObservations();
        for(const auto& obs:observations)
        {
            // not in running scans, set to fix
            if(std::find(associate_scans.begin(), associate_scans.end(), obs.first) == associate_scans.end())
            {
                associate_scans.push_back(obs.first);
//                std::cout << "add landmark scan vertex " << obs.first->getId() << std::endl;

                addVertex(optimizer, obs.first);
                g2o::OptimizableGraph::Vertex* v = optimizer.vertex(obs.first->getId());
                v->setFixed(true);
            }
            landmarkEdge le;
            le.laser_ = obs.first;
            le.landmark_ = running_landmarks[i];
            le.mes_ = obs.second;
            le.info_  = sigma_reflector_.inverse();
            landmarkEdges.push_back(le);
        }

    }


    // find all the scans associate with running scans
    // set scan not in running scans to be fixed
    // add scans to vertices, poseEdge to edges
    for (size_t i = 0; i < poseEdges_.size(); ++i) {
        auto from = std::find(running_scan_.begin(), running_scan_.end(), poseEdges_[i].from_scan_);
        auto to = std::find(running_scan_.begin(), running_scan_.end(), poseEdges_[i].to_scan_);

        // from and to both in running scans
        // as running scans are all in vectices, so just add edge
//        if(from != running_scan_.end() && to != running_scan_.end())
//        {
//            addEdge(optimizer, poseEdges_[i]);
//        }
        // from in running scan, but to not in running scan
        // if to not in associate scans, add it to associate scans and vertices and set to fix
        if(from != running_scan_.end() && to == running_scan_.end())
        {
            if(std::find(associate_scans.begin(), associate_scans.end(), poseEdges_[i].to_scan_) == associate_scans.end())
            {
                associate_scans.push_back(poseEdges_[i].to_scan_);
//                std::cout << "add to scan vertex " << poseEdges_[i].to_scan_->getId() << std::endl;

                addVertex(optimizer, poseEdges_[i].to_scan_);
                g2o::OptimizableGraph::Vertex* v = optimizer.vertex(poseEdges_[i].to_scan_->getId());
                v->setFixed(true);
            }
        }

        // from not in running scan, but to in running scan
        // if to not in associate scans, add it to associate scans and vertices and set to fix
        if(from == running_scan_.end() && to != running_scan_.end())
        {
            if(std::find(associate_scans.begin(), associate_scans.end(), poseEdges_[i].from_scan_) == associate_scans.end())
            {
                associate_scans.push_back(poseEdges_[i].from_scan_);
//                std::cout << "add from scan vertex " << poseEdges_[i].from_scan_->getId() << std::endl;

                addVertex(optimizer, poseEdges_[i].from_scan_);
                g2o::OptimizableGraph::Vertex* v = optimizer.vertex(poseEdges_[i].from_scan_->getId());
                v->setFixed(true);
            }
        }

        addEdge(optimizer, poseEdges_[i]);

    }

    // add constraints between pose and landmark
    for(size_t i = 0; i < landmarkEdges.size(); ++i) {
        addEdge(optimizer, landmarkEdges[i], static_cast<int>(scans_.size()));
    }


    optimizer.initializeOptimization();

    auto t1 = std::chrono::steady_clock::now();
    int iter = optimizer.optimize(100);
    auto t2 = std::chrono::steady_clock::now();
    auto delta_t = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);

    if (iter > 0) {
        std::cout << "local optimization finished after " << iter << " iterations. Cost time " <<
                  delta_t.count() * 1000.0 << "ms." << std::endl;
    }
    else {
        std::cout << "local optimization failed, result might be invalid!" << std::endl;
        return false;
    }

    g2o::SparseOptimizer::VertexContainer nodes = optimizer.activeVertices();

    if(nodes.size() != (associate_scans.size() + running_landmarks.size()))
    {
        std::cout << "node " << nodes.size() << " scan size " << associate_scans.size() << " landmark size " << running_landmarks.size() << std::endl;
        std::cout << "local optimization parameters wrong, thread will shut down!" << std::endl;
    }
    assert(nodes.size() == (associate_scans.size() + running_landmarks.size()));
    for(size_t i = 0; i < nodes.size(); ++i) {
        if(nodes[i]->id() < static_cast<int>(scans_.size()))
        {
            // update poses
            double estimate[3];
            if(nodes[i]->getEstimateData(estimate))
            {
                Eigen::Vector3f pose(static_cast<float>(estimate[0]),
                                     static_cast<float>(estimate[1]),
                                     static_cast<float>(estimate[2]));
                scans_[nodes[i]->id()]->setPose(pose);
                scans_[nodes[i]->id()]->transformPointCloud();
            }
            else {
                std::cout <<"Could not get estimated pose from Optimizer!" << std::endl;
            }
        } else
        {
            // update landmarks
            double estimate[2];
            std::shared_ptr<LandMark> landmark;
//            bool got_data = nodes[i]->getEstimateData(estimate);
//            bool got_landmark = getLandmarkByID(nodes[i]->id() - scans_.size(), landmark);
            if(nodes[i]->getEstimateData(estimate) && getLandmarkByID(nodes[i]->id() - scans_.size(), landmark))
//            if(got_data && got_landmark)
            {
                Eigen::Vector2f pos(estimate[0], estimate[1]);
                landmark->setMapPos(pos);
            } else
            {
                std::cout <<"Could not get estimated landmark from Optimizer!" << std::endl;
            }
        }
    }
//    std::cout <<"local optimization finished" << std::endl;
    return true;


}

void MapBuilder::alignLandmarkToMap()
{
    for (size_t i = 0; i < landmarks_.size(); ++i) {
        std::shared_ptr<LaserScan> ref_scan = landmarks_[i]->getReferenceScan();
        Eigen::Vector2f landmark_in_scan;
        if(landmarks_[i]->getObservation(ref_scan, landmark_in_scan))
        {
            Eigen::Vector2f landmark_in_map = ref_scan->transformLandmark(landmark_in_scan);
            landmarks_[i]->setMapPos(landmark_in_map);
        } else
        {
            std::cout << "ref scan not in observations, that's bad" << std::endl;
        }
    }
}


void MapBuilder::mergeLandmarks(const std::vector<std::shared_ptr<LandMark>>& running_landmarks)
{

    // get rid of landmarks which has few observations
//    for (int i = 0; i < landmarks_.size(); ++i) {
//        std::shared_ptr<LaserScan> ref_scan = landmarks_[i]->getReferenceScan();
//        if(std::find(running_scan_.begin(), running_scan_.end(), ref_scan) == running_scan_.end())
//        {
//            if(landmarks_[i]->getObsNum() <= erase_landmark_nobs_threshold_)
//            {
//                landmarks_[i]->setState(false);
//                std::map<std::shared_ptr<LaserScan>, Eigen::Vector2f> observations = landmarks_[i]->getObservations();
//                for(const auto& obs:observations)
//                {
//                    obs.first->deleteLandmark(landmarks_[i]);
//                }
//            }
//        }
//    }

    // merge obs into old landmarks
    for (size_t i = 0; i < landmarks_.size(); ++i) {
        size_t j = landmarks_.size() - i - 1;
//        std::shared_ptr<LaserScan> ref_scan = landmarks_[j]->getReferenceScan();
        if(std::find(running_landmarks.begin(), running_landmarks.end(), landmarks_[j]) == running_landmarks.end())
        {
            // get rid of landmark which has few observations
            if(landmarks_[j]->getObsNum() <= erase_landmark_nobs_threshold_)
            {
                landmarks_[j]->setState(false);
                continue;
            }

            // landmark has enough observations, try to merge it
            for (size_t k = 0; k < j; ++k) {
                float dis = (landmarks_[k]->getMapPos().head<2>() - landmarks_[j]->getMapPos().head<2>()).norm();
                if(dis < merge_dis_threshold_)
                {
                    // merge landmark j into landmark k
                    std::map<std::shared_ptr<LaserScan>, Eigen::Vector2f> merge_obs = landmarks_[j]->getObservations();
                    for(auto const& ob:merge_obs)
                    {
                        // landmark k add observation
                        landmarks_[k]->addObservations(ob.first, ob.second);
                        // scan change its association
//                        ob.first->replaceLandmark(landmarks_[j], landmarks_[k]);
                        ob.first->addLandmark(ob.second, landmarks_[k]);
                    }
                    landmarks_[j]->setState(false);
                }
            }
        }

    }



    // resize landmarks
    int j = 0;
    for(size_t i = 0; i < landmarks_.size(); ++i) {
        if(landmarks_[i]->getState())
        {
            landmarks_[j++] =  landmarks_[i];
        } else
        {
            std::map<std::shared_ptr<LaserScan>, Eigen::Vector2f> observations = landmarks_[i]->getObservations();
            for(const auto& obs:observations)
            {
                obs.first->deleteLandmark(landmarks_[i]);
            }
        }
    }
    landmarks_.resize(j);
}

void MapBuilder::updateLandmarks()
{
// merge obs into old landmarks
    for (size_t i = 0; i < landmarks_.size(); ++i) {
        size_t j = landmarks_.size() - i - 1;

        // get rid of landmark which has few observations
        if(landmarks_[j]->getObsNum() <= erase_landmark_nobs_threshold_)
        {
            landmarks_[j]->setState(false);
            continue;
        }
        // landmark has enough observations, try to merge it
        for (size_t k = 0; k < j; ++k) {
            float dis = (landmarks_[k]->getMapPos().head<2>() - landmarks_[j]->getMapPos().head<2>()).norm();
            if(dis < merge_dis_threshold_)
            {
                // merge landmark j into landmark k
                std::map<std::shared_ptr<LaserScan>, Eigen::Vector2f> merge_obs = landmarks_[j]->getObservations();
                for(auto const& ob:merge_obs)
                {
                    // landmark k add observation
                    landmarks_[k]->addObservations(ob.first, ob.second);
                    // scan change its association
//                        ob.first->replaceLandmark(landmarks_[j], landmarks_[k]);
                    ob.first->addLandmark(ob.second, landmarks_[k]);
                }
                landmarks_[j]->setState(false);
            }
        }
    }

    // resize landmarks
    size_t j = 0;
    for (size_t i = 0; i < landmarks_.size(); ++i) {
        if(landmarks_[i]->getState())
        {
            landmarks_[j++] =  landmarks_[i];
        } else
        {
            std::map<std::shared_ptr<LaserScan>, Eigen::Vector2f> observations = landmarks_[i]->getObservations();
            for(const auto& obs:observations)
            {
                obs.first->deleteLandmark(landmarks_[i]);
            }
        }
    }
    landmarks_.resize(j);
}

void MapBuilder::discardLandmarks()
{
    // get rid of landmarks which has few observations and merge ladnmarks
    for (size_t i = 0; i < landmarks_.size(); ++i) {
        std::shared_ptr<LaserScan> ref_scan = landmarks_[i]->getReferenceScan();
        if(std::find(running_scan_.begin(), running_scan_.end(), ref_scan) == running_scan_.end())
        {
            if(landmarks_[i]->getObsNum() <= erase_landmark_nobs_threshold_)
            {
                landmarks_[i]->setState(false);
                std::map<std::shared_ptr<LaserScan>, Eigen::Vector2f> observations = landmarks_[i]->getObservations();
                for(const auto& obs:observations)
                {
                    obs.first->deleteLandmark(landmarks_[i]);
                }
            }
        }
    }


    // merge obs into old landmarks
    for (size_t i = 0; i < landmarks_.size(); ++i) {
        size_t j = landmarks_.size() - i - 1;
        for (size_t k = 0; k < j; ++k) {
            float dis = (landmarks_[k]->getMapPos().head<2>() - landmarks_[j]->getMapPos().head<2>()).norm();
            if(dis < merge_dis_threshold_)
            {
                // merge landmark j into landmark k
                std::map<std::shared_ptr<LaserScan>, Eigen::Vector2f> merge_obs = landmarks_[j]->getObservations();
                for(auto const& ob:merge_obs)
                {
                    // landmark k add observation
                    landmarks_[k]->addObservations(ob.first, ob.second);
                    // scan change its association
                    ob.first->replaceLandmark(landmarks_[j], landmarks_[k]);
                }
                landmarks_[j]->setState(false);
                std::map<std::shared_ptr<LaserScan>, Eigen::Vector2f> observations = landmarks_[j]->getObservations();
                for(const auto& obs:observations)
                {
                    obs.first->deleteLandmark(landmarks_[j]);
                }
            }
        }
    }

    // resize landmarks
    int j = 0;
    for (size_t i = 0; i < landmarks_.size(); ++i) {
        if(landmarks_[i]->getState())
        {
            landmarks_[j++] =  landmarks_[i];
        }
    }
    landmarks_.resize(j);

}


void MapBuilder::discardScan()
{
    assert(probability_grid_map_ != nullptr);
    auto t1 = std::chrono::steady_clock::now();
    probability_grid_map_->discardDynamicObject(scans_);

    int discard_num = scans_.size() -  map_size_;
    for (int i = 0; i < discard_num; ++i) {
        int discard_index = probability_grid_map_->selectUninformativeScan(scans_);
        probability_grid_map_->updateProbMap(scans_[discard_index]);
        scans_.erase(scans_.begin() + discard_index);
        std::cout << "discard scan num " << discard_index << std::endl;
    }
    auto t2 = std::chrono::steady_clock::now();
    auto delta_t = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
    std::cout << "discard scan cost time = " << delta_t.count() * 1000.0 << "ms." << std::endl;
}

void MapBuilder::getGraph(std::vector<Eigen::Vector2f>& nodes,
                      std::vector<std::pair<Eigen::Vector2f, Eigen::Vector2f> >& edges)
{

    std::vector<landmarkEdge> landmarkEdges;

    for(auto & scan : scans_) {
        Eigen::Vector3f pose = scan->getPose();
        nodes.emplace_back(pose.head<2>());
        std::map<std::shared_ptr<LandMark>, Eigen::Vector2f> landmarks = scan->getLandmarks();
        for(auto const& ld:landmarks)
        {
            landmarkEdge le;
            le.laser_ = scan;
            le.landmark_ = ld.first;
            le.mes_ = ld.second;
            landmarkEdges.push_back(le);
        }
    }

    for(auto & poseEdge : poseEdges_) {
        Eigen::Vector3f pose1 = poseEdge.from_scan_->getPose();
        Eigen::Vector3f pose2 = poseEdge.to_scan_->getPose();
        edges.emplace_back(pose1.head<2>(), pose2.head<2>());
    }

    for(auto & landmark : landmarks_) {
        Eigen::Vector2f position = landmark->getMapPos();
        nodes.push_back(position);
    }

    for(auto & landmarkEdge : landmarkEdges) {
        Eigen::Vector3f pose = landmarkEdge.laser_->getPose();
        Eigen::Vector2f position = landmarkEdge.landmark_->getMapPos();
        edges.emplace_back(pose.head<2>(), position);
    }
}

void MapBuilder::addOdom(const Eigen::Vector3f& pose)
{
    if(!got_first_odom_) {
        last_odom_pose_ = pose;
        got_first_odom_ = true;
    }

    odom_pose_ = pose;
//    std::cout << "odom pose " << std::endl << odom_pose_ << std::endl;
}

void MapBuilder::addOdom(const float v, const float omega, const float dt)
{
    if(!got_first_odom_) {
        last_odom_pose_ = Eigen::Vector3f(0,0,0);
        odom_pose_ = Eigen::Vector3f(0,0,0);
        got_first_odom_ = true;
        return;
    }
    float angular_half_delta = odom_pose_[2] + omega * dt / 2.f;
    Eigen::Vector3f increment;
    if(omega > 0.01f)
    {
        increment[0] = (v / omega) * (std::sin(odom_pose_[2] + omega * dt) - std::sin(odom_pose_[2]));
        increment[1] = (v / omega) * std::cos(odom_pose_[2])  - (v / omega) * std::cos(odom_pose_[2] + omega * dt);
        increment[2] =  omega * dt;
    } else
    {
        increment[0] = v * dt * std::cos(angular_half_delta);
        increment[1] = v * dt * std::sin(angular_half_delta);
        increment[2] =  omega * dt;
    }


    Eigen::Matrix3d J_t_minus_1_  = Eigen::Matrix3d::Identity();
    J_t_minus_1_(0, 2) = -v * dt * std::sin(angular_half_delta);
    J_t_minus_1_(1, 2) = v * dt * std::cos(angular_half_delta);

    Eigen::Matrix<double, 3, 2> J_u = Eigen::MatrixXd::Zero(3, 2);
    J_u << dt * std::cos(angular_half_delta), -v * dt * dt * std::sin(angular_half_delta) / 2,
            dt * std::sin(angular_half_delta), v * dt * dt * std::cos(angular_half_delta) / 2,
            0, dt;

    sigma_u_ <<
             (0.1 * v + 0.1) * (0.1 * v + 0.1), 0,
             0, (0.05 * omega + 0.05) * (0.05 * omega + 0.05);

//    std::cout << "sigma u " << std::endl << sigma_u_  << std::endl;
//    std::cout << "J_u " << std::endl << J_u << std::endl;


//    std::cout << "sigma add speed" << std::endl << J_u * sigma_u_ * J_u.transpose() << std::endl;
    // update covariance of t
    sigma_t_ = J_t_minus_1_ * sigma_t_minus_1_ * J_t_minus_1_.transpose() +
            J_u * sigma_u_ * J_u.transpose() + sigma_model_;

//    Eigen::AngleAxisf rotation(odom_pose_[2], Eigen::Vector3f(0, 0, 1));

    odom_pose_ = odom_pose_ + increment;

//    std::cout << "increment pose " << std::endl << increment << std::endl;
//    std::cout << "last odom pose " << std::endl << last_odom_pose_ << std::endl;
//    std::cout << "odom pose " << std::endl << odom_pose_ << std::endl;
}

void MapBuilder::addLaserScan(std::shared_ptr<LaserScan> laser_scan, std::vector<Eigen::Vector2f>& reflective_markers)
{
    if(!got_first_scan_) {
        laser_scan->setId(scans_.size());
        laser_scan->setPose(Eigen::Vector3f(0,0,0));
        laser_scan->transformPointCloud();
        scans_.push_back(laser_scan);
        running_scan_.push_back(laser_scan);

        // add landmarks in first scan into the map
        for (size_t i = 0; i < reflective_markers.size(); ++i) {
            Eigen::Vector2f landmark_in_map = laser_scan->transformLandmark(reflective_markers[i]);
            std::shared_ptr<LandMark> landmark(new LandMark(landmark_in_map, reflective_markers[i], laser_scan));
            laser_scan->addLandmark(reflective_markers[i], landmark);
            landmarks_.push_back(landmark);
        }

        last_scan_ = laser_scan;
        got_first_scan_ = true;
//        addVertex(laser_scan);
        return;
    }

//    Eigen::Vector3f last_scan_pose = scans_.back()->getPose();
    Eigen::Vector3f last_scan_pose = last_scan_->getPose();
    Eigen::AngleAxisf rotation(last_scan_pose[2] - last_odom_pose_[2], Eigen::Vector3f(0, 0, 1));
    Eigen::Vector3f odom_increment = odom_pose_ - last_odom_pose_;
    Eigen::Vector3f update_pose = last_scan_pose + rotation * odom_increment;
//    std::cout << "last_scan_pose pose " << std::endl << last_scan_pose << std::endl;
//    std::cout << "update pose " << std::endl << update_pose << std::endl;
    laser_scan->setPose(update_pose);

    std::cout << "landmarks size " << landmarks_.size() << std::endl;
    last_odom_pose_ = odom_pose_;


    // 1. find association for landmarks

    // construct a running landmarks
    std::vector<std::shared_ptr<LandMark>> running_landmarks;
    for (size_t i = 0; i < running_scan_.size(); ++i) {
        std::map<std::shared_ptr<LandMark>, Eigen::Vector2f> landmarks_associations = running_scan_[i]->getLandmarks();
        for(auto const& asso:landmarks_associations)
        {
            if(std::find(running_landmarks.begin(), running_landmarks.end(), asso.first) == running_landmarks.end())
            {
                running_landmarks.push_back(asso.first);
            }
        }
    }


    std::map<int, std::shared_ptr<LandMark>> associations_all, associations_ransac;

    // find reflective markers association in landmarks
    associateLandmarksRANSAC(laser_scan->getPose(), reflective_markers, running_landmarks, associations_all, associations_ransac);




    auto t1 = std::chrono::steady_clock::now();
    if(use_correlative_scan_matcher_) {
        Eigen::Vector3f csm_pose;
        Eigen::Matrix3f covariance;

        correlative_scan_matcher_.lowResMatchScan(laser_scan, running_scan_, csm_pose, covariance);
        laser_scan->setPose(csm_pose);

        Eigen::Vector3f gnsm_pose  = gauss_newten_scan_matcher_.matchScanToMap(csm_pose, correlative_scan_matcher_.getCorrelativeGrid(),
                                                                               laser_scan->getRawPointCloud());

        float delta_trans = (gnsm_pose.head<2>() - csm_pose.head<2>()).norm();
        float delta_angle = normalizeAngle(gnsm_pose[2] - csm_pose[2]);
        if(fabs(delta_trans) < 0.04 && fabs(delta_angle) < 0.03) {
            laser_scan->setPose(gnsm_pose);
        }
        else {
            laser_scan->setPose(csm_pose);
            std::cout << "Gauss newten scan matcher result got large jump with delta_trans = " << delta_trans
                      << " delta_angle = " << delta_angle << std::endl;
        }
    }
    else {


        Eigen::Matrix3d sigma_update = sigma_t_;
        Eigen::Vector3f gnsm_pose = gauss_newten_scan_matcher_.matchScan(laser_scan, running_scan_, reflective_markers, associations_ransac,
                sigma_t_, sigma_reflector_, sigma_scan_,sigma_update);

        float delta_trans = (gnsm_pose.head<2>() - update_pose.head<2>()).norm();
        float delta_angle = normalizeAngle(gnsm_pose[2] - update_pose[2]);
        if(fabs(delta_trans) < 0.15 && fabs(delta_angle) < 0.15) {
            laser_scan->setPose(gnsm_pose);

        }
        else {
            std::cout << "Gauss newten scan matcher result got large jump with delta_trans = " << delta_trans
                      << " delta_angle = " << delta_angle << std::endl;
        }
        std::stringstream ss;
        ss << "gnsm_pose pose " << std::endl << gnsm_pose << std::endl;
        std::cout << ss.str();
        sigma_t_ = sigma_update;
        sigma_t_minus_1_ = sigma_t_;
    }

    auto t2 = std::chrono::steady_clock::now();
    auto delta_t = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
//    std::cout << "gauss newton cost time = " << delta_t.count() * 1000.0 << "ms." << std::endl;

    last_scan_ = laser_scan;
    laser_scan->setId(-1);
//    laser_scan->transformPointCloud();

    if(!checkPose(laser_scan)) {
        return;
    }


    if(associations_ransac.size() < 2 && keyframe_select_cnt_ < max_keyframe_select_cnt_)
    {
        keyframe_select_cnt_++;
        return;
    }

    keyframe_select_cnt_ = 0;
    laser_scan->setId(scans_.size());
    laser_scan->transformPointCloud();
    scans_.push_back(laser_scan);
    addRunningScan(laser_scan);


    // 4. add new landmarks to map(except for outliers), update existing landmarks
    for (size_t k = 0; k < reflective_markers.size(); ++k) {
        // if reflective can find association in map
        if(associations_ransac.find(k) != associations_ransac.end())
        {
            associations_ransac[k]->addObservations(laser_scan, reflective_markers[k]);
            laser_scan->addLandmark(reflective_markers[k],associations_ransac[k]);
        }
        // if reflective cannot find association in map
        else if(associations_all.find(k) == associations_all.end())
//        else
        {
            Eigen::Vector2f landmark_in_map = laser_scan->transformLandmark(reflective_markers[k]);
//            for (const auto& asso:associations_ransac)
            for (size_t i = 0; i < running_landmarks.size(); ++i)
            {
                float dis = (landmark_in_map.head<2>() - running_landmarks[i]->getMapPos().head<2>()).norm();
                if(dis < merge_dis_threshold_ && running_landmarks[i]->getObsNum() >= erase_landmark_nobs_threshold_)
                {
                    std::cout << "do not add new landmark" << std::endl;
                    continue;
                }
            }

            std::shared_ptr<LandMark> landmark(new LandMark(landmark_in_map, reflective_markers[k], laser_scan));
            laser_scan->addLandmark(reflective_markers[k], landmark);
            landmarks_.push_back(landmark);
            // used when merging landmarks
            running_landmarks.push_back(landmark);
        }
        // else: reflective found association, but get rejected, we don't want these points
    }

    // 5. discard bad landmarks
//    discardLandmarks();

    // push measurement into pose measurements
    poseEdge pe;
    pe.from_scan_ = scans_[laser_scan->getId() - 1];
    pe.to_scan_ = laser_scan;
    Eigen::AngleAxisf rotation_constraint(-scans_[laser_scan->getId() - 1]->getPose()[2], Eigen::Vector3f(0, 0, 1));
    Eigen::Vector3f delta_pose = rotation_constraint * (laser_scan->getPose() - scans_[laser_scan->getId() - 1]->getPose());
    pe.mes_ = delta_pose;
//    float dis = (delta_pose.head<2>()).norm();
//    float dis_info = 1.f / std::pow(0.05f * dis,2.f);
//    if(dis < 0.5f)
//    {
//        dis_info = 1.f / std::pow(0.05f * 0.05,2.f);
//    }
//    float angle_info = 1.f / std::pow(0.03f * std::abs(delta_pose(2)),2.f);
//    if(std::abs(delta_pose(2)) < 0.05f)
//    {
//        angle_info = 1.f / std::pow(0.03f * 0.05f,2.f);
//    }
    pe.info_  = (sigma_t_ + sigma_k_minus_1_).inverse();
    std::cout << "odom info " << std::endl << pe.info_ << std::endl;
    poseEdges_.push_back(pe);

    sigma_k_minus_1_ = sigma_t_;


    detectLoopClosure(laser_scan, reflective_markers);
    bool do_global_adjustment = false;
    if(new_loop_constraint_count_ >= optimize_every_n_constraint_) {
        do_global_adjustment = doPoseAdjustment();
//        alignLandmarkToMap();
//        mergeLandmarks();
        optimize_time_ = std::chrono::steady_clock::now();
        new_loop_constraint_count_ = 0;
    }

    if(new_loop_constraint_count_ > 0) {
        auto delta_t = std::chrono::duration_cast<std::chrono::duration<double>>(
                       std::chrono::steady_clock::now() - optimize_time_);
        std::cout << "global optimization time = " << delta_t.count() * 1000.0 << "ms." << std::endl;

        if(delta_t.count() > 10.0) {
            do_global_adjustment = doPoseAdjustment();
//            alignLandmarkToMap();
//            mergeLandmarks();
            optimize_time_ = std::chrono::steady_clock::now();
            new_loop_constraint_count_ = 0;
        }
    }

    mergeLandmarks(running_landmarks);

    if(!do_global_adjustment && running_scan_.size() > local_optimize_min_scans_)
    {
        doLocalPoseAdjustment(running_landmarks);
    }

}


bool MapBuilder::localizeLaserScan(const Eigen::Vector3f& initial_pose, std::shared_ptr<LaserScan> laser_scan, std::vector<Eigen::Vector2f>& reflective_markers)
{
    if(!got_first_scan_) {
        laser_scan->setId(-1);
        laser_scan->setPose(initial_pose);
        laser_scan->transformPointCloud();
        running_scan_.push_back(laser_scan);
        got_first_scan_ = true;
        last_scan_ = laser_scan;

//        addVertex(laser_scan);
        return true;
    }


    auto t1 = std::chrono::steady_clock::now();
    Eigen::Vector3f last_scan_pose = last_scan_->getPose();
    Eigen::AngleAxisf rotation(last_scan_pose[2] - last_odom_pose_[2], Eigen::Vector3f(0, 0, 1));
    Eigen::Vector3f odom_increment = odom_pose_ - last_odom_pose_;
    Eigen::Vector3f update_pose = last_scan_pose + rotation * odom_increment;
//    std::cout << "last_scan_pose pose " << std::endl << last_scan_pose << std::endl;
//    std::cout << "update pose " << std::endl << update_pose << std::endl;
    laser_scan->setPose(update_pose);
    last_odom_pose_ = odom_pose_;


    std::map<int, std::shared_ptr<LandMark>> associations_all, associations_ransac;

    // find reflective markers association in landmarks
    associateLandmarksRANSAC(laser_scan->getPose(), reflective_markers, landmarks_, associations_all, associations_ransac);
    cur_landmarks_.clear();
    for(const auto& asso:associations_ransac)
    {
        cur_landmarks_.push_back(reflective_markers[asso.first]);
    }

    // get closest scans
    std::vector<std::shared_ptr<LaserScan>> map_scans;
    if(!getClosestScans(laser_scan, map_scans))
    {
        std::cout << "cannot get closest scans" << std::endl;
        return false;
    }

    // localize
    Eigen::Matrix3d sigma_update = sigma_t_;
    double cost;
    Eigen::Vector3f gnsm_pose = gauss_newten_scan_matcher_.matchScan(laser_scan, running_scan_, map_scans,
            reflective_markers, associations_ransac, sigma_t_, sigma_reflector_, sigma_scan_, sigma_update, cost);

//    if(prev_cost_ < 10)
//    {
//        prev_cost_ = cost;
//    }
    std::stringstream ss;
    cost = cost / laser_scan->getRawPointCloud().size();
    ss << "laser scan size: " << laser_scan->getRawPointCloud().size();
    ss << "   cost: " << cost << std::endl;

    if(cost < err_cost_)
    {
        ss << "  map should be rebuilt  ";
        std::cerr << "map should be rebuilt "  << gnsm_pose.x() << " " <<
                  gnsm_pose.y() << " " << gnsm_pose.z() << std::endl;
    }else if(cost < alert_cost_)
    {
        ss << "  map changes a lot  ";
        std::cerr << "map changes a lot around map point "  << gnsm_pose.x() << " " <<
                  gnsm_pose.y() << " " << gnsm_pose.z() << std::endl;
    }


    float delta_trans = (gnsm_pose.head<2>() - update_pose.head<2>()).norm();
    float delta_angle = normalizeAngle(gnsm_pose[2] - update_pose[2]);
//    if(fabs(delta_trans) < 0.2 && fabs(delta_angle) < 0.15 && std::abs(cost - prev_cost_) / prev_cost_ < 0.5 )
    if(fabs(delta_trans) < 0.2 && fabs(delta_angle) < 0.15)
    {
        prev_cost_ = cost;
        laser_scan->setPose(gnsm_pose);

        ss << "  gnsm_pose pose "  << gnsm_pose.x() << " " <<
        gnsm_pose.y() << " " << gnsm_pose.z() << std::endl;
        std::cout << ss.str();
        sigma_t_ = sigma_update;
        sigma_t_minus_1_ = sigma_t_;

        last_scan_ = laser_scan;
        laser_scan->setId(-1);
//    laser_scan->transformPointCloud();

        auto t2 = std::chrono::steady_clock::now();
        auto delta_t = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
//        std::cout << "localization cost time = " << delta_t.count() * 1000.0 << "ms." << std::endl;


        if(!checkPose(laser_scan)) {
            return true;
        }

        laser_scan->transformPointCloud();
//    if(probability_grid_map_->calculateCrossEntropy(laser_scan, 1500))
//    {
//        std::cout << "add one scan " << scans_.size() << std::endl;
//        laser_scan->setId(scans_.size());
//        scans_.push_back(laser_scan);
//    }
        addRunningScan(laser_scan);
        return true;
    }
    else {
        ss << "   Gauss newten scan matcher result either got large jump with delta_trans = " << delta_trans
                  << " delta_angle = " << delta_angle << std::endl << "or got large jump with cost " << cost - prev_cost_ <<
                  std::endl;
        std::cout << ss.str();
        return false;
    }


}


std::vector<Eigen::Vector3f> MapBuilder::getPath()
{
    std::vector<Eigen::Vector3f> path;
    for(const std::shared_ptr<LaserScan>& scan : scans_) {
        path.push_back(scan->getPose());
    }

    return path;
}




std::vector<std::shared_ptr<LandMark>> MapBuilder::getLandmarks()
{
    return landmarks_;
}

std::vector<Eigen::Vector2f> MapBuilder::getCurLandmarks()
{
    return cur_landmarks_;
}



std::shared_ptr<CorrelativeGrid> MapBuilder::getCorrelativeGrid()
{
    if(use_correlative_scan_matcher_) {
        return correlative_scan_matcher_.getCorrelativeGrid();
    }
    else {
        return gauss_newten_scan_matcher_.getCorrelativeGrid();
    }
}

std::shared_ptr<OccupancyGridMap> MapBuilder::getOccupancyGridMap()
{
    Range map_range;
    for(const std::shared_ptr<LaserScan>& scan : scans_) {
        map_range.addRange(scan->getRange());
    }

    const Eigen::Vector2f& max = map_range.getMax();
    const Eigen::Vector2f& min = map_range.getMin();
    int width = ceil((max[0] - min[0]) / occupancy_grid_map_resolution_);
    int height = ceil((max[1] - min[1]) / occupancy_grid_map_resolution_);

    std::shared_ptr<OccupancyGridMap> occupancy_grid_map(new OccupancyGridMap(width, height, occupancy_grid_map_resolution_));
    occupancy_grid_map->setOrigin(min);

    std::vector<std::shared_ptr<LaserScan>> all_scan = scans_;
    occupancy_grid_map->createFromScan(all_scan);

    return occupancy_grid_map;
}

std::shared_ptr<ProbabilityGridMap> MapBuilder::getProbabilityGridMap()
{
    Range map_range;
    for(const std::shared_ptr<LaserScan>& scan : scans_) {
        map_range.addRange(scan->getRange());
    }

    const Eigen::Vector2f& max = map_range.getMax();
    const Eigen::Vector2f& min = map_range.getMin();
    int width = std::ceil((max[0] - min[0]) / occupancy_grid_map_resolution_) + 50;
    int height = std::ceil((max[1] - min[1]) / occupancy_grid_map_resolution_) + 50;

    std::shared_ptr<ProbabilityGridMap> probability_grid_map(new ProbabilityGridMap(width, height, occupancy_grid_map_resolution_));
    probability_grid_map->setOrigin(min);

    std::vector<std::shared_ptr<LaserScan>> all_scan = scans_;
    probability_grid_map->createFromScan(all_scan);


    // assign probability_grid_map to probability_grid_map_ for lifelong mapping purpose
    probability_grid_map_ = probability_grid_map;
    return probability_grid_map;
}


int MapBuilder::pointRansac(const float& disTh, const std::vector<Eigen::Vector2f>& cands, std::vector<int>& flags)
{
    int best_inlier_total = 0;
    for (size_t i = 0; i < cands.size(); ++i) {
        Eigen::Vector2f refPose = cands[i];
        int inlier_total = 0;
        std::vector<int> inlier_flags(cands.size());
        for (size_t j = 0; j < cands.size(); ++j) {
            Eigen::Vector2f propPose = cands[j];
            Eigen::Vector2f relPose = refPose - propPose;
            float dis = sqrt(relPose.x() * relPose.x() + relPose.y() * relPose.y()) ;
            if (dis < disTh)
            {
                inlier_flags[j] = 1;
                inlier_total++;
            } else
            {
                inlier_flags[j] = 0;
            }
        }
        if(inlier_total > best_inlier_total)
        {
            best_inlier_total = inlier_total;
            flags = inlier_flags;
        }
    }
    return best_inlier_total;
}



bool MapBuilder::poseRansac(const float& disTh, const std::vector<Eigen::Vector2f>& mes, const std::vector<Eigen::Vector2f>& est,
        std::vector<int>& flags, Eigen::Vector3f& pose_ransac)
{
    const float confidence_level = 0.99;

    int best_inlier_total = 0;
    int asso_total = est.size();
    assert(asso_total > 2);
//        cout << "pose_total " << pose_total << endl;
    int sampling_total = asso_total;
    int sampling_cnt = 0;
    std::map<std::pair<int, int>, int> draw_table;
    while (sampling_cnt < sampling_total){
        int draw_ind[2];
        draw_ind[0] = std::rand() % asso_total;
        draw_ind[1] = std::rand() % asso_total;
        Eigen::Matrix4f A;
        int i = 0;
        while(i < 50)
        {
            if(draw_ind[1] != draw_ind[0] && draw_table.find(std::make_pair(draw_ind[0], draw_ind[1])) == draw_table.end() && draw_table.find(std::make_pair(draw_ind[0], draw_ind[1])) == draw_table.end())
            {

                A <<
                  est[draw_ind[0]].x(), -est[draw_ind[0]].y(), 1, 0,
                        est[draw_ind[0]].y(),  est[draw_ind[0]].x(), 0, 1,
                        est[draw_ind[1]].x(), -est[draw_ind[1]].y(), 1, 0,
                        est[draw_ind[1]].y(),  est[draw_ind[1]].x(), 0, 1;
                if(A.determinant() > 1e-2)
                {
                    break;
                }
            }
            draw_ind[1] = std::rand() % asso_total;
            i++;
        }
        if(i >= 50)
        {
            std::cout << "pose ransac fail" << std::endl;
            return false;
        }
        draw_table[std::make_pair(draw_ind[0], draw_ind[1])] = 1;
        // count the number of inliers
        int inlier_total = 0;
        std::vector<int> inlier_flags(asso_total);

        Eigen::Vector4f b;
        b << mes[draw_ind[0]].x(), mes[draw_ind[0]].y(),mes[draw_ind[1]].x(),mes[draw_ind[1]].y();
//            std::cout << "A = " << A << std::endl;
//            std::cout << "b = " << b << std::endl;

        Eigen::Vector4f X = A.inverse() * b;
//            std::cout << "x = " << X(2) << " y = " << X(3) << " theta = " << atan2(X(1), X(0)) << std::endl;

        Eigen::Matrix2f rot;
        float det = std::sqrt(X(0)*X(0) + X(1)*X(1));
        rot << X(0), -X(1),
                X(1), X(0);
        rot = rot / det;
//        Eigen::Matrix3d trans;
//        trans <<
//                rot(0,0), rot(0,1), X(2),
//                rot(1,0), rot(1,1), X(3),
//                0, 0, 1;


        for (int idx = 0; idx < asso_total; idx++){
            Eigen::Vector2f pred = rot * est[idx] + Eigen::Vector2f(X(2), X(3));
            float dist_diff = (pred.x() - mes[idx].x())*(pred.x() - mes[idx].x())
                               + (pred.y() - mes[idx].y())*(pred.y() - mes[idx].y());
            if (dist_diff > disTh*disTh){
                inlier_flags[idx] = 0;
                continue;
            }
            else{
                ++inlier_total;
                inlier_flags[idx] = 1;
            }

        }
//        std::cout << "inliers_total " << inlier_total << std::endl;
        if (inlier_total > best_inlier_total){
            best_inlier_total = inlier_total;
            flags = inlier_flags;
            pose_ransac << X(2), X(3), atan2(X(1), X(0));

            sampling_total = ComputeRequiredSamplingTotal(2, best_inlier_total,
                                                          asso_total, sampling_total, confidence_level);
        }
//            cout << "sampling_cnt " << sampling_cnt << " sampling_total " << sampling_total << endl;
        ++sampling_cnt;
    }
    return true;
}


int MapBuilder::ComputeRequiredSamplingTotal(const int draw_total, const int inlier_total,
                                             const int pose_total, const int current_sampling_total, const float confidence_level){
    float ep = 1 - static_cast<float>(inlier_total)/static_cast<float>(pose_total);
    if (ep == 1.0) {
        ep = 0.5;
    }
    int newSamplingTotal = (int)(std::log(1 - confidence_level)/std::log(1 - pow(1 - ep, draw_total)) + 0.5);
    if (newSamplingTotal < current_sampling_total) {
        return newSamplingTotal;
    } else {
        return current_sampling_total;
    }
}





void MapBuilder::associateLandmarksRANSAC(const Eigen::Vector3f& prior_pose,
        const std::vector<Eigen::Vector2f>& reflective_markers,
      const std::vector<std::shared_ptr<LandMark>>& reference_landmarks,
      std::map<int, std::shared_ptr<LandMark>>& associations_all,
      std::map<int, std::shared_ptr<LandMark>>& associations_ransac)
{
    const float detect_threshold = max_landmark_detect_err_ * max_landmark_detect_err_;
//    std::vector<std::pair<std::shared_ptr<LandMark>, int>> associations;
    std::map<int, std::shared_ptr<LandMark>> associations;
    std::vector<Eigen::Vector2f> reflective_markers_in_map;
    Eigen::Affine2f transform(Eigen::Translation2f(prior_pose[0], prior_pose[1]) * Eigen::Rotation2Df(prior_pose[2]));

    for (size_t i = 0; i < reflective_markers.size(); ++i) {
        Eigen::Vector2f landmark_in_map = transform * reflective_markers[i];
        reflective_markers_in_map.push_back(landmark_in_map);
        float min_dis = detect_threshold;
        std::shared_ptr<LandMark> best_associate;
        for (size_t j = 0; j < reference_landmarks.size(); ++j) {
            Eigen::Vector2f landmark_candidate = reference_landmarks[j]->getMapPos();
            float dis = (landmark_in_map.x() - landmark_candidate.x()) * (landmark_in_map.x() - landmark_candidate.x())
                        + (landmark_in_map.y() - landmark_candidate.y()) * (landmark_in_map.y() - landmark_candidate.y());
            if(dis < min_dis)
            {
                min_dis = dis;
                best_associate = reference_landmarks[j];
            }
        }
        if(min_dis < (detect_threshold - 0.0001f))
        {
            associations[i] = best_associate;
        }
    }


    // 2. reverse association(only consistent one stays)
    std::vector<Eigen::Vector2f> mes, est;
    for(auto const& asso:associations)
    {
        Eigen::Vector2f landmark_candidate = asso.second->getMapPos();
        float min_dis = detect_threshold;
        int best_id = -1;
        for (size_t i = 0; i < reflective_markers_in_map.size(); ++i) {
            float dis = (reflective_markers_in_map[i].x() - landmark_candidate.x()) * (reflective_markers_in_map[i].x() - landmark_candidate.x())
                        + (reflective_markers_in_map[i].y() - landmark_candidate.y()) * (reflective_markers_in_map[i].y() - landmark_candidate.y());
            if(dis < min_dis)
            {
                min_dis = dis;
                best_id = i;
            }
        }
        if(best_id == asso.first)
        {
            associations_all[best_id] = asso.second;
            mes.push_back(landmark_candidate);
            est.push_back(reflective_markers[best_id]);
        }
    }



    // 3. ransac (reject outliers)
    std::vector<int> flags;
    bool pose_ransac_refine = false;
    if(mes.size() > 2)
    {
        Eigen::Vector3f pose_ransac;
        pose_ransac_refine = poseRansac(ransac_dis_threshold_, mes, est, flags, pose_ransac);
        Eigen::Vector3f jump = pose_ransac - prior_pose;
        if(pose_ransac_refine && std::accumulate(flags.begin(), flags.end(), 0) > 2 &&
        std::sqrt(jump.x() * jump.x() + jump.y() * jump.y()) < ransac_jump_dis_threshold_ &&
        std::abs(jump.z()) < ransac_jump_angle_threshold_)
        {
            assert(flags.size() == associations_all.size());
            utility::reduceVector(mes, flags);
            utility::reduceVector(est, flags);
            int id = 0;
            for(auto const& asso:associations_all)
            {
                if(flags[id])
                {
                    associations_ransac[asso.first] = asso.second;
                }
                id++;
            }
        }

    }
}



bool MapBuilder::getLandmarkByID(const size_t id, std::shared_ptr<LandMark>& landmark)
{
    for(size_t i = 0; i < landmarks_.size(); ++i) {
        if(id == landmarks_[i]->getId())
        {
            landmark = landmarks_[i];
            return true;
        }
    }
    return false;
}


void MapBuilder::saveScans(const std::string& save_path, const Eigen::Vector3f& lidar_to_odom)
{
    std::vector<landmarkEdge> landmarkEdges;
    std::string POSE_GRAPH_SAVE_PATH = save_path + "keyframes/";

    // Check if the directory exists, if so, clear it; otherwise, create it
    std::filesystem::path path_to_save(POSE_GRAPH_SAVE_PATH);
    if (std::filesystem::exists(path_to_save)) {
        for (const auto& entry : std::filesystem::directory_iterator(path_to_save)) {
            std::filesystem::remove_all(entry.path());
        }
    } else {
        std::filesystem::create_directories(path_to_save);
    }

    // Save scans to "pose_graph.txt"
    std::ofstream pose_file(POSE_GRAPH_SAVE_PATH + "pose_graph.txt");
    if (!pose_file.is_open()) {
        std::cerr << "Failed to open pose_graph.txt for writing." << std::endl;
        return;
    }
    pose_file << lidar_to_odom(0) << " " << lidar_to_odom(1) << " " << lidar_to_odom(2) << "\n";

    for (const auto& scan : scans_) {
        Eigen::Vector3f pose = scan->getPose();
        PointCloud scan_data = scan->getRawPointCloud();

        pose_file << scan->getId() << " " << pose.x() << " " << pose.y() << " " << pose.z() << " " << scan_data.size() << "\n";

        // Save individual scan data
        std::ofstream scan_file(POSE_GRAPH_SAVE_PATH + std::to_string(scan->getId()) + "_scan.txt");
        if (!scan_file.is_open()) {
            std::cerr << "Failed to open scan file for ID " << scan->getId() << std::endl;
            continue;
        }
        for (const auto& point : scan_data) {
            scan_file << point.x() << " " << point.y() << "\n";
        }

        // Collect landmark edges
        auto landmarks = scan->getLandmarks();
        for (const auto& ld : landmarks) {
            landmarkEdge le;
            le.laser_ = scan;
            le.landmark_ = ld.first;
            le.mes_ = ld.second;
            le.info_ = sigma_reflector_.inverse();
            landmarkEdges.push_back(le);
        }
    }

    // Save landmarks to "landmarks.txt"
    std::ofstream landmarks_file(POSE_GRAPH_SAVE_PATH + "landmarks.txt");
    if (!landmarks_file.is_open()) {
        std::cerr << "Failed to open landmarks.txt for writing." << std::endl;
        return;
    }
    for (const auto& landmark : landmarks_) {
        Eigen::Vector2f pose = landmark->getMapPos();
        int id = landmark->getId();
        landmarks_file << id << " " << pose.x() << " " << pose.y() << "\n";
    }

    // Save pose constraints to "pose_constraints.txt"
    std::ofstream pose_constraints_file(POSE_GRAPH_SAVE_PATH + "pose_constraints.txt");
    if (!pose_constraints_file.is_open()) {
        std::cerr << "Failed to open pose_constraints.txt for writing." << std::endl;
        return;
    }
    for (const auto& edge : poseEdges_) {
        Eigen::Vector3f pose = edge.mes_;
        Eigen::Matrix3d info = edge.info_;
        int front_scan = edge.from_scan_->getId();
        int to_scan = edge.to_scan_->getId();
        pose_constraints_file << front_scan << " " << to_scan << " "
                              << pose.x() << " " << pose.y() << " " << pose.z() << " "
                              << info(0, 0) << " " << info(1, 1) << " " << info(2, 2) << "\n";
    }

    // Save landmark constraints to "landmark_constraints.txt"
    std::ofstream landmark_constraints_file(POSE_GRAPH_SAVE_PATH + "landmark_constraints.txt");
    if (!landmark_constraints_file.is_open()) {
        std::cerr << "Failed to open landmark_constraints.txt for writing." << std::endl;
        return;
    }
    for (const auto& edge : landmarkEdges) {
        Eigen::Vector2f pose = edge.mes_;
        Eigen::Matrix2d info = edge.info_;
        int front_scan = edge.laser_->getId();
        int to_landmark = edge.landmark_->getId();
        landmark_constraints_file << front_scan << " " << to_landmark << " "
                                  << pose.x() << " " << pose.y() << " "
                                  << info(0, 0) << " " << info(1, 1) << "\n";
    }

    std::cout << "Scans saved successfully to " << POSE_GRAPH_SAVE_PATH << std::endl;
}





bool MapBuilder::loadScans(const std::string& load_path)
{
    std::cout << "load database ... please wait..." << std::endl;
    std::string POSE_GRAPH_SAVE_PATH = load_path + "keyframes/";

    // 加载 pose graph
    std::string file_path = POSE_GRAPH_SAVE_PATH + "pose_graph.txt";
    std::ifstream pose_file(file_path);
    if (!pose_file.is_open()) {
        std::cerr << "Failed to load pose graph from: " << file_path << std::endl;
        return false;
    }

    Eigen::Vector3f lidar_to_odom;
    pose_file >> lidar_to_odom[0] >> lidar_to_odom[1] >> lidar_to_odom[2];

    int index, size;
    float x, y, theta;
    while (pose_file >> index >> x >> y >> theta >> size) {
        std::string scan_path = POSE_GRAPH_SAVE_PATH + std::to_string(index) + "_scan.txt";
        std::ifstream scan_file(scan_path);
        if (!scan_file.is_open()) {
            std::cerr << "Failed to load scan from: " << scan_path << std::endl;
            return false;
        }

        PointCloud scan_points;
        for (int i = 0; i < size; ++i) {
            float p_x, p_y;
            if (!(scan_file >> p_x >> p_y)) {
                std::cerr << "Failed to read scan point at index " << i << " in file: " << scan_path << std::endl;
                return false;
            }
            scan_points.emplace_back(p_x, p_y);
        }

        auto laser_scan = std::make_shared<RLML::LaserScan>(scan_points);
        laser_scan->setId(index);
        laser_scan->setCalibration(lidar_to_odom);
        laser_scan->setPose(Eigen::Vector3f(x, y, theta));
        laser_scan->transformPointCloud();
        scans_.push_back(laser_scan);
    }

    // 加载 pose constraints
    file_path = POSE_GRAPH_SAVE_PATH + "pose_constraints.txt";
    std::ifstream pose_constraints_file(file_path);
    if (!pose_constraints_file.is_open()) {
        std::cerr << "Failed to load pose constraints from: " << file_path << std::endl;
        return false;
    }

    int index_source, index_to;
    float info1, info2, info3;
    while (pose_constraints_file >> index_source >> index_to >> x >> y >> theta >> info1 >> info2 >> info3) {
        poseEdge pe;
        pe.from_scan_ = scans_[index_source];
        pe.to_scan_ = scans_[index_to];
        pe.mes_ = Eigen::Vector3f(x, y, theta);
        pe.info_ << info1, 0, 0, 0, info2, 0, 0, 0, info3;
        poseEdges_.push_back(pe);
    }

    // 加载 landmarks
    file_path = POSE_GRAPH_SAVE_PATH + "landmarks.txt";
    std::ifstream landmarks_file(file_path);
    if (!landmarks_file.is_open()) {
        std::cerr << "Failed to load landmarks from: " << file_path << std::endl;
        return false;
    }

    while (landmarks_file >> index >> x >> y) {
        auto landmark = std::make_shared<LandMark>(index, Eigen::Vector2f(x, y));
        landmarks_.push_back(landmark);
    }

    // 加载 landmark constraints
    file_path = POSE_GRAPH_SAVE_PATH + "landmark_constraints.txt";
    std::ifstream landmark_constraints_file(file_path);
    if (!landmark_constraints_file.is_open()) {
        std::cerr << "Failed to load landmark constraints from: " << file_path << std::endl;
        return false;
    }

    size_t index_laser, index_landmark;
    while (landmark_constraints_file >> index_laser >> index_landmark >> x >> y >> info1 >> info2) {
        auto it = std::find_if(landmarks_.begin(), landmarks_.end(), [index_landmark](const std::shared_ptr<LandMark>& lm) {
            return lm->getId() == index_landmark;
        });
        if (it != landmarks_.end()) {
            scans_[index_laser]->addLandmark(Eigen::Vector2f(x, y), *it);
        }
    }

    // 计算 landmarks 距离表
    landmark_dis_tab_.reserve(landmarks_.size());
    for (size_t i = 0; i < landmarks_.size(); ++i) {
        std::map<int, float> landmark_dis;
        for (size_t j = 0; j < landmarks_.size(); ++j) {
            if (i == j) continue;
            float dis = (landmarks_[i]->getMapPos() - landmarks_[j]->getMapPos()).norm();
            if (dis < 60.f) {
                landmark_dis[j] = dis;
            }
        }
        landmark_dis_tab_.push_back(landmark_dis);
    }

    map_size_ = scans_.size();
    std::cout << "load database done" << std::endl;
    std::cout << "keyframes number: " << map_size_ << " reflector size: " << landmarks_.size() << std::endl;
    return true;
}



} // namespace RLML
