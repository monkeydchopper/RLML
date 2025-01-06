#ifndef MAP_BUILDER_H
#define MAP_BUILDER_H


#include <chrono>
#include <numeric>
#include <vector>
#include <Eigen/QR>
#include <random>
#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iostream>

#include "math_func.h"
#include "utility.h"
#include "laser_scan.h"
#include "occupancy_grid_map.h"
#include "probability_grid_map.h"
#include "correlative_scan_matcher.h"
#include "gauss_newten_scan_matcher.h"
#include "g2o/core/sparse_optimizer.h"
#include "g2o/core/block_solver.h"
#include "g2o/core/robust_kernel_impl.h"
#include "g2o/solvers/cholmod/linear_solver_cholmod.h"
#include "g2o/core/optimization_algorithm_levenberg.h"
#include "g2o/types/slam2d/types_slam2d.h"

namespace RLML
{

class landmark
{
public:
    int id_;
    Eigen::Vector2f pos_;
    std::vector<std::pair<int, Eigen::Vector2f>> associate_in_scan_;
};

struct landmarkEdge
{
    std::shared_ptr<LaserScan> laser_;
    std::shared_ptr<LandMark> landmark_;
    Eigen::Vector2f mes_;
    Eigen::Matrix2d info_;
};
struct poseEdge
{
    std::shared_ptr<LaserScan> from_scan_;
    std::shared_ptr<LaserScan> to_scan_;
    Eigen::Vector3f mes_;
    Eigen::Matrix3d info_;
};


class MapBuilder
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    MapBuilder();
    ~MapBuilder() {}

    void initialize();
    void addOdom(const Eigen::Vector3f& pose);
    void addOdom(const float v, const float omega, const float t);
    void addLaserScan(std::shared_ptr<LaserScan> laser_scan, std::vector<Eigen::Vector2f>& reflective_markers);
    bool localizeLaserScan(const Eigen::Vector3f& initial_pose,
            std::shared_ptr<LaserScan> laser_scan,
            std::vector<Eigen::Vector2f>& reflective_markers);
    bool relocalize(const std::shared_ptr<LaserScan>& scan, const std::vector<Eigen::Vector2f>& reflective_markers, const bool use_reflector_relocalize);
    bool localRelocalize(const std::shared_ptr<LaserScan>& scan);

    std::shared_ptr<OccupancyGridMap> getOccupancyGridMap();
    std::shared_ptr<ProbabilityGridMap> getProbabilityGridMap();
    std::shared_ptr<CorrelativeGrid> getCorrelativeGrid();
    void getGraph(std::vector<Eigen::Vector2f>& nodes,
                  std::vector<std::pair<Eigen::Vector2f, Eigen::Vector2f> > &edges);
    std::vector<Eigen::Vector3f> getPath();
    std::vector<Eigen::Vector4d> getLandmarkAreas();
    std::vector<std::shared_ptr<LandMark>> getLandmarks();
    std::vector<Eigen::Vector2f> getCurLandmarks();
    std::vector<Eigen::Vector2f> getMatchedKeypoints();
    std::vector<Eigen::Vector2f> getCurKeypoints();
    bool doPoseAdjustment();
    bool doLocalPoseAdjustment(const std::vector<std::shared_ptr<LandMark>>& running_landmarks);
    void alignLandmarkToMap();
    void mergeLandmarks(const std::vector<std::shared_ptr<LandMark>>& running_landmarks);
    void updateLandmarks();
    void discardLandmarks();
    void discardScan();

    void setCostThreshold(float err_cost, float alert_cost) {
        err_cost_ = err_cost;
        alert_cost_ = alert_cost;
    }
    void setOccupancyGridMapResolution(float res) { occupancy_grid_map_resolution_ = res; }
    void setMinUpdateDistance(float distance) { min_update_distance_ = distance; }
    void setMinUpdateOrientation(float angle) { min_update_orientation_ = angle; }
    void setScanBufferSize(int num) { scan_buffer_size_ = num; }
    void setLoopScanSearchDistance(float distance) { loop_scan_search_distance_ = distance; }
    void setLoopMatchMinChainSize(int size) { loop_match_min_chain_size_ = size; }
    void setLoopClosureMinResponse(float response) { loop_closure_min_response_ = response; }
    void setLoopClosureXYVarianceThreshold(float variance) { loop_closure_xy_variance_threshold_ = variance; }
    void setLoopClosureAngleVarianceThreshold(float variance) { loop_closure_angle_variance_threshold_ = variance; }
    void setOptimizeEveryNConstraints(int n) { optimize_every_n_constraint_ = n; }

    void setLoopClosureXYSearchRange(float range) { loop_closure_xy_search_range_ = range; }
    void setLoopClosureAngleSearchRange(float range) { loop_closure_angle_search_range_ = range; }
    void setLoopClosureGridResolution(float res) { loop_closure_grid_resolution_ = res; }
    void setLoopClosureCoarseXYSearchResolution(float res) { loop_closure_coarse_xy_search_resolution_ = res; }
    void setLoopClosureFineXYSearchRange(float range) { loop_closure_fine_xy_search_range_ = range; }
    void setLoopClosureCoarseAngleSearchResolution(float res) { loop_closure_coarse_angle_search_resolution_ = res; }
    void setLoopClosureFineAngleSearchRange(float range) { loop_closure_fine_angle_search_range_ = range; }
    void setLoopClosureFineAngleSearchResolution(float res) { loop_closure_fine_angle_search_resolution_ = res; }

    void useCorrelativeScanMatcher(bool flag) { use_correlative_scan_matcher_ = flag; }
    void setCSMXYSearchRange(float range) { csm_xy_search_range_ = range; }
    void setCSMAngleSearchRange(float range) { csm_angle_search_range_ = range; }
    void setCSMGridResolution(float res) { csm_grid_resolution_ = res; }
    void setCSMXYSearchResolution(float res) { csm_xy_search_resolution_ = res; }
    void setCSMAngleSearchResolution(float res) { csm_angle_search_resolution_ = res; }
    void saveScans(const std::string& save_path, const Eigen::Vector3f& lidar_to_odom);
    void saveLandmarkArea(const std::string& save_path);
    bool loadScans(const std::string& load_path);
    void setInitialCovariance(const Eigen::Matrix3d& initial_pose_covariance)
    {
        sigma_t_minus_1_ = initial_pose_covariance;
        sigma_t_ = sigma_t_minus_1_;
        std::cout << "set initial pose covariance" << std::endl << sigma_t_minus_1_ << std::endl;
    }
    void setOdomCovariance(const Eigen::Matrix3d& odom_model_covariance)
    {
        sigma_model_ = odom_model_covariance;
        std::cout << "set odom model covariance " << std::endl << sigma_model_ << std::endl;
    }
    void setReflectorCovariance(const Eigen::Matrix2d& reflector_covariance)
    {
        sigma_reflector_ = reflector_covariance;
        std::cout << "set reflector covariance" << std::endl << sigma_reflector_ << std::endl;
    }

    void setRelocalizeCovariance(const Eigen::Matrix3d& relocalize_covariance)
    {
        sigma_relocalize_ = relocalize_covariance;
        std::cout << "set relocalize covariance " << std::endl << sigma_relocalize_ << std::endl;
    }
    void setVelocityCovarianceParam(const double v_covariance_param)
    {
        v_covariance_param_ = v_covariance_param;
        std::cout << "set velocity covariance param " <<  v_covariance_param_ << std::endl;
    }
    void setOmegaCovarianceParam(const double omega_covariance_param)
    {
        omega_covariance_param_ = omega_covariance_param;
        std::cout << "set omega covariance param " <<  omega_covariance_param_ << std::endl;
    }
    void setScanCovariance(const double scan_probability_covariance_param)
    {
        sigma_scan_ = scan_probability_covariance_param * scan_probability_covariance_param;
        std::cout << "set scan covariance "  << sigma_scan_ << std::endl;
    }
    void setRelocalizeMinResponse(const float relocalize_min_response)
    {
        relocalize_min_response_ = relocalize_min_response;
    }

    int getLoopCnt()
    {
        return loop_constraint_count_;
    }

private:
    bool getClosestScans(const std::shared_ptr<LaserScan>& laser_scan, std::vector<std::shared_ptr<LaserScan>>& map_scans);
    bool checkPose(const std::shared_ptr<LaserScan>& laser_scan);
    void addRunningScan(std::shared_ptr<LaserScan> laser_scan);
    void addVertex(g2o::SparseOptimizer& optimizer, const std::shared_ptr<LaserScan>& scan);
    void addVertex(g2o::SparseOptimizer& optimizer, const std::shared_ptr<LandMark>& landmark, const int offset_id);
    void addEdge(g2o::SparseOptimizer& optimizer, std::shared_ptr<LaserScan> source_scan, const Eigen::Vector3f& source_pose,
                 std::shared_ptr<LaserScan> target_scan, const Eigen::Vector3f& target_pose,
                 const Eigen::Matrix3d& information);
    void addEdge(g2o::SparseOptimizer& optimizer,
                             const poseEdge& pose_edge);
    void addEdge(g2o::SparseOptimizer& optimizer,
                             const landmarkEdge& landmark_edge,
                             const int offset_id);
    std::shared_ptr<LaserScan> getClosestScan(const std::shared_ptr<LaserScan>& base_scan,
                                             const std::vector<std::shared_ptr<LaserScan>>& chain);
    void detectLoopClosure(const std::shared_ptr<LaserScan>& scan, const std::vector<Eigen::Vector2f>& reflective_markers);
    int pointRansac(const float& disTh, const std::vector<Eigen::Vector2f>& cands, std::vector<int>& flags);
    bool poseRansac(const float& disTh, const std::vector<Eigen::Vector2f>& mes, const std::vector<Eigen::Vector2f>& est,
            std::vector<int>& flags, Eigen::Vector3f& pose_ransac);
    int ComputeRequiredSamplingTotal(const int draw_total, const int inlier_total,
                                                 const int pose_total, const int current_sampling_total, const float confidence_level);

    void associateLandmarksRANSAC(const Eigen::Vector3f& prior_pose,
            const std::vector<Eigen::Vector2f>& reflective_markers,
            const std::vector<std::shared_ptr<LandMark>>& reference_landmarks,
          std::map<int, std::shared_ptr<LandMark>>& associations_all,
          std::map<int, std::shared_ptr<LandMark>>& associations_ransac);
//    void associateKeypointsRANSAC(const Eigen::Vector3f& prior_pose, const std::vector<std::shared_ptr<LaserScan>>& map_scans, const std::vector<falkolib::FALKO>& cur_keypoints,
//            const std::vector<falkolib::CGH>& cur_desc, std::vector<Eigen::Vector2f>& est, std::vector<Eigen::Vector2f>& mes, std::vector<float>& weightsRansac);
    bool getLandmarkByID(const size_t id, std::shared_ptr<LandMark>& landmark);

//    bool getLandmarkID(const std::shared_ptr<LandMark>& landmark, int& id);

private:


    std::vector<poseEdge> poseEdges_;
    // parameter
    float occupancy_grid_map_resolution_;
    float min_update_distance_;
    float min_update_orientation_;
    size_t scan_buffer_size_;

    float loop_scan_search_distance_;
    float loop_closure_min_response_;
    float loop_closure_xy_variance_threshold_;
    float loop_closure_angle_variance_threshold_;
    size_t loop_match_min_chain_size_;

    float loop_closure_xy_search_range_;
    float loop_closure_angle_search_range_;
    float loop_closure_grid_resolution_;
    float loop_closure_coarse_xy_search_resolution_;
    float loop_closure_fine_xy_search_range_;
    float loop_closure_coarse_angle_search_resolution_;
    float loop_closure_fine_angle_search_range_;
    float loop_closure_fine_angle_search_resolution_;

    int optimize_every_n_constraint_;
    float gn_scan_matcher_grid_resolution_;

    bool use_correlative_scan_matcher_;
    float csm_xy_search_range_;
    float csm_angle_search_range_;
    float csm_grid_resolution_;
    float csm_xy_search_resolution_;
    float csm_angle_search_resolution_;

    float search_space_standard_deviation_;


    int new_loop_constraint_count_, loop_constraint_count_;
    bool got_first_scan_;
    bool got_first_odom_;
    Eigen::Vector3f last_odom_pose_;
    Eigen::Vector3f odom_pose_;

    std::shared_ptr<LaserScan> last_scan_;
    std::vector<std::shared_ptr<LandMark>> landmarks_;
    std::vector<std::map<int, float>> landmark_dis_tab_;
    std::vector<Eigen::Vector2f> cur_landmarks_;
    std::vector<std::shared_ptr<LaserScan>> scans_;
    std::vector<std::shared_ptr<LaserScan>> running_scan_;
    std::shared_ptr<ProbabilityGridMap> probability_grid_map_;

    CorrelativeScanMatcher correlative_scan_matcher_;
    CorrelativeScanMatcher loop_closure_scan_matcher_;
    GaussNewtenScanMatcher gauss_newten_scan_matcher_;

    std::chrono::steady_clock::time_point optimize_time_;


    const float ransac_jump_dis_threshold_ = 0.2f;
    const float ransac_jump_angle_threshold_ = 0.1f;
    const float max_landmark_detect_err_ = 0.5f;
    const float ransac_dis_threshold_ = 0.08f;
    const float merge_dis_threshold_ = 0.2f;
    const size_t erase_landmark_nobs_threshold_ = 6;
    const size_t max_keyframe_select_cnt_ = 4;
    const size_t local_optimize_min_scans_ = 20;



    size_t map_size_;
    size_t edge_count_;
    size_t keyframe_select_cnt_;

    Eigen::Matrix3d sigma_t_minus_1_, sigma_t_;
    Eigen::Matrix2d sigma_u_;
    Eigen::Matrix2d sigma_reflector_, sigma_keypoint_;
    double sigma_scan_, v_covariance_param_, omega_covariance_param_;
    Eigen::Matrix3d sigma_k_minus_1_;
    Eigen::Matrix3d sigma_relocalize_;
    Eigen::Matrix3d sigma_model_;

    float relocalize_min_response_;


    std::vector<Eigen::Vector2f> cur_keypoints_, matched_keypoints_;

    double prev_cost_;
    float err_cost_, alert_cost_;
};

} // namespace RLML

#endif // MAP_BUILDER_H
