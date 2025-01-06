#include "probability_grid_map.h"

namespace RLML
{

void ProbabilityGridMap::createFromScan(const std::vector<std::shared_ptr<LaserScan>>& scans)
{
    for(const std::shared_ptr<LaserScan>& scan : scans) {
//        Eigen::Vector2f start = getMapCoords(scan->getPose());
        Eigen::Vector2f start = getMapCoords(scan->getScanPose());
        const PointCloud& point_cloud = scan->getTransformedPointCloud();
        for(const Eigen::Vector2f& point : point_cloud) {
            Eigen::Vector2f end = getMapCoords(point);
            std::vector<Eigen::Vector2i> points;
            bresenham(start[0],  start[1], end[0], end[1], points);

            int n = points.size();
            if(n == 0) {
                continue;
            }

            for(int j = 0; j < n - 1; ++j) {
                int index = getIndex(points[j][0], points[j][1]);
                if(value_[index] + log_odds_free_ >= log_odds_min_) {
                    value_[index] += log_odds_free_;
                    log_odds_[index] += log_odds_free_;
                }
            }

            int index = getIndex(points[n - 1][0], points[n - 1][1]);
            log_odds_[index] += log_odds_occupied_;
            if(value_[index] + log_odds_occupied_ <= log_odds_max_) {
                value_[index] += log_odds_occupied_;
            }
        }
    }
}

void ProbabilityGridMap::bresenham(int x0, int y0, int x1, int y1, std::vector<Eigen::Vector2i>& cells)
{
    int dx = abs(x1 - x0), sx = x0 < x1 ? 1 : -1;
    int dy = abs(y1 - y0), sy = y0 < y1 ? 1 : -1;
    int err = (dx > dy ? dx : -dy) / 2;

    Eigen::Vector2i point(x0, y0);
    cells.push_back(point);

    while (x0 != x1 || y0 != y1) {
        int e2 = err;
        if (e2 > -dx) { err -= dy; x0 += sx; }
        if (e2 <  dy) { err += dx; y0 += sy; }
        point[0] = x0;
        point[1] = y0;
        cells.push_back(point);
    }
}



//double ProbabilityGridMap::calculateEntropyLoss(const std::shared_ptr<LaserScan>& scan, double min_entropy_loss)
//{
//    const PointCloud& point_cloud = scan->getTransformedPointCloud();
//    double mutual_info = 0;
//    for(const Eigen::Vector2f& point : point_cloud) {
//        Eigen::Vector2f end = getMapCoords(point);
//        int index = getIndex(end[0], end[1]);
//        if(log_odds_[index] <=  log_theshold_)
//        {
//            assert(log_odds_[index] > 50);
//            int N = (log_odds_[index] - LogOdds_Unknown);
//            double post_prob = exp(N) / (exp(N) + 1);
//            N = (log_odds_[index] - log_odds_occupied_ - LogOdds_Unknown);
//            double prior_prob = exp(N) / (exp(N) + 1);
//            double info_gain = prior_prob * log(1./prior_prob) + (1. - prior_prob) * log(1./(1-prior_prob))
//                    - post_prob * log(1./post_prob) - (1. - post_prob) * log(1./(1-post_prob));
//            mutual_info += info_gain;
//            if(mutual_info > min_entropy_loss)
//            {
//                return mutual_info;
//            }
//        }
//    }
//    return mutual_info;
//}



float ProbabilityGridMap::calculateEntropyLoss(const std::shared_ptr<LaserScan>& scan, const float min_entropy_loss)
{
    const PointCloud& point_cloud = scan->getTransformedPointCloud();
    float mutual_info = 0;
    for(const Eigen::Vector2f& point : point_cloud) {
        Eigen::Vector2f end = getMapCoords(point);
        int index = getIndex(end[0], end[1]);
        if(log_odds_[index] <=  log_theshold_ && log_odds_[index] > LogOdds_Unknown)
        {
//            assert(log_odds_[index] > 50);
            // N 初始值应该为0, log(P_o / (1 - P_o)) = N
            int N = (log_odds_[index] - LogOdds_Unknown);
            float post_prob = expf(N) / (expf(N) + 1);
            N = (log_odds_[index] - log_odds_occupied_ - LogOdds_Unknown);
            float prior_prob = expf(N) / (expf(N) + 1);
            float info_gain = prior_prob * logf(1.f/prior_prob) + (1.f - prior_prob) * logf(1.f/(1.f-prior_prob))
                               - post_prob * logf(1.f/post_prob) - (1.f - post_prob) * logf(1.f/(1.f-post_prob));
            mutual_info += info_gain;
            if(mutual_info > min_entropy_loss)
            {
                return mutual_info;
            }
        }
    }
    return mutual_info;
}

bool ProbabilityGridMap::calculateCrossEntropy(const std::shared_ptr<LaserScan>& scan, const float cross_entropy_threshold)
{

    Eigen::Vector2f start = getMapCoords(scan->getScanPose());
    const PointCloud& point_cloud = scan->getTransformedPointCloud();
    std::unordered_map<int, int> lookup_tab;

//    auto t1 = std::chrono::steady_clock::now();
    for(const Eigen::Vector2f& point : point_cloud) {
        Eigen::Vector2f end = getMapCoords(point);
        std::vector<Eigen::Vector2i> points;
        bresenham(start[0],  start[1], end[0], end[1], points);

        int n = points.size();
        if(n == 0) {
            continue;
        }

        for(int j = 0; j < n - 1; ++j) {
            int index = getIndex(points[j][0], points[j][1]);
//            if(log_odds_[index] <= LogOdds_Unknown) continue;
            if(lookup_tab.find(index) != lookup_tab.end())
            {
//                if(lookup_tab[index].second + log_odds_free_ >= log_odds_min_) lookup_tab[index].second += log_odds_free_;
                lookup_tab[index] += log_odds_free_ * 8;

            } else
            {
                lookup_tab[index] = log_odds_free_ * 8;
            }
        }

        int index = getIndex(points[n - 1][0], points[n - 1][1]);
//        if(log_odds_[index] > LogOdds_Unknown) continue;
        if(lookup_tab.find(index) != lookup_tab.end())
        {
            lookup_tab[index] += log_odds_occupied_ * 8;
        } else
        {
            lookup_tab[index] = log_odds_occupied_ * 8;
        }
    }

//    auto t2 = std::chrono::steady_clock::now();
//    auto delta_t = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
//    std::cout << "look up table time = " << delta_t.count() * 1000.0 << "ms." << std::endl;


    float total_relative_entropy = 0.f;
    for(const auto& kv:lookup_tab)
    {

        // N 初始值应该为0, log(P_o / (1 - P_o)) = N
        int cur_value = log_odds_[kv.first] < log_odds_max_? log_odds_[kv.first]:log_odds_max_;
        int N = (cur_value - LogOdds_Unknown);
        float prior_prob = 0.f;
        if(N > log_theshold_ - LogOdds_Unknown) prior_prob = 1.0f;
        else if (N > (log_odds_max_ - log_theshold_) - LogOdds_Unknown) prior_prob = expf(N) / (expf(N) + 1);
        N = (cur_value + kv.second - LogOdds_Unknown);
        float post_prob = 0.f;
        if(N > log_theshold_ - LogOdds_Unknown) post_prob = 1.0f;
        else if (N > (log_odds_max_ - log_theshold_) - LogOdds_Unknown) post_prob = expf(N) / (expf(N) + 1);
        if(std::abs(prior_prob - post_prob) < 0.00001f) continue;
        post_prob = std::min(std::max(post_prob, 0.00001f), 0.99999f);
        prior_prob = std::min(std::max(prior_prob, 0.00001f), 0.99999f);
        float relative_entropy =
                post_prob * logf(post_prob / prior_prob) +
                        (1 - post_prob) * logf((1 - post_prob) / (1 - prior_prob));
        total_relative_entropy += relative_entropy;
    }

//    t1 = std::chrono::steady_clock::now();
//    delta_t = std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t2);
//    std::cout << "iterate look up table time = " << delta_t.count() * 1000.0 << "ms." << std::endl;
    std::cout << "total_relative_entropy = " << total_relative_entropy << std::endl;
    if(total_relative_entropy > cross_entropy_threshold)
    {
        for(const auto& kv:lookup_tab)
        {
            log_odds_[kv.first] = log_odds_[kv.first] + kv.second / 8;
            if(log_odds_[kv.first] < log_odds_min_) log_odds_[kv.first] = log_odds_min_;
        }
    }

    return total_relative_entropy > cross_entropy_threshold;
}


//float ProbabilityGridMap::calculateEntropyLoss(const std::shared_ptr<LaserScan>& scan, float min_entropy_loss)
//{
//    const PointCloud& point_cloud = scan->getTransformedPointCloud();
//    float mutual_info = 0;
//    std::unordered_map<int, std::pair<int,int>> lookup_tab;
//    for(const Eigen::Vector2f& point : point_cloud) {
//        Eigen::Vector2f end = getMapCoords(point);
//        int index = getIndex(end[0], end[1]);
//        if(log_odds_[index] > LogOdds_Unknown)
//        {
//            if(lookup_tab.find(index) != lookup_tab.end())
//            {
//                if(lookup_tab[index].second > LogOdds_Unknown)
//                {
//                    lookup_tab[index].second -= log_odds_occupied_;
//                }
//            } else
//            {
//                lookup_tab[index] = std::make_pair(log_odds_[index], log_odds_[index] - log_odds_occupied_);
//            }
//        }
//    }
//    for(const auto& kv:lookup_tab)
//    {
//        if(kv.second.first < log_theshold_)
//        {
//            // N 初始值应该为0, log(P_o / (1 - P_o)) = N
//            int N = (kv.second.first - LogOdds_Unknown);
//            float post_prob = expf(N) / (expf(N) + 1);
//            N = (kv.second.second - LogOdds_Unknown);
//            float prior_prob = expf(N) / (expf(N) + 1);
//            float info_gain = prior_prob * logf(1.f/prior_prob) + (1.f - prior_prob) * logf(1.f/(1.f-prior_prob))
//                               - post_prob * logf(1.f/post_prob) - (1.f - post_prob) * logf(1.f/(1.f-post_prob));
//            mutual_info += info_gain;
//            // 如果mutual_info太大，证明这一帧信息量很大了，不会被选出来丢弃，直接跳过
//            if(mutual_info > min_entropy_loss)
//            {
//                return mutual_info;
//            }
//        }
//    }
//    return mutual_info;
//}

void ProbabilityGridMap::updateProbMap(const std::shared_ptr<LaserScan>& scan)
{
    const PointCloud& point_cloud = scan->getTransformedPointCloud();
    for(const Eigen::Vector2f& point : point_cloud) {
        Eigen::Vector2f end = getMapCoords(point);
        int index = getIndex(end[0], end[1]);
        if(log_odds_[index] - log_odds_occupied_ > log_odds_min_)
        {
            log_odds_[index] -= log_odds_occupied_;
        }
    }
}

// 找出mutual_info最小的一帧，然后丢掉
int ProbabilityGridMap::selectUninformativeScan(const std::vector<std::shared_ptr<LaserScan>>& scans)
{
    float min_mutual_info = std::numeric_limits<float>::max();
    int best_index = 0;
    int index = 0;
    std::shared_ptr<LaserScan> scan_to_discard;
    for(const std::shared_ptr<LaserScan>& scan : scans) {
        float mutual_info = calculateEntropyLoss(scan, min_mutual_info);
        if (mutual_info > 0.) {
            if (mutual_info < min_mutual_info) {
                min_mutual_info = mutual_info;
                best_index = index;
            }
            index++;
        } else {
            return index;
        }
    }
    return best_index;
}


void ProbabilityGridMap::discardDynamicObject(const std::vector<std::shared_ptr<LaserScan>>& scans)
{
    for(const std::shared_ptr<LaserScan>& scan : scans) {
        const PointCloud& point_cloud = scan->getTransformedPointCloud();
        std::vector<int> status;
        status.reserve(point_cloud.size());
        for(const Eigen::Vector2f& point : point_cloud) {
            Eigen::Vector2f end = getMapCoords(point);
            int index = getIndex(end[0], end[1]);
            if(log_odds_[index] > LogOdds_Unknown)
            {
                status.push_back(1);
            } else
            {
                status.push_back(0);
            }
        }
        scan->updatePointCloud(status);
    }
}


} // namespace RLML
