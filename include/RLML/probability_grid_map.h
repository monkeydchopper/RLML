#ifndef PROBABILITY_GRID_MAP_H
#define PROBABILITY_GRID_MAP_H

#include <chrono>
#include <memory>
#include <unordered_map>
#include "grid_map.h"
#include "laser_scan.h"

namespace RLML
{

const uint8_t LogOdds_Unknown = 50;


class ProbabilityGridMap : public GridMap<uint8_t>
{
public:
    ProbabilityGridMap(int width, int height, float resolution) :
        GridMap(width, height, resolution), log_odds_occupied_(2), log_odds_free_(-1)
    {
        log_odds_min_ = 0;
        log_odds_max_ = 100;
        log_theshold_ = 65;
        log_odds_ = new int16_t[size_];
        for(int i = 0; i < size_; ++i) {
            value_[i] = LogOdds_Unknown;
            log_odds_[i] = LogOdds_Unknown;
        }
    }
    ~ProbabilityGridMap() {
        if(log_odds_) {
            delete[] log_odds_;
        }    }

    void createFromScan(const std::vector<std::shared_ptr<LaserScan>>& scans);
//    void updateMap(const std::shared_ptr<LaserScan>& scan);
    bool calculateCrossEntropy(const std::shared_ptr<LaserScan>& scan, const float cross_entropy_threshold);
    float calculateEntropyLoss(const std::shared_ptr<LaserScan>& scan, const float min_entropy_loss);
    void updateProbMap(const std::shared_ptr<LaserScan>& scan);
    int selectUninformativeScan(const std::vector<std::shared_ptr<LaserScan>>& scans);
    void discardDynamicObject(const std::vector<std::shared_ptr<LaserScan>>& scans);

private:
    void bresenham(int x0, int y0, int x1, int y1, std::vector<Eigen::Vector2i>& cells);

private:
    int8_t log_odds_occupied_;
    int8_t log_odds_free_;
    int8_t log_odds_min_;
    int8_t log_odds_max_;

    // for the usage of lifelong mapping, it can be larger than 100 and smaller than 0
    int16_t* log_odds_;
    int8_t log_theshold_;
};

} // namespace RLML

#endif // PROBABILITY_GRID_MAP_H
