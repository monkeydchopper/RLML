//
// Created by ep1 on 11/30/20.
//

#ifndef RLML_LANDMARK_H
#define RLML_LANDMARK_H

#include "laser_scan.h"
#include "map"
#include <memory>

namespace RLML
{

    class LaserScan;


    class LandMark
    {
    public:
        LandMark(const Eigen::Vector2f& pos, const Eigen::Vector2f& pos_in_scan, const std::shared_ptr<LaserScan>& refScan);
        LandMark(const long unsigned int id, const Eigen::Vector2f& pos);

        void setMapPos(const Eigen::Vector2f& pos);
        Eigen::Vector2f getMapPos() const;
        std::shared_ptr<LaserScan> getReferenceScan() const;

        void addObservations(const std::shared_ptr<LaserScan>& scan, const Eigen::Vector2f& pos_in_scan);
        std::map<std::shared_ptr<LaserScan>, Eigen::Vector2f> getObservations() const;
        bool getObservation(const std::shared_ptr<LaserScan>& scan, Eigen::Vector2f& pos_in_scan) const;

        bool getState() const;
        void setState(const bool good);

        size_t getId() const;
        size_t getObsNum() const;

    protected:
        static size_t next_id_;
        size_t id_;
        size_t n_obs_;
        Eigen::Vector2f pos_;
        bool good_;
        std::shared_ptr<LaserScan> ref_scan_;
        std::map<std::shared_ptr<LaserScan>, Eigen::Vector2f> observations_;

    };
}




#endif //RLML_LANDMARK_H
