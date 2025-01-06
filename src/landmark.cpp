//
// Created by ep1 on 11/30/20.
//

#include "landmark.h"



namespace RLML
{
    long unsigned int LandMark::next_id_ = 0;

    LandMark::LandMark(const Eigen::Vector2f& pos, const Eigen::Vector2f& pos_in_scan, const std::shared_ptr<LaserScan>& ref_scan)
            : n_obs_(0), pos_(pos), good_(true), ref_scan_(ref_scan)
    {
        id_ = next_id_++;
        addObservations(ref_scan, pos_in_scan);
    }

    LandMark::LandMark(const long unsigned int id, const Eigen::Vector2f& pos)
            : id_(id), pos_(pos)
    {
    }

    void LandMark::setMapPos(const Eigen::Vector2f& pos)
    {
        pos_ = pos;
    }

    Eigen::Vector2f LandMark::getMapPos() const
    {
        return pos_;
    }

    std::shared_ptr<LaserScan> LandMark::getReferenceScan() const
    {
        return ref_scan_;
    }

    void LandMark::addObservations(const std::shared_ptr<LaserScan>& scan, const Eigen::Vector2f& pos_in_scan)
    {
        if (observations_.count(scan))
        {
            return;
        }

        observations_[scan] = pos_in_scan;
        n_obs_++;
    }

    std::map<std::shared_ptr<LaserScan>, Eigen::Vector2f> LandMark::getObservations() const
    {
        return observations_;
    }

    bool LandMark::getObservation(const std::shared_ptr<LaserScan>& scan, Eigen::Vector2f& pos_in_scan) const
    {
        if (observations_.count(scan))
        {
            pos_in_scan = observations_.at(scan);
            return true;
        }
        else
        {
            return false;
        }
    }

    bool LandMark::getState() const
    {
        return good_;
    }

    void LandMark::setState(const bool good)
    {
        good_ = good;
    }

    size_t LandMark::getId() const
    {
        return id_;
    }

    size_t LandMark::getObsNum() const
    {
        return n_obs_;
    }
}
