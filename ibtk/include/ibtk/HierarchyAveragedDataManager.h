// ---------------------------------------------------------------------
//
// Copyright (c) 2022 - 2022 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

/////////////////////////////// INCLUDE GUARD ////////////////////////////////

#ifndef included_IBTK_HierarchyAveragedDataManager
#define included_IBTK_HierarchyAveragedDataManager

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/SnapshotCache.h>
#include <ibtk/ibtk_utilities.h>

#include "CellVariable.h"
#include "HierarchySideDataOpsReal.h"
#include "IntVector.h"
#include "SideVariable.h"
#include "VisItDataWriter.h"
#include "tbox/Database.h"
#include "tbox/Pointer.h"

#include <string>
#include <vector>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class HierarchyTimeInterpolator provides a method of determining and storing average fields over periodic
 * intervals. This class can also be used to determine intermediate velocity fields when at a periodic steady state.
 */
class HierarchyAveragedDataManager
{
public:
    /*!
     * The constructor for class HierarchyTimeInterpolator sets some default values and determines data centering. The
     * expected period and number of snapshots must be set in the input database. This class assumes the snapshots will
     * be taken at equidistant points along the periodic interval.
     */
    HierarchyAveragedDataManager(std::string object_name,
                                 SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > var,
                                 SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                                 SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                                 SAMRAI::tbox::Pointer<SAMRAI::hier::GridGeometry<NDIM> > grid_geom,
                                 bool register_for_restart = true);

    /*!
     * The constructor for class HierarchyTimeInterpolator sets some default values and determines data centering. In
     * this constructor, the times at which snapshots are taken are set, as well as the period length.
     */
    HierarchyAveragedDataManager(std::string object_name,
                                 SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > var,
                                 SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                                 SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                                 std::set<double> snapshot_time_points,
                                 double t_start,
                                 double t_end,
                                 SAMRAI::tbox::Pointer<SAMRAI::hier::GridGeometry<NDIM> > grid_geom,
                                 bool register_for_restart = true);

    /*!
     * \brief The destructor for class HierarchyTimeInterpolator deallocates patch data as needed.
     */
    ~HierarchyAveragedDataManager() = default;

    /*!
     * Clear the snapshots taken with this class.
     */
    void clearSnapshots();

    /*!
     * Update the time averaged snapshot. The update is computed through the relation
     * ``` u_avg = u_avg + 1/N * (u - u_avg)```
     * in which N is the number of updates to the snapshot. If this the first update after clearSnapshots() is called,
     * we simply copy data.
     *
     * This function returns true if the snapshot is at a steady state, or that 1/N*||u - u_avg|| < threshold. Note that
     * a weight patch index should be supplied to accurately compute the norm.
     */
    inline bool updateTimeAveragedSnapshot(int u_idx,
                                           double time,
                                           SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                                           const int wgt_idx = IBTK::invalid_index,
                                           double tol = 1.0e-8)
    {
        return updateTimeAveragedSnapshot(u_idx, time, hierarchy, d_mean_refine_type, wgt_idx, tol);
    }
    bool updateTimeAveragedSnapshot(int u_idx,
                                    double time,
                                    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                                    const std::string& mean_refine_type,
                                    const int wgt_idx = IBTK::invalid_index,
                                    double tol = 1.0e-8);

    /*!
     * Return true if all the tracked mean fields are at a steady state.
     */
    inline bool isAtPeriodicSteadyState()
    {
        bool steady_state = true;
        for (const auto& idx_steady_state : d_idx_steady_state_map)
        {
            steady_state = steady_state && idx_steady_state.second;
        }
        return steady_state;
    }

    /*!
     * Return whether the point at the specified time is at a periodic steady state.
     */
    inline bool isAtPeriodicSteadyState(double time, const double tol)
    {
        time = getTimePt(time, tol);
        return d_idx_steady_state_map.at(time);
    }

    /*!
     * Return the time points this class is storing.
     */
    inline const std::set<double>& getSnapshotTimePts()
    {
        return d_snapshot_time_pts;
    }

    /*!
     * Returns the index of the snapshot, if that time point is stored. Otherwise returns -1 (the largest unsigned
     * integer).
     */
    double getTimePt(double time, double tol);

    /*!
     * Get the SnapshotCache object
     */
    const std::unique_ptr<SnapshotCache>& getSnapshotCache()
    {
        return d_snapshot_cache;
    }

    /*!
     * Set the threshold for achieving a steady state.
     */
    void setSteadyStateThreshold(const double threshold)
    {
        d_periodic_thresh = threshold;
    }

private:
    /*!
     * \brief Registers a scratch variable with the variable database. This data is allocated and deallocated as needed.
     */
    void commonConstructor(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                           SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy);
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    HierarchyAveragedDataManager() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    HierarchyAveragedDataManager(const HierarchyAveragedDataManager& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    HierarchyAveragedDataManager& operator=(const HierarchyAveragedDataManager& that) = delete;

    std::string d_object_name;

    /*
     * Data for tracking mean flow quantities.
     */
    unsigned int d_num_averaging_cycles = 0;
    SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > d_var;
    int d_scratch_idx = IBTK::invalid_index;
    std::string d_mean_refine_type = "CONSERVATIVE_LINEAR_REFINE";

    /*
     * Interval to keep flow quantities
     */
    double d_t_start = std::numeric_limits<double>::quiet_NaN(), d_t_end = std::numeric_limits<double>::quiet_NaN(),
           d_t_period = std::numeric_limits<double>::quiet_NaN(),
           d_periodic_thresh = std::numeric_limits<double>::quiet_NaN();

    /*
     * Map between time point and the number of updates. We use this value to update average snapshots.
     *
     * Note we have to be careful when retrieving values from this map. Values should ONLY be retrieved and stored using
     * the values stored in d_snapshot_time_pts.
     */
    std::map<double, unsigned int> d_idx_num_updates_map;

    /*
     * Map between time point and whether it's at a steady state.
     *
     * Note we have to be careful when retrieving values from this map. Values should ONLY be retrieved and stored using
     * the values stored in d_snapshot_time_pts.
     */
    std::map<double, bool> d_idx_steady_state_map;

    /*
     * Set of time points that snapshots are determined. These values are between d_t_start and d_t_end.
     */
    std::set<double> d_snapshot_time_pts;

    std::unique_ptr<SnapshotCache> d_snapshot_cache;

    SAMRAI::tbox::Pointer<SAMRAI::math::HierarchyDataOpsReal<NDIM, double> > d_hier_data_ops;

    // Drawing stuff
    std::unique_ptr<SAMRAI::appu::VisItDataWriter<NDIM> > d_visit_data_writer;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_mean_var;
    int d_mean_idx = IBTK::invalid_index;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_dev_var;
    int d_dev_idx = IBTK::invalid_index;
    int d_visit_ts = 0;
    bool d_output_data = true;
    bool d_enable_logging = true;
};

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBTK_HierarchyAveragedDataManager