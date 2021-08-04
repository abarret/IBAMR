// ---------------------------------------------------------------------
//
// Copyright (c) 2011 - 2019 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#ifndef included_IBTK_FaceDataSynchronization
#define included_IBTK_FaceDataSynchronization

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "ibtk/ibtk_utilities.h"

#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/xfer/CoarsenAlgorithm.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/xfer/RefineAlgorithm.h"



#include <string>
#include <vector>

namespace SAMRAI
{
namespace xfer
{

class CoarsenSchedule;

class RefineSchedule;
} // namespace xfer
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class FaceDataSynchronization encapsulates the operations required to
 * "synchronize" face-centered values defined at patch boundaries.
 */
class FaceDataSynchronization : public 
{
public:
    /*!
     * \brief Class FaceDataSynchronization::SynchronizationTransactionComponent
     * encapsulates options for filling ghost cell values via class
     * FaceDataSynchronization.
     */
    class SynchronizationTransactionComponent
    {
    public:
        friend class FaceDataSynchronization;

        /*!
         * \brief Default constructor.
         */
        inline SynchronizationTransactionComponent(int data_idx = -1, const std::string& coarsen_op_name = "NONE")
            : d_data_idx(data_idx), d_coarsen_op_name(coarsen_op_name)
        {
            // intentionally blank
            return;
        } // SynchronizationTransactionComponent

        /*!
         * \brief Copy constructor.
         *
         * \param from The value to copy to this object.
         */
        inline SynchronizationTransactionComponent(const SynchronizationTransactionComponent& from)
            : d_data_idx(from.d_data_idx), d_coarsen_op_name(from.d_coarsen_op_name)
        {
            // intentionally blank
            return;
        } // SynchronizationTransactionComponent

        /*!
         * \brief Assignment operator.
         *
         * \param that The value to assign to this object.
         *
         * \return A reference to this object.
         */
        inline SynchronizationTransactionComponent& operator=(const SynchronizationTransactionComponent& that)
        {
            if (this != &that)
            {
                d_data_idx = that.d_data_idx;
                d_coarsen_op_name = that.d_coarsen_op_name;
            }
            return *this;
        } // operator=

        /*!
         * \brief Destructor.
         */
        inline ~SynchronizationTransactionComponent()
        {
            // intentionally blank
            return;
        } // ~SynchronizationTransactionComponent

    private:
        // Data.
        int d_data_idx;
        std::string d_coarsen_op_name;
    };

    /*!
     * \brief Default constructor.
     */
    FaceDataSynchronization() = default;

    /*!
     * \brief Destructor.
     */
    virtual ~FaceDataSynchronization();

    /*!
     * \brief Setup the hierarchy synchronization operator to perform the
     * specified synchronization transactions on the specified patch hierarchy.
     */
    void initializeOperatorState(const SynchronizationTransactionComponent& transaction_comp,
                                 std::shared_ptr<SAMRAI::hier::PatchHierarchy > hierarchy);

    /*!
     * \brief Setup the hierarchy synchronization operator to perform the
     * specified collection of synchronization transactions on the specified
     * patch hierarchy.
     */
    void initializeOperatorState(const std::vector<SynchronizationTransactionComponent>& transaction_comps,
                                 std::shared_ptr<SAMRAI::hier::PatchHierarchy > hierarchy);

    /*!
     * \brief Reset transaction component with the synchronization operator.
     */
    void resetTransactionComponent(const SynchronizationTransactionComponent& transaction_comps);

    /*!
     * \brief Reset transaction components with the synchronization operator.
     */
    void resetTransactionComponents(const std::vector<SynchronizationTransactionComponent>& transaction_comps);

    /*!
     * \brief Clear all cached data.
     */
    void deallocateOperatorState();

    /*!
     * \brief Synchronize the data on all levels of the patch hierarchy.
     */
    void synchronizeData(double fill_time);

protected:
private:
    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    FaceDataSynchronization(const FaceDataSynchronization& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    FaceDataSynchronization& operator=(const FaceDataSynchronization& that) = delete;

    // Boolean indicating whether the operator is initialized.
    bool d_is_initialized = false;

    // The component synchronization operations to perform.
    std::vector<SynchronizationTransactionComponent> d_transaction_comps;

    // Hierarchy configuration.
    std::shared_ptr<SAMRAI::hier::PatchHierarchy > d_hierarchy;
    std::shared_ptr<SAMRAI::geom::CartesianGridGeometry > d_grid_geom;
    int d_coarsest_ln = IBTK::invalid_level_number, d_finest_ln = IBTK::invalid_level_number;

    // Cached communications algorithms and schedules.
    std::shared_ptr<SAMRAI::xfer::CoarsenAlgorithm > d_coarsen_alg;
    std::vector<std::shared_ptr<SAMRAI::xfer::CoarsenSchedule > > d_coarsen_scheds;

    std::shared_ptr<SAMRAI::xfer::RefineAlgorithm > d_refine_alg;
    std::vector<std::shared_ptr<SAMRAI::xfer::RefineSchedule > > d_refine_scheds;
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBTK_FaceDataSynchronization
