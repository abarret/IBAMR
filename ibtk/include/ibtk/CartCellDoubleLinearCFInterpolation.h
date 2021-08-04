// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2019 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#ifndef included_IBTK_CartCellDoubleLinearCFInterpolation
#define included_IBTK_CartCellDoubleLinearCFInterpolation

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "ibtk/CoarseFineBoundaryRefinePatchStrategy.h"

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/geom/CartesianCellDoubleLinearRefine.h"
#include "SAMRAI/hier/ComponentSelector.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/RefineOperator.h"


#include <set>
#include <vector>

namespace SAMRAI
{
namespace hier
{

class BoxContainer;

class CoarseFineBoundary;

class Patch;
} // namespace hier
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class CartCellDoubleLinearCFInterpolation is a concrete
 * SAMRAI::xfer::RefinePatchStrategy which sets coarse-fine interface ghost cell
 * values for cell-centered double precision patch data via linear interpolation
 * in the normal and tangential directions at coarse-fine interfaces.
 */
class CartCellDoubleLinearCFInterpolation : public CoarseFineBoundaryRefinePatchStrategy
{
public:
    /*!
     * \brief Constructor.
     */
    CartCellDoubleLinearCFInterpolation() = default;

    /*!
     * \brief Destructor.
     */
    ~CartCellDoubleLinearCFInterpolation();

    /*!
     * \name SAMRAI::xfer::RefinePatchStrategy interface.
     */
    //\{

    /*!
     * Function to set data associated with the given list of patch data indices
     * at patch boundaries that intersect the physical domain boundary.  The
     * patch data components set in this routine correspond to the "scratch"
     * components specified in calls to the registerRefine() function in the
     * SAMRAI::xfer::RefineAlgorithm class.
     *
     * Presently, the implementation does nothing.
     *
     * \param patch                Patch on which to fill boundary data.
     * \param fill_time            Double simulation time for boundary filling.
     * \param ghost_width_to_fill  Integer vector describing maximum ghost width to fill over
     *all
     *registered scratch components.
     */
    void setPhysicalBoundaryConditions(SAMRAI::hier::Patch& patch,
                                       double fill_time,
                                       const SAMRAI::hier::IntVector& ghost_width_to_fill) override;

    /*!
     * Function to return maximum stencil width needed over user-defined data
     * interpolation operations.  This is needed to determine the correct
     * interpolation data dependencies.
     */
    SAMRAI::hier::IntVector getRefineOpStencilWidth(const SAMRAI::tbox::Dimension& dim) const override;

    /*!
     * Function to perform user-defined preprocess data refine operations.  This
     * member function is called before standard refine operations (expressed
     * using concrete subclasses of the SAMRAI::hier::RefineOperator base
     * class).  The preprocess function must refine data from the scratch
     * components of the coarse patch into the scratch components of the fine
     * patch on the specified fine box region.  Recall that the scratch
     * components are specified in calls to the registerRefine() function in the
     * SAMRAI::xfer::RefineAlgorithm class.
     *
     * Presently, the implementation does nothing.
     *
     * \param fine      Fine patch containing destination data.
     * \param coarse    Coarse patch containing source data.
     * \param fine_box  Box region on fine patch into which data is refined.
     * \param ratio     Integer vector containing ratio relating index space between coarse and
     *fine
     *patches.
     */
    void preprocessRefine(SAMRAI::hier::Patch& fine,
                          const SAMRAI::hier::Patch& coarse,
                          const SAMRAI::hier::Box& fine_box,
                          const SAMRAI::hier::IntVector& ratio) override;

    /*!
     * Function to perform user-defined postprocess data refine operations.
     * This member function is called after standard refine operations
     * (expressed using concrete subclasses of the SAMRAI::hier::RefineOperator
     * base class).  The postprocess function must refine data from the scratch
     * components of the coarse patch into the scratch components of the fine
     * patch on the specified fine box region.  Recall that the scratch
     * components are specified in calls to the registerRefine() function in the
     * SAMRAI::xfer::RefineAlgorithm class.
     *
     * \param fine      Fine patch containing destination data.
     * \param coarse    Coarse patch containing source data.
     * \param fine_box  Box region on fine patch into which data is refined.
     * \param ratio     Integer vector containing ratio relating index space between coarse and
     *fine
     *patches.
     */
    void postprocessRefine(SAMRAI::hier::Patch& fine,
                           const SAMRAI::hier::Patch& coarse,
                           const SAMRAI::hier::Box& fine_box,
                           const SAMRAI::hier::IntVector& ratio) override;

    //\}

    /*!
     * \name Extension of SAMRAI::xfer::RefinePatchStrategy interface to support more
     * complex coarse-fine interface discretizations.
     */
    //\{

    /*!
     * Whether or not to employ a consistent interpolation scheme at "Type 2"
     * coarse-fine interface ghost cells.
     */
    void setConsistentInterpolationScheme(bool consistent_type_2_bdry) override;

    /*!
     * \brief Reset the patch data index operated upon by this class.
     */
    void setPatchDataIndex(int patch_data_index) override;

    /*!
     * \brief Reset the patch data indices operated upon by this class.
     */
    void setPatchDataIndices(const std::set<int>& patch_data_indices) override;

    /*!
     * \brief Reset the patch data indices operated upon by this class.
     */
    void setPatchDataIndices(const SAMRAI::hier::ComponentSelector& patch_data_indices) override;

    /*!
     * Set the patch hierarchy used in constructing coarse-fine interface
     * boundary boxes.
     */
    void setPatchHierarchy(std::shared_ptr<SAMRAI::hier::PatchHierarchy > hierarchy) override;

    /*!
     * Clear the patch hierarchy used in constructing coarse-fine interface
     * boundary boxes.
     */
    void clearPatchHierarchy() override;

    /*!
     * Compute the normal extension of fine data at coarse-fine interfaces.
     */
    void computeNormalExtension(SAMRAI::hier::Patch& patch,
                                const SAMRAI::hier::IntVector& ratio,
                                const SAMRAI::hier::IntVector& ghost_width_to_fill) override;

    //\}

protected:
private:
    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    CartCellDoubleLinearCFInterpolation(const CartCellDoubleLinearCFInterpolation& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    CartCellDoubleLinearCFInterpolation& operator=(const CartCellDoubleLinearCFInterpolation& that) = delete;

    /*!
     * The patch data indices corresponding to the "scratch" patch data that is
     * operated on by this class.
     */
    std::set<int> d_patch_data_indices;

    /*!
     * Boolean value indicating whether we are enforcing a consistent
     * interpolation scheme at "Type 2" coarse-fine interface ghost cells.
     */
    bool d_consistent_type_2_bdry = false;

    /*!
     * Cached hierarchy-related information.
     */
    std::shared_ptr<SAMRAI::hier::PatchHierarchy > d_hierarchy;
    std::vector<SAMRAI::hier::CoarseFineBoundary > d_cf_boundary;
    std::vector<SAMRAI::hier::BoxContainer > d_domain_boxes;
    std::vector<SAMRAI::hier::IntVector > d_periodic_shift;
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBTK_CartCellDoubleLinearCFInterpolation
