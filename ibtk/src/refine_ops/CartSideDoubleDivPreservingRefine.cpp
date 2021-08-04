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

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <IBTK_config.h>

#include "ibtk/CartSideDoubleDivPreservingRefine.h"
#include "ibtk/IndexUtilities.h"
#include "ibtk/namespaces.h" // IWYU pragma: keep

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/PatchDescriptor.h"
#include "SAMRAI/hier/PatchGeometry.h"
#include "SAMRAI/hier/RefineOperator.h"
#include "SAMRAI/xfer/RefinePatchStrategy.h"
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/pdat/SideGeometry.h"
#include "SAMRAI/pdat/SideIndex.h"
#include "SAMRAI/tbox/Array.h"

#include <cmath>
#include <limits>
#include <string>

// FORTRAN ROUTINES
#if (NDIM == 2)
#define DIV_PRESERVING_CORRECTION_FC IBTK_FC_FUNC_(div_preserving_correction2d, DIV_PRESERVING_CORRECTION2D)
#endif

#if (NDIM == 3)
#define DIV_PRESERVING_CORRECTION_FC IBTK_FC_FUNC_(div_preserving_correction3d, DIV_PRESERVING_CORRECTION3D)
#endif

extern "C"
{
    void DIV_PRESERVING_CORRECTION_FC(double* u0,
                                      double* u1,
#if (NDIM == 3)
                                      double* u2,
#endif
                                      const int& u_gcw,
                                      const int& ilower0,
                                      const int& iupper0,
                                      const int& ilower1,
                                      const int& iupper1,
#if (NDIM == 3)
                                      const int& ilower2,
                                      const int& iupper2,
#endif
                                      const int& correction_box_ilower0,
                                      const int& correction_box_iupper0,
                                      const int& correction_box_ilower1,
                                      const int& correction_box_iupper1,
#if (NDIM == 3)
                                      const int& correction_box_ilower2,
                                      const int& correction_box_iupper2,
#endif
                                      const int* ratio,
                                      const double* dx_fine);
}

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

CartSideDoubleDivPreservingRefine::CartSideDoubleDivPreservingRefine(const int u_dst_idx,
                                                                     const int u_src_idx,
                                                                     const int indicator_idx,
                                                                     std::shared_ptr<RefineOperator > refine_op,
                                                                     std::shared_ptr<CoarsenOperator > coarsen_op,
                                                                     const double fill_time,
                                                                     RefinePatchStrategy* const phys_bdry_op)
    : d_u_dst_idx(u_dst_idx),
      d_u_src_idx(u_src_idx),
      d_indicator_idx(indicator_idx),
      d_fill_time(fill_time),
      d_phys_bdry_op(phys_bdry_op),
      d_refine_op(refine_op),
      d_coarsen_op(coarsen_op)
{
    // intentionally blank
    return;
} // CartSideDoubleDivPreservingRefine

void
CartSideDoubleDivPreservingRefine::setPhysicalBoundaryConditions(Patch& patch,
                                                                 const double fill_time,
                                                                 const IntVector& ghost_width_to_fill)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(MathUtilities<double>::equalEps(fill_time, d_fill_time));
#else
    NULL_USE(d_fill_time);
#endif
    if (d_phys_bdry_op) d_phys_bdry_op->setPhysicalBoundaryConditions(patch, fill_time, ghost_width_to_fill);
    return;
} // setPhysicalBoundaryConditions

IntVector
CartSideDoubleDivPreservingRefine::getRefineOpStencilWidth() const
{
    return REFINE_OP_STENCIL_WIDTH;
} // getRefineOpStencilWidth

void
CartSideDoubleDivPreservingRefine::preprocessRefine(Patch& /*fine*/,
                                                    const Patch& /*coarse*/,
                                                    const Box& /*fine_box*/,
                                                    const IntVector& /*ratio*/)
{
    // intentionally blank
    return;
} // preprocessRefine

void
CartSideDoubleDivPreservingRefine::postprocessRefine(Patch& fine,
                                                     const Patch& coarse,
                                                     const Box& unrestricted_fine_box,
                                                     const IntVector& ratio)
{
    // NOTE: This operator cannot fill the full ghost cell width of the
    // destination data.  We instead restrict the size of the fine box to ensure
    // that we have adequate data to apply the divergence- and curl-preserving
    // corrections.
    const Box fine_box = unrestricted_fine_box * Box::grow(fine.getBox(), 2);

#if !defined(NDEBUG)
    for (int d = 0; d < NDIM; ++d)
    {
        if (ratio(d) % 2 != 0)
        {
            TBOX_ERROR("CartSideDoubleDivPreservingRefine::postprocessRefine():\n"
                       << "  refinement ratio must be a power of 2 for divergence- and "
                          "curl-preserving refinement operator."
                       << std::endl);
        }
    }
#endif

    std::shared_ptr<SideData<double> > fdata = std::static_pointer_cast<SideData<double> >(fine.getPatchData(d_u_dst_idx));
#if !defined(NDEBUG)
    TBOX_ASSERT(fdata);
#endif
    const int fdata_ghosts = fdata->getGhostCellWidth().max();
#if !defined(NDEBUG)
    TBOX_ASSERT(fdata_ghosts == fdata->getGhostCellWidth().min());
#endif
    const int fdata_depth = fdata->getDepth();

    std::shared_ptr<SideData<double> > cdata = std::static_pointer_cast<SideData<double> >(coarse.getPatchData(d_u_dst_idx));
#if !defined(NDEBUG)
    TBOX_ASSERT(cdata);
    const int cdata_ghosts = cdata->getGhostCellWidth().max();
    TBOX_ASSERT(cdata_ghosts == cdata->getGhostCellWidth().min());
    const int cdata_depth = cdata->getDepth();
    TBOX_ASSERT(cdata_depth == fdata_depth);
#endif

    if (ratio == IntVector(Dimension(NDIM), 2))
    {
        // Perform (limited) conservative prolongation of the coarse grid data.
        d_refine_op->refine(fine, coarse, d_u_dst_idx, d_u_dst_idx, fine_box, ratio);

        std::shared_ptr<SideData<double> > u_src_data = std::static_pointer_cast<SideData<double> >(fine.getPatchData(d_u_src_idx));
        std::shared_ptr<SideData<double> > indicator_data = std::static_pointer_cast<SideData<double> >(fine.getPatchData(d_indicator_idx));

        // Ensure that we do not modify any of the data from the old level by
        // setting the value of the fine grid data to equal u_src wherever the
        // indicator data equals "1".
        if (u_src_data && indicator_data)
        {
            for (unsigned int axis = 0; axis < NDIM; ++axis)
            {
                for (Box::Iterator b(SideGeometry::toSideBox(fine_box, axis)); b; b++)
                {
                    const hier::Index& i = b();
                    const SideIndex i_s(i, axis, 0);
                    if (std::abs((*indicator_data)(i_s)-1.0) < 1.0e-12)
                    {
                        for (int depth = 0; depth < fdata_depth; ++depth)
                        {
                            (*fdata)(i_s, depth) = (*u_src_data)(i_s, depth);
                        }
                    }
                }
            }
        }

        // Reinterpolate data in the normal direction in the newly refined part
        // of the level wherever the indicator data does NOT equal "1".  Notice
        // that this loop actually modifies only data that is NOT covered by an
        // overlying coarse grid cell face.
        if (indicator_data)
        {
            for (unsigned int axis = 0; axis < NDIM; ++axis)
            {
                for (Box::Iterator b(SideGeometry::toSideBox(fine_box, axis)); b; b++)
                {
                    const hier::Index& i = b();
                    const SideIndex i_s(i, axis, 0);
                    if (!(std::abs((*indicator_data)(i_s)-1.0) < 1.0e-12))
                    {
                        const hier::Index i_coarse_lower = IndexUtilities::coarsen(i, ratio);
                        const hier::Index i_lower = IndexUtilities::refine(i_coarse_lower, ratio);
                        if (i(axis) == i_lower(axis)) continue;

                        hier::Index i_coarse_upper = i_coarse_lower;
                        i_coarse_upper(axis) += 1;
                        const hier::Index i_upper = IndexUtilities::refine(i_coarse_upper, ratio);

                        const double w1 =
                            static_cast<double>(i(axis) - i_lower(axis)) / static_cast<double>(ratio(axis));
                        const double w0 = 1.0 - w1;

                        const SideIndex i_s_lower(i_lower, axis, 0);
                        const SideIndex i_s_upper(i_upper, axis, 0);
                        for (int depth = 0; depth < fdata_depth; ++depth)
                        {
                            (*fdata)(i_s, depth) = w0 * (*fdata)(i_s_lower, depth) + w1 * (*fdata)(i_s_upper, depth);
                        }
                    }
                }
            }
        }

        // Determine the box on which we need to compute the divergence- and
        // curl-preserving correction.
        const Box correction_box = Box::refine(Box::coarsen(fine_box, 2), 2);
#if !defined(NDEBUG)
        TBOX_ASSERT(fdata->getGhostBox().contains(correction_box));
#endif
        // Apply the divergence- and curl-preserving correction to the fine grid
        // data.
        std::shared_ptr<CartesianPatchGeometry > pgeom_fine = std::static_pointer_cast<CartesianPatchGeometry >(fine.getPatchGeometry());
        const double* const dx_fine = pgeom_fine->getDx();
        for (int d = 0; d < fdata_depth; ++d)
        {
            DIV_PRESERVING_CORRECTION_FC(fdata->getPointer(0, d),
                                         fdata->getPointer(1, d),
#if (NDIM == 3)
                                         fdata->getPointer(2, d),
#endif
                                         fdata_ghosts,
                                         fdata->getBox().lower()(0),
                                         fdata->getBox().upper()(0),
                                         fdata->getBox().lower()(1),
                                         fdata->getBox().upper()(1),
#if (NDIM == 3)
                                         fdata->getBox().lower()(2),
                                         fdata->getBox().upper()(2),
#endif
                                         correction_box.lower()(0),
                                         correction_box.upper()(0),
                                         correction_box.lower()(1),
                                         correction_box.upper()(1),
#if (NDIM == 3)
                                         correction_box.lower()(2),
                                         correction_box.upper()(2),
#endif
                                         ratio,
                                         dx_fine);
        }
    }
    else
    {
        // Setup an intermediate patch.
        const Box intermediate_patch_box = Box::refine(coarse.getBox(), 2);
        Patch intermediate(intermediate_patch_box, coarse.getPatchDescriptor());
        intermediate.allocatePatchData(d_u_dst_idx);

        // Setup a patch geometry object for the intermediate patch.
        std::shared_ptr<CartesianPatchGeometry > pgeom_coarse = std::static_pointer_cast<CartesianPatchGeometry >(coarse.getPatchGeometry());
        const IntVector& ratio_to_level_zero_coarse = pgeom_coarse->getRatio();
        std::vector<Array<bool> > touches_regular_bdry(NDIM), touches_periodic_bdry(NDIM);
        for (int axis = 0; axis < NDIM; ++axis)
        {
            touches_regular_bdry[axis].resizeArray(2);
            touches_periodic_bdry[axis].resizeArray(2);
            for (int upperlower = 0; upperlower < 2; ++upperlower)
            {
                touches_regular_bdry[axis][upperlower] = pgeom_coarse->getTouchesRegularBoundary(axis, upperlower);
                touches_periodic_bdry[axis][upperlower] = pgeom_coarse->getTouchesPeriodicBoundary(axis, upperlower);
            }
        }
        const double* const dx_coarse = pgeom_coarse->getDx();

        const IntVector ratio_to_level_zero_intermediate = ratio_to_level_zero_coarse * 2;
        double dx_intermediate[NDIM], x_lower_intermediate[NDIM], x_upper_intermediate[NDIM];
        for (int d = 0; d < NDIM; ++d)
        {
            dx_intermediate[d] = 0.5 * dx_coarse[d];
            x_lower_intermediate[d] = pgeom_coarse->getXLower()[d];
            x_upper_intermediate[d] = pgeom_coarse->getXUpper()[d];
        }
        intermediate.setPatchGeometry(new CartesianPatchGeometry(ratio_to_level_zero_intermediate,
                                                                       touches_regular_bdry,
                                                                       touches_periodic_bdry,
                                                                       dx_intermediate,
                                                                       x_lower_intermediate,
                                                                       x_upper_intermediate));

        // The intermediate box where we need to fill data must be large enough
        // to provide ghost cell values for the fine fill box.
        const Box intermediate_box = Box::grow(Box::coarsen(fine_box, ratio / 2), 2);

        // Setup the original velocity and indicator data.
        if (fine.checkAllocated(d_u_src_idx) && fine.checkAllocated(d_indicator_idx))
        {
            intermediate.allocatePatchData(d_u_src_idx);
            intermediate.allocatePatchData(d_indicator_idx);
            std::shared_ptr<SideData<double> > u_src_idata = std::static_pointer_cast<SideData<double> >(intermediate.getPatchData(d_u_src_idx));
            std::shared_ptr<SideData<double> > indicator_idata = std::static_pointer_cast<SideData<double> >(intermediate.getPatchData(d_indicator_idx));
            u_src_idata->fillAll(std::numeric_limits<double>::quiet_NaN());
            indicator_idata->fillAll(-1.0);
#if !defined(NDEBUG)
            std::shared_ptr<SideData<double> > u_src_fdata = std::static_pointer_cast<SideData<double> >(fine.getPatchData(d_u_src_idx));
            std::shared_ptr<SideData<double> > indicator_fdata = std::static_pointer_cast<SideData<double> >(fine.getPatchData(d_indicator_idx));
            TBOX_ASSERT(u_src_fdata->getGhostBox().contains(Box::refine(intermediate_box, ratio / 2)));
            TBOX_ASSERT(indicator_fdata->getGhostBox().contains(Box::refine(intermediate_box, ratio / 2)));
#endif
            d_coarsen_op->coarsen(intermediate, fine, d_u_src_idx, d_u_src_idx, intermediate_box, ratio / 2);
            d_coarsen_op->coarsen(intermediate, fine, d_indicator_idx, d_indicator_idx, intermediate_box, ratio / 2);
        }

        // Recursively refine from the coarse patch to the fine patch.
        postprocessRefine(intermediate, coarse, intermediate_box, 2);
        postprocessRefine(fine, intermediate, fine_box, ratio / 2);

        // Deallocate any allocated patch data.
        intermediate.deallocatePatchData(d_u_dst_idx);
        if (fine.checkAllocated(d_u_src_idx) && fine.checkAllocated(d_indicator_idx))
        {
            intermediate.deallocatePatchData(d_u_src_idx);
            intermediate.deallocatePatchData(d_indicator_idx);
        }
    }
    return;
} // postprocessRefine

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
