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

#include "ibtk/CartCellDoubleQuadraticCFInterpolation.h"
#include "ibtk/namespaces.h" // IWYU pragma: keep

#include "SAMRAI/hier/BoundaryBox.h"
#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/hier/CoarseFineBoundary.h"
#include "SAMRAI/geom/GridGeometry.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/hier/RefineOperator.h"
#include "SAMRAI/tbox/Array.h"

#include <array>
#include <memory>
#include <set>
#include <string>
#include <vector>

// FORTRAN ROUTINES
#if (NDIM == 2)
#define CC_QUAD_TANGENTIAL_INTERPOLATION_FC                                                                            \
    IBTK_FC_FUNC(ccquadtangentialinterpolation2d, CCQUADTANGENTIALINTERPOLATION2D)
#define CC_QUAD_NORMAL_INTERPOLATION_FC IBTK_FC_FUNC(ccquadnormalinterpolation2d, CCQUADNORMALINTERPOLATION2D)
#endif
#if (NDIM == 3)
#define CC_QUAD_TANGENTIAL_INTERPOLATION_FC                                                                            \
    IBTK_FC_FUNC(ccquadtangentialinterpolation3d, CCQUADTANGENTIALINTERPOLATION3D)
#define CC_QUAD_NORMAL_INTERPOLATION_FC IBTK_FC_FUNC(ccquadnormalinterpolation3d, CCQUADNORMALINTERPOLATION3D)
#endif

// Function interfaces
extern "C"
{
    void CC_QUAD_TANGENTIAL_INTERPOLATION_FC(double* U_fine,
                                             const int& U_fine_gcw,
                                             const double* U_coarse,
                                             const int& U_crse_gcw,
                                             const int& ilowerf0,
                                             const int& iupperf0,
                                             const int& ilowerf1,
                                             const int& iupperf1,
#if (NDIM == 3)
                                             const int& ilowerf2,
                                             const int& iupperf2,
#endif
                                             const int& ilowerc0,
                                             const int& iupperc0,
                                             const int& ilowerc1,
                                             const int& iupperc1,
#if (NDIM == 3)
                                             const int& ilowerc2,
                                             const int& iupperc2,
#endif
                                             const int& loc_index,
                                             const int* ratio_to_coarser,
                                             const int* blower,
                                             const int* bupper);

    void CC_QUAD_NORMAL_INTERPOLATION_FC(double* U,
                                         const int& U_gcw,
                                         const int& ilower0,
                                         const int& iupper0,
                                         const int& ilower1,
                                         const int& iupper1,
#if (NDIM == 3)
                                         const int& ilower2,
                                         const int& iupper2,
#endif
                                         const int& loc_index,
                                         const int* ratio_to_coarser,
                                         const int* blower,
                                         const int* bupper);
}

// Note that there are two versions of this code:
//
//    - The expensive version uses only C++ constructs.
//    - The optimized version uses hand-coded Fortran routines.
//
// These two versions of the code may produce different values since they employ
// different treatments at "Type 2" coarse-fine boundary ghost cells.  The
// optimized Fortran version does not presently set values in "Type 0"
// coarse-fine interface ghost cells.
//
// The version of the code to be employed is determined at compile time by the
// flag --enable-expensive-cf-interpolation.

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
static const int REFINE_OP_STENCIL_WIDTH = 1;
static const int GHOST_WIDTH_TO_FILL = 1;

inline int
coarsen(const int& index, const int& ratio)
{
    return (index < 0 ? (index + 1) / ratio - 1 : index / ratio);
} // coarsen

inline hier::Index
coarsen(const hier::Index& index, const IntVector& ratio)
{
    hier::Index coarse_index(Dimension(NDIM));
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        coarse_index(d) = coarsen(index(d), ratio(d));
    }
    return coarse_index;
} // coarsen

inline bool
bdry_boxes_contain_index(const hier::Index& i, const std::vector<const BoundaryBox*>& patch_cf_bdry_boxes)
{
    for (const auto& patch_cf_bdry_box : patch_cf_bdry_boxes)
    {
        const BoundaryBox& bdry_box = *patch_cf_bdry_box;
        if (bdry_box.getBox().contains(i)) return true;
    }
    return false;
} // bdry_boxes_contain_index

inline bool
is_corner_point(const hier::Index& i,
                const unsigned int bdry_normal_axis,
                const bool is_lower,
                const Box& patch_box,
                const std::vector<const BoundaryBox*>& patch_cf_bdry_boxes,
                const IntVector& periodic_shift,
                const BoxContainer& domain_boxes)
{
    // Check to see if the index is adjacent to the patch boundary.  If not, it
    // cannot be a corner point.
    if (is_lower)
    {
        if (i(bdry_normal_axis) + 1 != patch_box.lower()(bdry_normal_axis)) return false;
    }
    else
    {
        if (i(bdry_normal_axis) - 1 != patch_box.upper()(bdry_normal_axis)) return false;
    }

    // Check to see if the adjacent points in the tangential directions are
    // contained in the coarse-fine interface.  If not, the point is a corner
    // point.
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        if (axis != bdry_normal_axis)
        {
            const bool periodic = periodic_shift(axis) > 0;

            hier::Index i_lower(i);
            i_lower(axis) = i(axis) - 1;
            if ((periodic || domain_boxes.contains(i_lower, BlockId(0))) && !bdry_boxes_contain_index(i_lower, patch_cf_bdry_boxes))
                return true;

            hier::Index i_upper(i);
            i_upper(axis) = i(axis) + 1;
            if ((periodic || domain_boxes.contains(i_upper, BlockId(0))) && !bdry_boxes_contain_index(i_upper, patch_cf_bdry_boxes))
                return true;
        }
    }
    return false;
} // is_corner_point
} // namespace

/////////////////////////////// PUBLIC ///////////////////////////////////////

CartCellDoubleQuadraticCFInterpolation::~CartCellDoubleQuadraticCFInterpolation()
{
    clearPatchHierarchy();
    return;
} // ~CartCellDoubleQuadraticCFInterpolation

void
CartCellDoubleQuadraticCFInterpolation::setPhysicalBoundaryConditions(Patch& /*patch*/,
                                                                      const double /*fill_time*/,
                                                                      const IntVector& /*ghost_width_to_fill*/)
{
    // intentionally blank
    return;
} // setPhysicalBoundaryConditions

IntVector
CartCellDoubleQuadraticCFInterpolation::getRefineOpStencilWidth(const SAMRAI::tbox::Dimension& dim) const
{
    return IntVector(dim, REFINE_OP_STENCIL_WIDTH);
} // getRefineOpStencilWidth

void
CartCellDoubleQuadraticCFInterpolation::preprocessRefine(Patch& /*fine*/,
                                                         const Patch& /*coarse*/,
                                                         const Box& /*fine_box*/,
                                                         const IntVector& /*ratio*/)
{
    // intentionally blank
    return;
} // preprocessRefine

void
CartCellDoubleQuadraticCFInterpolation::postprocessRefine(Patch& fine,
                                                          const Patch& coarse,
                                                          const Box& fine_box,
                                                          const IntVector& ratio)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_hierarchy);
#endif
#if !defined(NDEBUG)
    if (fine.inHierarchy())
    {
        // Ensure the fine patch corresponds to the expected patch in the cached
        // patch hierarchy.
        const GlobalId patch_num = fine.getGlobalId();
        const int fine_patch_level_num = fine.getPatchLevelNumber();
        std::shared_ptr<PatchLevel > fine_level = d_hierarchy->getPatchLevel(fine_patch_level_num);
        TBOX_ASSERT(&fine == fine_level->getPatch(patch_num).get());
    }
#endif
    postprocessRefine_optimized(fine, coarse, ratio);
    return;
} // postprocessRefine

void
CartCellDoubleQuadraticCFInterpolation::setConsistentInterpolationScheme(const bool consistent_type_2_bdry)
{
    d_consistent_type_2_bdry = consistent_type_2_bdry;
    return;
} // setConsistentInterpolationScheme

void
CartCellDoubleQuadraticCFInterpolation::setPatchDataIndex(const int patch_data_index)
{
    std::set<int> patch_data_indices;
    patch_data_indices.insert(patch_data_index);
    setPatchDataIndices(patch_data_indices);
    return;
} // setPatchDataIndex

void
CartCellDoubleQuadraticCFInterpolation::setPatchDataIndices(const std::set<int>& patch_data_indices)
{
    d_patch_data_indices.clear();
    d_patch_data_indices = patch_data_indices;
    return;
} // setPatchDataIndices

void
CartCellDoubleQuadraticCFInterpolation::setPatchDataIndices(const ComponentSelector& patch_data_indices)
{
    std::set<int> patch_data_index_set;
    for (int l = 0; l < patch_data_indices.getSize(); ++l)
    {
        if (patch_data_indices.isSet(l))
        {
            const int patch_data_index = l;
            patch_data_index_set.insert(patch_data_index);
        }
    }
    setPatchDataIndices(patch_data_index_set);
    return;
} // setPatchDataIndices

void
CartCellDoubleQuadraticCFInterpolation::setPatchHierarchy(std::shared_ptr<PatchHierarchy > hierarchy)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(hierarchy);
#endif
    if (d_hierarchy) clearPatchHierarchy();
    d_hierarchy = hierarchy;
    const int finest_level_number = d_hierarchy->getFinestLevelNumber();

    d_cf_boundary.resize(finest_level_number + 1, CoarseFineBoundary(Dimension(NDIM)));
    const IntVector& max_ghost_width = getRefineOpStencilWidth(Dimension(NDIM));
    for (int ln = 0; ln <= finest_level_number; ++ln)
    {
        d_cf_boundary[ln] = CoarseFineBoundary(*d_hierarchy, ln, max_ghost_width);
    }

    auto grid_geom = std::static_pointer_cast<GridGeometry>(d_hierarchy->getGridGeometry());
    const BoxContainer& domain_boxes = grid_geom->getPhysicalDomain();

    CoarseFineBoundary bdry1(Dimension(NDIM)), bdry2(Dimension(NDIM));
    bdry1 = bdry2;

    d_domain_boxes.resize(finest_level_number + 1);
    d_periodic_shift.resize(finest_level_number + 1, IntVector::getZero(Dimension(NDIM)));
    for (int ln = 0; ln <= finest_level_number; ++ln)
    {
        std::shared_ptr<PatchLevel > level = d_hierarchy->getPatchLevel(ln);
        const IntVector& ratio = level->getRatioToCoarserLevel();
        d_domain_boxes[ln] = BoxContainer(domain_boxes);
        d_domain_boxes[ln].refine(ratio);
        d_periodic_shift[ln] = grid_geom->getPeriodicShift(ratio);
    }
    return;
} // setPatchHierarchy

void
CartCellDoubleQuadraticCFInterpolation::clearPatchHierarchy()
{
    d_hierarchy = nullptr;
    d_cf_boundary.clear();
    d_domain_boxes.clear();
    d_periodic_shift.clear();
    return;
} // clearPatchHierarchy

void
CartCellDoubleQuadraticCFInterpolation::computeNormalExtension(Patch& patch,
                                                               const IntVector& ratio,
                                                               const IntVector& /*ghost_width_to_fill*/)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_hierarchy);
#endif
    // Ensure that the fine patch is located on the expected destination level;
    // if not, we are not guaranteed to have appropriate coarse-fine interface
    // boundary box information.
    if (!patch.inHierarchy())
    {
        return;
    }
#if !defined(NDEBUG)
    else
    {
        const GlobalId patch_num = patch.getGlobalId();
        const int patch_level_num = patch.getPatchLevelNumber();
        std::shared_ptr<PatchLevel > level = d_hierarchy->getPatchLevel(patch_level_num);
        TBOX_ASSERT(&patch == level->getPatch(patch_num).get());
    }
#endif
    computeNormalExtension_optimized(patch, ratio);
    return;
} // computeNormalExtension

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

void
CartCellDoubleQuadraticCFInterpolation::postprocessRefine_expensive(Patch& fine,
                                                                    const Patch& coarse,
                                                                    const IntVector& ratio)
{
    // Get the cf boundary boxes.
    const GlobalId patch_num = fine.getGlobalId();
    const int fine_patch_level_num = fine.getPatchLevelNumber();
    std::vector<const BoundaryBox*> patch_cf_bdry_boxes;
    {
        const std::vector<BoundaryBox >& cf_bdry_codim1_boxes =
            d_cf_boundary[fine_patch_level_num].getBoundaries(patch_num, 1);
        for (size_t k = 0; k < cf_bdry_codim1_boxes.size(); ++k)
        {
            patch_cf_bdry_boxes.push_back(&(cf_bdry_codim1_boxes[k]));
        }
    }
#if (NDIM > 1)
    {
        const std::vector<BoundaryBox >& cf_bdry_codim2_boxes =
            d_cf_boundary[fine_patch_level_num].getBoundaries(patch_num, 2);
        for (size_t k = 0; k < cf_bdry_codim2_boxes.size(); ++k)
        {
            patch_cf_bdry_boxes.push_back(&(cf_bdry_codim2_boxes[k]));
        }
    }
#if (NDIM > 2)
    {
        const std::vector<BoundaryBox >& cf_bdry_codim3_boxes =
            d_cf_boundary[fine_patch_level_num].getBoundaries(patch_num, 3);
        for (size_t k = 0; k < cf_bdry_codim3_boxes.size(); ++k)
        {
            patch_cf_bdry_boxes.push_back(&(cf_bdry_codim3_boxes[k]));
        }
    }
#endif
#endif

    // Check to see if there are any coarse-fine boundary boxes associated with
    // the patch; if not, there is nothing to do.
    if (patch_cf_bdry_boxes.empty()) return;

    // Get the patch data.
    for (int patch_data_index : d_patch_data_indices)
    {
        auto fdata = std::static_pointer_cast<CellData<double>>(fine.getPatchData(patch_data_index));
        auto cdata = std::static_pointer_cast<CellData<double>>(coarse.getPatchData(patch_data_index));
#if !defined(NDEBUG)
        TBOX_ASSERT(fdata);
        TBOX_ASSERT(cdata);
        TBOX_ASSERT(cdata->getDepth() == fdata->getDepth());
#endif
        const int data_depth = fdata->getDepth();
        const IntVector ghost_width_to_fill(Dimension(NDIM), GHOST_WIDTH_TO_FILL);

        const Box& patch_box_fine = fine.getBox();
        const hier::Index& patch_lower_fine = patch_box_fine.lower();
        auto pgeom_fine = std::static_pointer_cast<CartesianPatchGeometry>(fine.getPatchGeometry());
        const double* const XLower_fine = pgeom_fine->getXLower();
        const double* const dx_fine = pgeom_fine->getDx();

        const Box& patch_box_crse = coarse.getBox();
        const hier::Index& patch_lower_crse = patch_box_crse.lower();
        auto pgeom_crse = std::static_pointer_cast<CartesianPatchGeometry>(coarse.getPatchGeometry());
        const double* const XLower_crse = pgeom_crse->getXLower();
        const double* const dx_crse = pgeom_crse->getDx();

        // Reset values that lie in the coarse-fine interface boundary boxes
        // (located on the *coarse* side of the coarse-fine interface).
        //
        // For co-dimension 1 boundary boxes away from corners in the locally
        // refined grid, we only perform interpolation in the tangential
        // direction at the boundary.  For co-dimension 2 and 3 boundary boxes,
        // and for ghost cells at corners in the co-dimension 1 boundary boxes,
        // we perform coarse interpolation in all directions.
        const BoxContainer& domain_boxes = d_domain_boxes[fine_patch_level_num];
        const IntVector& periodic_shift = d_periodic_shift[fine_patch_level_num];
        for (const auto& patch_cf_bdry_box : patch_cf_bdry_boxes)
        {
            const BoundaryBox& bdry_box = *patch_cf_bdry_box;
            const Box bc_fill_box = pgeom_fine->getBoundaryFillBox(bdry_box, patch_box_fine, ghost_width_to_fill);

            const int bdry_type = bdry_box.getBoundaryType();

            // NOTE: The following values are only used for co-dimension 1
            // boundary boxes.
            const unsigned int location_index = bdry_box.getLocationIndex();
            const unsigned int bdry_normal_axis = location_index / 2;
            const bool is_lower = location_index % 2 == 0;

            for (const auto& i_fine : bc_fill_box)
            {
                const hier::Index i_crse = coarsen(i_fine, ratio);
                const bool corner_point = bdry_type != 1 ? false :
                                                           is_corner_point(i_fine,
                                                                           bdry_normal_axis,
                                                                           is_lower,
                                                                           patch_box_fine,
                                                                           patch_cf_bdry_boxes,
                                                                           periodic_shift,
                                                                           domain_boxes);

                // Determine the interpolation stencil in the coarse index
                // space.
                Box stencil_box_crse(i_crse, i_crse, BlockId::zero());
                stencil_box_crse.grow(IntVector::getOne(Dimension(NDIM)));
                if (bdry_type == 1 && !corner_point)
                {
                    stencil_box_crse.setLower(bdry_normal_axis, i_crse(bdry_normal_axis));
                    stencil_box_crse.setUpper(bdry_normal_axis, i_crse(bdry_normal_axis));
                }

                // Determine the interpolation degrees and weights.
                std::array<int, NDIM> interp_degree;
                for (unsigned int axis = 0; axis < NDIM; ++axis)
                {
                    interp_degree[axis] = stencil_box_crse.upper()(axis) - stencil_box_crse.lower()(axis);
                }
                std::array<std::vector<double>, NDIM> wgts;
                for (unsigned int axis = 0; axis < NDIM; ++axis)
                {
                    const int& degree = interp_degree[axis];
                    const double X = XLower_fine[axis] +
                                     dx_fine[axis] * (static_cast<double>(i_fine(axis) - patch_lower_fine(axis)) + 0.5);
                    std::vector<double> X_crse(degree + 1, 0.0);
                    for (int i_crse = stencil_box_crse.lower()(axis), k = 0; i_crse <= stencil_box_crse.upper()(axis);
                         ++i_crse, ++k)
                    {
                        X_crse[k] = XLower_crse[axis] +
                                    dx_crse[axis] * (static_cast<double>(i_crse - patch_lower_crse(axis)) + 0.5);
                    }
                    wgts[axis].resize(degree + 1, 0.0);
                    switch (degree)
                    {
                    case 0:
                        wgts[axis][0] = 1.0;
                        break;
                    case 1:
                        wgts[axis][0] = (X - X_crse[1]) / (X_crse[0] - X_crse[1]);
                        wgts[axis][1] = (X - X_crse[0]) / (X_crse[1] - X_crse[0]);
                        break;
                    case 2:
                        wgts[axis][0] =
                            ((X - X_crse[1]) * (X - X_crse[2])) / ((X_crse[0] - X_crse[1]) * (X_crse[0] - X_crse[2]));
                        wgts[axis][1] =
                            ((X - X_crse[0]) * (X - X_crse[2])) / ((X_crse[1] - X_crse[0]) * (X_crse[1] - X_crse[2]));
                        wgts[axis][2] =
                            ((X - X_crse[0]) * (X - X_crse[1])) / ((X_crse[2] - X_crse[0]) * (X_crse[2] - X_crse[1]));
                        break;
                    }
                }

                // Interpolate from the coarse grid to the fine grid.
                hier::Index i_intrp(Dimension(NDIM));
                for (int d = 0; d < data_depth; ++d)
                {
                    (*fdata)(CellIndex(i_fine), d) = 0.0;
#if (NDIM > 2)
                    for (int i2 = 0; i2 <= interp_degree[2]; ++i2)
                    {
                        const double& wgt2 = wgts[2][i2];
                        i_intrp(2) = stencil_box_crse.lower()(2) + i2;
#else
                    const double wgt2 = 1.0;
#endif
#if (NDIM > 1)
                        for (int i1 = 0; i1 <= interp_degree[1]; ++i1)
                        {
                            const double& wgt1 = wgts[1][i1];
                            i_intrp(1) = stencil_box_crse.lower()(1) + i1;
#else
                    const double wgt1 = 1.0;
#endif
                            for (int i0 = 0; i0 <= interp_degree[0]; ++i0)
                            {
                                const double& wgt0 = wgts[0][i0];
                                i_intrp(0) = stencil_box_crse.lower()(0) + i0;

                                (*fdata)(CellIndex(i_fine), d) += wgt0 * wgt1 * wgt2 * (*cdata)(CellIndex(i_intrp), d);
                            }
#if (NDIM > 1)
                        }
#endif
#if (NDIM > 2)
                    }
#endif
                }
            }
        }
    }
    return;
} // postprocessRefine_expensive

void
CartCellDoubleQuadraticCFInterpolation::postprocessRefine_optimized(Patch& fine,
                                                                    const Patch& coarse,
                                                                    const IntVector& ratio)
{
    // Get the co-dimension 1 cf boundary boxes.
    const GlobalId patch_num = fine.getGlobalId();
    const int fine_patch_level_num = fine.getPatchLevelNumber();
    const std::vector<BoundaryBox >& cf_bdry_codim1_boxes =
        d_cf_boundary[fine_patch_level_num].getBoundaries(patch_num, 1);
    if (cf_bdry_codim1_boxes.size() == 0) return;

    // Get the patch data.
    for (const auto& patch_data_index : d_patch_data_indices)
    {
        auto fdata = std::static_pointer_cast<CellData<double>>(fine.getPatchData(patch_data_index));
        auto cdata = std::static_pointer_cast<CellData<double>>(coarse.getPatchData(patch_data_index));
#if !defined(NDEBUG)
        TBOX_ASSERT(fdata);
        TBOX_ASSERT(cdata);
        TBOX_ASSERT(cdata->getDepth() == fdata->getDepth());
#endif
        const int U_fine_ghosts = (fdata->getGhostCellWidth()).max();
        const int U_crse_ghosts = (cdata->getGhostCellWidth()).max();
#if !defined(NDEBUG)
        if (U_fine_ghosts != (fdata->getGhostCellWidth()).min())
        {
            TBOX_ERROR("CartCellDoubleQuadraticCFInterpolation::postprocessRefine():\n"
                       << "   patch data does not have uniform ghost cell widths" << std::endl);
        }
        if (U_crse_ghosts != (cdata->getGhostCellWidth()).min())
        {
            TBOX_ERROR("CartCellDoubleQuadraticCFInterpolation::postprocessRefine():\n"
                       << "   patch data does not have uniform ghost cell widths" << std::endl);
        }
#endif
        const int data_depth = fdata->getDepth();
        const IntVector ghost_width_to_fill(Dimension(NDIM), GHOST_WIDTH_TO_FILL);
        auto pgeom_fine = std::static_pointer_cast<CartesianPatchGeometry>(fine.getPatchGeometry());
        const Box& patch_box_fine = fine.getBox();
        const Box& patch_box_crse = coarse.getBox();
        for (size_t k = 0; k < cf_bdry_codim1_boxes.size(); ++k)
        {
            const BoundaryBox& bdry_box = cf_bdry_codim1_boxes[k];
            const Box bc_fill_box = pgeom_fine->getBoundaryFillBox(bdry_box, patch_box_fine, ghost_width_to_fill);
            const unsigned int location_index = bdry_box.getLocationIndex();
            for (int depth = 0; depth < data_depth; ++depth)
            {
                double* const U_fine = fdata->getPointer(depth);
                const double* const U_crse = cdata->getPointer(depth);
                CC_QUAD_TANGENTIAL_INTERPOLATION_FC(U_fine,
                                                    U_fine_ghosts,
                                                    U_crse,
                                                    U_crse_ghosts,
                                                    patch_box_fine.lower(0),
                                                    patch_box_fine.upper(0),
                                                    patch_box_fine.lower(1),
                                                    patch_box_fine.upper(1),
#if (NDIM == 3)
                                                    patch_box_fine.lower(2),
                                                    patch_box_fine.upper(2),
#endif
                                                    patch_box_crse.lower(0),
                                                    patch_box_crse.upper(0),
                                                    patch_box_crse.lower(1),
                                                    patch_box_crse.upper(1),
#if (NDIM == 3)
                                                    patch_box_crse.lower(2),
                                                    patch_box_crse.upper(2),
#endif
                                                    location_index,
                                                    &(ratio(0)),
                                                    &(bc_fill_box.lower()(0)),
                                                    &(bc_fill_box.upper()(0)));
            }
        }
    }
    return;
} // postprocessRefine_optimized

void
CartCellDoubleQuadraticCFInterpolation::computeNormalExtension_expensive(Patch& patch,
                                                                         const IntVector& ratio,
                                                                         const IntVector& ghost_width_to_fill)
{
    // Get the co-dimension 1 cf boundary boxes.
    const GlobalId patch_num = patch.getGlobalId();
    const int patch_level_num = patch.getPatchLevelNumber();
    const std::vector<BoundaryBox >& cf_bdry_codim1_boxes = d_cf_boundary[patch_level_num].getBoundaries(patch_num, 1);
    const int n_cf_bdry_codim1_boxes = cf_bdry_codim1_boxes.size();

    // Check to see if there are any co-dimension 1 coarse-fine boundary boxes
    // associated with the patch; if not, there is nothing to do.
    if (n_cf_bdry_codim1_boxes == 0) return;

    // Collect pointers to all of the cf boundary boxes.
    std::vector<const BoundaryBox*> patch_cf_bdry_boxes;
    {
        const std::vector<BoundaryBox >& cf_bdry_codim1_boxes =
            d_cf_boundary[patch_level_num].getBoundaries(patch_num, 1);
        for (size_t k = 0; k < cf_bdry_codim1_boxes.size(); ++k)
        {
            patch_cf_bdry_boxes.push_back(&(cf_bdry_codim1_boxes[k]));
        }
    }
#if (NDIM > 1)
    {
        const std::vector<BoundaryBox >& cf_bdry_codim2_boxes =
            d_cf_boundary[patch_level_num].getBoundaries(patch_num, 2);
        for (size_t k = 0; k < cf_bdry_codim2_boxes.size(); ++k)
        {
            patch_cf_bdry_boxes.push_back(&(cf_bdry_codim2_boxes[k]));
        }
    }
#if (NDIM > 2)
    {
        const std::vector<BoundaryBox >& cf_bdry_codim3_boxes =
            d_cf_boundary[patch_level_num].getBoundaries(patch_num, 3);
        for (size_t k = 0; k < cf_bdry_codim3_boxes.size(); ++k)
        {
            patch_cf_bdry_boxes.push_back(cf_bdry_codim3_boxes.getPointer(k));
        }
    }
#endif
#endif

    // Get the patch data.
    for (const auto& patch_data_index : d_patch_data_indices)
    {
        auto data = std::static_pointer_cast<CellData<double>>(patch.getPatchData(patch_data_index));
#if !defined(NDEBUG)
        TBOX_ASSERT(data);
#endif
        const int data_depth = data->getDepth();

        const Box& patch_box = patch.getBox();
        const hier::Index& patch_lower = patch_box.lower();
        const hier::Index& patch_upper = patch_box.upper();
        auto pgeom = std::static_pointer_cast<CartesianPatchGeometry>(patch.getPatchGeometry());

        // Use quadratic interpolation in the normal direction to reset values
        // that lie in the co-dimension 1 coarse-fine interface boundary boxes
        // (located on the *coarse* side of the coarse-fine interface).
        const BoxContainer& domain_boxes = d_domain_boxes[patch_level_num];
        const IntVector& periodic_shift = d_periodic_shift[patch_level_num];
        for (int k = 0; k < n_cf_bdry_codim1_boxes; ++k)
        {
            const BoundaryBox& bdry_box = cf_bdry_codim1_boxes[k];
            const Box bc_fill_box = pgeom->getBoundaryFillBox(bdry_box, patch_box, ghost_width_to_fill);
            const unsigned int location_index = bdry_box.getLocationIndex();
            const unsigned int bdry_normal_axis = location_index / 2;
            const bool is_lower = location_index % 2 == 0;
            for (const auto& i_bdry : bc_fill_box)
            {
                if (!is_corner_point(i_bdry,
                                     bdry_normal_axis,
                                     is_lower,
                                     patch_box,
                                     patch_cf_bdry_boxes,
                                     periodic_shift,
                                     domain_boxes))
                {
                    hier::Index i_intr0(i_bdry), i_intr1(i_bdry);
                    if (is_lower)
                    {
                        i_intr0(bdry_normal_axis) = patch_lower(bdry_normal_axis);
                        i_intr1(bdry_normal_axis) = patch_lower(bdry_normal_axis) + 1;
                    }
                    else
                    {
                        i_intr0(bdry_normal_axis) = patch_upper(bdry_normal_axis);
                        i_intr1(bdry_normal_axis) = patch_upper(bdry_normal_axis) - 1;
                    }

                    const auto r = static_cast<double>(ratio(bdry_normal_axis));
                    const int i = (is_lower ? i_bdry(bdry_normal_axis) - patch_lower(bdry_normal_axis) :
                                              i_bdry(bdry_normal_axis) - patch_upper(bdry_normal_axis));

                    const double X = (is_lower ? static_cast<double>(i) + 0.5 :
                                                 static_cast<double>(i) - 0.5); // the position of the c-f boundary cell
                    const double X0 = (is_lower ? +0.5 : -0.5);         // the position of the first  interior cell
                    const double X1 = (is_lower ? +1.5 : -1.5);         // the position of the second interior cell
                    const double X2 = (is_lower ? -0.5 * r : +0.5 * r); // the position of the
                                                                        // coarse-grid cell

                    const double wgt0 = (X - X1) * (X - X2) / ((X0 - X1) * (X0 - X2));
                    const double wgt1 = (X - X0) * (X - X2) / ((X1 - X0) * (X1 - X2));
                    const double wgt2 = (X - X0) * (X - X1) / ((X2 - X0) * (X2 - X1));

                    for (int d = 0; d < data_depth; ++d)
                    {
                        (*data)(CellIndex(i_bdry), d) =
                            wgt0 * (*data)(CellIndex(i_intr0), d) + wgt1 * (*data)(CellIndex(i_intr1), d) + wgt2 * (*data)(CellIndex(i_bdry), d);
                    }
                }
            }
        }
    }
    return;
} // computeNormalExtension_expensive

void
CartCellDoubleQuadraticCFInterpolation::computeNormalExtension_optimized(Patch& patch,
                                                                         const IntVector& ratio)
{
    // Get the co-dimension 1 cf boundary boxes.
    const GlobalId patch_num = patch.getGlobalId();
    const int patch_level_num = patch.getPatchLevelNumber();
    const std::vector<BoundaryBox >& cf_bdry_codim1_boxes = d_cf_boundary[patch_level_num].getBoundaries(patch_num, 1);
    const int n_cf_bdry_codim1_boxes = cf_bdry_codim1_boxes.size();

    // Check to see if there are any co-dimension 1 coarse-fine boundary boxes
    // associated with the patch; if not, there is nothing to do.
    if (n_cf_bdry_codim1_boxes == 0) return;

    // Get the patch data.
    for (int patch_data_index : d_patch_data_indices)
    {
        std::shared_ptr<CellData<double> > data = std::static_pointer_cast<CellData<double> >(patch.getPatchData(patch_data_index));
#if !defined(NDEBUG)
        TBOX_ASSERT(data);
#endif
        const int U_ghosts = (data->getGhostCellWidth()).max();
#if !defined(NDEBUG)
        if (U_ghosts != (data->getGhostCellWidth()).min())
        {
            TBOX_ERROR("CartCellDoubleQuadraticCFInterpolation::computeNormalExtension():\n"
                       << "   patch data does not have uniform ghost cell widths" << std::endl);
        }
#endif
        const int data_depth = data->getDepth();
        const IntVector ghost_width_to_fill(Dimension(NDIM), GHOST_WIDTH_TO_FILL);
        std::shared_ptr<CartesianPatchGeometry > pgeom = std::static_pointer_cast<CartesianPatchGeometry >(patch.getPatchGeometry());
        const Box& patch_box = patch.getBox();
        for (int k = 0; k < n_cf_bdry_codim1_boxes; ++k)
        {
            const BoundaryBox& bdry_box = cf_bdry_codim1_boxes[k];
            const Box bc_fill_box = pgeom->getBoundaryFillBox(bdry_box, patch_box, ghost_width_to_fill);
            const unsigned int location_index = bdry_box.getLocationIndex();
            for (int depth = 0; depth < data_depth; ++depth)
            {
                double* const U = data->getPointer(depth);
                CC_QUAD_NORMAL_INTERPOLATION_FC(U,
                                                U_ghosts,
                                                patch_box.lower(0),
                                                patch_box.upper(0),
                                                patch_box.lower(1),
                                                patch_box.upper(1),
#if (NDIM == 3)
                                                patch_box.lower(2),
                                                patch_box.upper(2),
#endif
                                                location_index,
                                                &(ratio(0)),
                                                &(bc_fill_box.lower()(0)),
                                                &(bc_fill_box.upper()(0)));
            }
        }
    }
    return;
} // computeNormalExtension_optimized

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
