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

#include "ibtk/CartCellDoubleLinearCFInterpolation.h"
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

#include <memory>
#include <set>
#include <string>
#include <vector>

// FORTRAN ROUTINES
#if (NDIM == 2)
#define CC_LINEAR_NORMAL_INTERPOLATION_FC IBTK_FC_FUNC(cclinearnormalinterpolation2d, CCLINEARNORMALINTERPOLATION2D)
#endif
#if (NDIM == 3)
#define CC_LINEAR_NORMAL_INTERPOLATION_FC IBTK_FC_FUNC(cclinearnormalinterpolation3d, CCLINEARNORMALINTERPOLATION3D)
#endif

// Function interfaces
extern "C"
{
    void CC_LINEAR_NORMAL_INTERPOLATION_FC(double* U,
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
                                           const int& ratio,
                                           const int* blower,
                                           const int* bupper);
}

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
static const int REFINE_OP_STENCIL_WIDTH = 1;
static const int GHOST_WIDTH_TO_FILL = 1;
} // namespace

/////////////////////////////// PUBLIC ///////////////////////////////////////

CartCellDoubleLinearCFInterpolation::~CartCellDoubleLinearCFInterpolation()
{
    clearPatchHierarchy();
    return;
} // ~CartCellDoubleLinearCFInterpolation

void
CartCellDoubleLinearCFInterpolation::setPhysicalBoundaryConditions(Patch& /*patch*/,
                                                                   const double /*fill_time*/,
                                                                   const IntVector& /*ghost_width_to_fill*/)
{
    // intentionally blank
    return;
} // setPhysicalBoundaryConditions

IntVector
CartCellDoubleLinearCFInterpolation::getRefineOpStencilWidth(const Dimension& dim) const
{
    return IntVector(dim, REFINE_OP_STENCIL_WIDTH);
} // getRefineOpStencilWidth

void
CartCellDoubleLinearCFInterpolation::preprocessRefine(Patch& /*fine*/,
                                                      const Patch& /*coarse*/,
                                                      const Box& /*fine_box*/,
                                                      const IntVector& /*ratio*/)
{
    // intentionally blank
    return;
} // preprocessRefine

void
CartCellDoubleLinearCFInterpolation::postprocessRefine(Patch& /*fine*/,
                                                       const Patch& coarse,
                                                       const Box& fine_box,
                                                       const IntVector& ratio)
{
    // intentionally blank
    return;
} // postprocessRefine

void
CartCellDoubleLinearCFInterpolation::setConsistentInterpolationScheme(const bool consistent_type_2_bdry)
{
    d_consistent_type_2_bdry = consistent_type_2_bdry;
    return;
} // setConsistentInterpolationScheme

void
CartCellDoubleLinearCFInterpolation::setPatchDataIndex(const int patch_data_index)
{
    std::set<int> patch_data_indices;
    patch_data_indices.insert(patch_data_index);
    setPatchDataIndices(patch_data_indices);
    return;
} // setPatchDataIndex

void
CartCellDoubleLinearCFInterpolation::setPatchDataIndices(const std::set<int>& patch_data_indices)
{
    d_patch_data_indices.clear();
    d_patch_data_indices = patch_data_indices;
    return;
} // setPatchDataIndices

void
CartCellDoubleLinearCFInterpolation::setPatchDataIndices(const ComponentSelector& patch_data_indices)
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
CartCellDoubleLinearCFInterpolation::setPatchHierarchy(std::shared_ptr<PatchHierarchy > hierarchy)
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

    std::shared_ptr<GridGeometry > grid_geom = std::static_pointer_cast<GridGeometry>(d_hierarchy->getGridGeometry());
    const BoxContainer& domain_boxes = grid_geom->getPhysicalDomain();

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
CartCellDoubleLinearCFInterpolation::clearPatchHierarchy()
{
    d_hierarchy = nullptr;
    d_cf_boundary.clear();
    d_domain_boxes.clear();
    d_periodic_shift.clear();
    return;
} // clearPatchHierarchy

void
CartCellDoubleLinearCFInterpolation::computeNormalExtension(Patch& patch,
                                                            const IntVector& ratio,
                                                            const IntVector& /*ghost_width_to_fill*/)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_hierarchy);
    TBOX_ASSERT(!d_consistent_type_2_bdry);
    TBOX_ASSERT(ratio.min() == ratio.max());
#endif
    // Ensure that the fine patch is located on the expected destination level;
    // if not, we are not guaranteed to have appropriate coarse-fine interface
    // boundary box information.
    if (!patch.inHierarchy()) return;

    // Get the co-dimension 1 cf boundary boxes.
    const GlobalId& patch_num = patch.getGlobalId();
    const int patch_level_num = patch.getPatchLevelNumber();
#if !defined(NDEBUG)
    std::shared_ptr<PatchLevel > level = d_hierarchy->getPatchLevel(patch_level_num);
    TBOX_ASSERT(&patch == level->getPatch(patch_num).get());
#endif
    const std::vector<BoundaryBox >& cf_bdry_codim1_boxes = d_cf_boundary[patch_level_num].getBoundaries(patch.getGlobalId(), 1);
    const int n_cf_bdry_codim1_boxes = cf_bdry_codim1_boxes.size();

    // Check to see if there are any co-dimension 1 coarse-fine boundary boxes
    // associated with the patch; if not, there is nothing to do.
    if (n_cf_bdry_codim1_boxes == 0) return;

    // Get the patch data.
    for (const auto& patch_data_index : d_patch_data_indices)
    {
        auto data = std::static_pointer_cast<CellData<double>>(patch.getPatchData(patch_data_index));
#if !defined(NDEBUG)
        TBOX_ASSERT(data);
#endif
        const int U_ghosts = (data->getGhostCellWidth()).max();
#if !defined(NDEBUG)
        if (U_ghosts != (data->getGhostCellWidth()).min())
        {
            TBOX_ERROR("CartCellDoubleLinearCFInterpolation::computeNormalExtension():\n"
                       << "   patch data does not have uniform ghost cell widths" << std::endl);
        }
#endif
        const int data_depth = data->getDepth();
        const IntVector ghost_width_to_fill(Dimension(NDIM), GHOST_WIDTH_TO_FILL);
        auto pgeom = std::static_pointer_cast<CartesianPatchGeometry>(patch.getPatchGeometry());
        const Box& patch_box = patch.getBox();
        for (int k = 0; k < n_cf_bdry_codim1_boxes; ++k)
        {
            const BoundaryBox& bdry_box = cf_bdry_codim1_boxes[k];
            const Box bc_fill_box = pgeom->getBoundaryFillBox(bdry_box, patch_box, ghost_width_to_fill);
            const unsigned int location_index = bdry_box.getLocationIndex();
            for (int depth = 0; depth < data_depth; ++depth)
            {
                double* const U = data->getPointer(depth);
                const int r = ratio.min();
                CC_LINEAR_NORMAL_INTERPOLATION_FC(U,
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
                                                  r,
                                                  &(bc_fill_box.lower()(0)),
                                                  &(bc_fill_box.upper()(0)));
            }
        }
    }
    return;
} // computeNormalExtension

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
