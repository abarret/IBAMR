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

#include "ibtk/PhysicalBoundaryUtilities.h"
#include "ibtk/StaggeredPhysicalBoundaryHelper.h"
#include "ibtk/namespaces.h" // IWYU pragma: keep

#include "SAMRAI/pdat/ArrayData.h"
#include "SAMRAI/hier/BoundaryBox.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/PatchGeometry.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/solv/RobinBcCoefStrategy.h"
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/hier/Variable.h"
#include "SAMRAI/tbox/Array.h"
#include "SAMRAI/tbox/MathUtilities.h"

#include <algorithm>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>
#include "SAMRAI/pdat/CellVariable.h"

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

void
StaggeredPhysicalBoundaryHelper::copyDataAtDirichletBoundaries(const int u_out_data_idx,
                                                               const int u_in_data_idx,
                                                               const int coarsest_ln,
                                                               const int finest_ln) const
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_hierarchy);
#endif
    const int finest_hier_level = d_hierarchy->getFinestLevelNumber();
    for (int ln = (coarsest_ln == -1 ? 0 : coarsest_ln); ln <= (finest_ln == -1 ? finest_hier_level : finest_ln); ++ln)
    {
        std::shared_ptr<PatchLevel > level = d_hierarchy->getPatchLevel(ln);
        for (PatchLevel::Iterator p = level->begin(); p != level->end(); p++)
        {
            std::shared_ptr<Patch > patch = *p;
            if (patch->getPatchGeometry()->getTouchesRegularBoundary())
            {
                auto u_out_data = std::static_pointer_cast<SideData<double> >(patch->getPatchData(u_out_data_idx));
                auto u_in_data = std::static_pointer_cast<SideData<double> >(patch->getPatchData(u_in_data_idx));
                copyDataAtDirichletBoundaries(u_out_data, u_in_data, patch);
            }
        }
    }
    return;
} // copyDataAtDirichletBoundaries

void
StaggeredPhysicalBoundaryHelper::copyDataAtDirichletBoundaries(std::shared_ptr<SideData<double> > u_out_data,
                                                               std::shared_ptr<SideData<double> > u_in_data,
                                                               std::shared_ptr<Patch > patch) const
{
    if (!patch->getPatchGeometry()->getTouchesRegularBoundary()) return;
    const int ln = patch->getPatchLevelNumber();
    const int& patch_num = patch->getLocalId().getValue();
    const std::vector<BoundaryBox >& physical_codim1_boxes = d_physical_codim1_boxes[ln].find(patch_num)->second;
    const int n_physical_codim1_boxes = physical_codim1_boxes.size();
    const std::vector<std::shared_ptr<ArrayData<int> > >& dirichlet_bdry_locs =
        d_dirichlet_bdry_locs[ln].find(patch_num)->second;
    for (int n = 0; n < n_physical_codim1_boxes; ++n)
    {
        const int bdry_normal_axis = physical_codim1_boxes[n].getLocationIndex() / 2;
        const ArrayData<int>& bdry_locs_data = *dirichlet_bdry_locs[n];
        for (const auto& i : bdry_locs_data.getBox())
        {
            if (bdry_locs_data(i, 0) == 1)
            {
                (*u_out_data)(SideIndex(i, bdry_normal_axis, SideIndex::Lower)) =
                    (*u_in_data)(SideIndex(i, bdry_normal_axis, SideIndex::Lower));
            }
        }
    }
    return;
} // copyDataAtDirichletBoundaries

void
StaggeredPhysicalBoundaryHelper::setupMaskingFunction(const int mask_data_idx,
                                                      const int coarsest_ln,
                                                      const int finest_ln) const
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_hierarchy);
#endif
    const int finest_hier_level = d_hierarchy->getFinestLevelNumber();
    for (int ln = (coarsest_ln == -1 ? 0 : coarsest_ln); ln <= (finest_ln == -1 ? finest_hier_level : finest_ln); ++ln)
    {
        std::shared_ptr<PatchLevel > level = d_hierarchy->getPatchLevel(ln);
        for (PatchLevel::Iterator p = level->begin(); p != level->end(); p++)
        {
            std::shared_ptr<Patch > patch = *p;
            auto mask_data = std::static_pointer_cast<SideData<int> >(patch->getPatchData(mask_data_idx));
            if (patch->getPatchGeometry()->getTouchesRegularBoundary())
            {
                setupMaskingFunction(mask_data, patch);
            }
            else
            {
                mask_data->fillAll(0);
            }
        }
    }
    return;
} // setupMaskingFunction

void
StaggeredPhysicalBoundaryHelper::setupMaskingFunction(std::shared_ptr<SideData<int> > mask_data,
                                                      std::shared_ptr<Patch > patch) const
{
    mask_data->fillAll(0);
    if (patch->getPatchGeometry()->getTouchesRegularBoundary()) return;
    const int ln = patch->getPatchLevelNumber();
    const int& patch_num = patch->getLocalId().getValue();
    const std::vector<BoundaryBox >& physical_codim1_boxes = d_physical_codim1_boxes[ln].find(patch_num)->second;
    const int n_physical_codim1_boxes = physical_codim1_boxes.size();
    const std::vector<std::shared_ptr<ArrayData<int> > >& dirichlet_bdry_locs =
        d_dirichlet_bdry_locs[ln].find(patch_num)->second;
    for (int n = 0; n < n_physical_codim1_boxes; ++n)
    {
        const int bdry_normal_axis = physical_codim1_boxes[n].getLocationIndex() / 2;
        const ArrayData<int>& bdry_locs_data = *dirichlet_bdry_locs[n];
        for (const auto& i : bdry_locs_data.getBox())
        {
            if (bdry_locs_data(i, 0) == 1) (*mask_data)(SideIndex(i, bdry_normal_axis, SideIndex::Lower)) = 1;
        }
    }
    return;
} // setupMaskingFunction

bool
StaggeredPhysicalBoundaryHelper::patchTouchesDirichletBoundary(std::shared_ptr<Patch > patch) const
{
    if (!patch->getPatchGeometry()->getTouchesRegularBoundary()) return false;
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        if (patchTouchesDirichletBoundaryAxis(patch, axis)) return true;
    }
    return false;
} // patchTouchesDirichletBoundary

bool
StaggeredPhysicalBoundaryHelper::patchTouchesDirichletBoundaryAxis(std::shared_ptr<Patch > patch,
                                                                   const unsigned int axis) const
{
    if (!patch->getPatchGeometry()->getTouchesRegularBoundary()) return false;
    const int ln = patch->getPatchLevelNumber();
    const int patch_num = patch->getLocalId().getValue();
    const std::vector<BoundaryBox >& physical_codim1_boxes = d_physical_codim1_boxes[ln].find(patch_num)->second;
    const int n_physical_codim1_boxes = physical_codim1_boxes.size();
    const std::vector<std::shared_ptr<ArrayData<int> > >& dirichlet_bdry_locs =
        d_dirichlet_bdry_locs[ln].find(patch_num)->second;
    for (int n = 0; n < n_physical_codim1_boxes; ++n)
    {
        const unsigned int bdry_normal_axis = physical_codim1_boxes[n].getLocationIndex() / 2;
        if (bdry_normal_axis == axis)
        {
            const ArrayData<int>& bdry_locs_data = *dirichlet_bdry_locs[n];
            for (const auto& i : bdry_locs_data.getBox())
            {
                if (bdry_locs_data(i, 0) == 1) return true;
            }
        }
    }
    return false;
} // patchTouchesDirichletBoundaryAxis

void
StaggeredPhysicalBoundaryHelper::cacheBcCoefData(const std::vector<RobinBcCoefStrategy*>& u_bc_coefs,
                                                 const double fill_time,
                                                 const std::shared_ptr<PatchHierarchy > hierarchy)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(u_bc_coefs.size() == NDIM);
    TBOX_ASSERT(hierarchy);
#endif
    if (d_hierarchy) clearBcCoefData();

    // Cache boundary values.
    d_hierarchy = hierarchy;
    const int finest_hier_level = d_hierarchy->getFinestLevelNumber();
    d_physical_codim1_boxes.resize(finest_hier_level + 1);
    d_dirichlet_bdry_locs.resize(finest_hier_level + 1);
    for (int ln = 0; ln <= finest_hier_level; ++ln)
    {
        std::shared_ptr<PatchLevel > level = d_hierarchy->getPatchLevel(ln);
        for (PatchLevel::Iterator p = level->begin(); p != level->end(); p++)
        {
            std::shared_ptr<Patch > patch = *p;
            const int patch_num = patch->getLocalId().getValue();
            auto pgeom = std::static_pointer_cast<CartesianPatchGeometry >(patch->getPatchGeometry());
            if (pgeom->getTouchesRegularBoundary())
            {
                std::vector<BoundaryBox >& physical_codim1_boxes = d_physical_codim1_boxes[ln][patch_num];
                physical_codim1_boxes = PhysicalBoundaryUtilities::getPhysicalBoundaryCodim1Boxes(*patch);
                const int n_physical_codim1_boxes = physical_codim1_boxes.size();
                std::vector<std::shared_ptr<ArrayData<int> > >& dirichlet_bdry_locs =
                    d_dirichlet_bdry_locs[ln][patch_num];
                dirichlet_bdry_locs.resize(n_physical_codim1_boxes);
                Box bc_coef_box(Dimension(NDIM));
                BoundaryBox trimmed_bdry_box(Dimension(NDIM));
                for (int n = 0; n < n_physical_codim1_boxes; ++n)
                {
                    const BoundaryBox& bdry_box = physical_codim1_boxes[n];
                    StaggeredPhysicalBoundaryHelper::setupBcCoefBoxes(bc_coef_box, trimmed_bdry_box, bdry_box, patch);
                    const unsigned int bdry_normal_axis = bdry_box.getLocationIndex() / 2;
                    auto acoef_data = std::make_shared<ArrayData<double>>(bc_coef_box, 1);
                    auto bcoef_data = std::make_shared<ArrayData<double>>(bc_coef_box, 1);
                    std::shared_ptr<ArrayData<double> > gcoef_data;
                    auto dummy_variable = std::make_shared<CellVariable<double>>(Dimension(NDIM), "dummy");
                    u_bc_coefs[bdry_normal_axis]->setBcCoefs(acoef_data,
                                                             bcoef_data,
                                                             gcoef_data,
                                                             dummy_variable,
                                                             *patch,
                                                             trimmed_bdry_box,
                                                             fill_time);
                    dirichlet_bdry_locs[n] = std::make_shared<ArrayData<int>>(bc_coef_box, 1);
                    std::shared_ptr<ArrayData<int>>& bdry_locs_data = dirichlet_bdry_locs[n];
                    auto temp = std::make_shared<ArrayData<double>>(bc_coef_box, 1);
                    for (const auto& i : bc_coef_box)
                    {
                        const double& alpha = (*acoef_data)(i, 0);
                        const double& beta = (*bcoef_data)(i, 0);
#if !defined(NDEBUG)
                        TBOX_ASSERT(MathUtilities<double>::equalEps(alpha + beta, 1.0));
                        TBOX_ASSERT(MathUtilities<double>::equalEps(alpha, 1.0) ||
                                    MathUtilities<double>::equalEps(beta, 1.0));
#endif
                        (*bdry_locs_data)(i, 0) = (MathUtilities<double>::equalEps(alpha, 1.0) &&
                                (beta == 0.0 || MathUtilities<double>::equalEps(beta, 0.0))) ? 1 : 0;
                    }
                }
            }
            else
            {
                d_physical_codim1_boxes[ln][patch_num].resize(0, BoundaryBox(Dimension(NDIM)));
                d_dirichlet_bdry_locs[ln][patch_num].clear();
            }
        }
    }
    return;
} // cacheBcCoefData

void
StaggeredPhysicalBoundaryHelper::clearBcCoefData()
{
    d_hierarchy = nullptr;
    d_physical_codim1_boxes.clear();
    d_dirichlet_bdry_locs.clear();
    return;
} // clearBcCoefData

/////////////////////////////// PROTECTED ////////////////////////////////////

void
StaggeredPhysicalBoundaryHelper::setupBcCoefBoxes(Box& bc_coef_box,
                                                  BoundaryBox& trimmed_bdry_box,
                                                  const BoundaryBox& bdry_box,
                                                  std::shared_ptr<Patch > patch)
{
    std::shared_ptr<PatchGeometry > pgeom = std::static_pointer_cast<PatchGeometry >(patch->getPatchGeometry());
    const Box& patch_box = patch->getBox();
    const unsigned int location_index = bdry_box.getLocationIndex();
    const unsigned int bdry_normal_axis = location_index / 2;
    Box bc_fill_box = pgeom->getBoundaryFillBox(bdry_box, patch_box, /* gcw_to_fill */ IntVector(Dimension(NDIM), 1));
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        if (d != bdry_normal_axis)
        {
            bc_fill_box.setLower(d, std::max(bc_fill_box.lower(d), patch_box.lower(d)));
            bc_fill_box.setUpper(d, std::min(bc_fill_box.upper(d), patch_box.upper(d)));
        }
    }
    trimmed_bdry_box = BoundaryBox(bdry_box.getBox() * bc_fill_box, /* codimension */ 1, location_index);
    bc_coef_box = PhysicalBoundaryUtilities::makeSideBoundaryCodim1Box(trimmed_bdry_box);
    return;
} // setupBcCoefBoxes

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
