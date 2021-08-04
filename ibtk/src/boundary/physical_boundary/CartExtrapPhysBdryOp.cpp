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

#include "ibtk/CartExtrapPhysBdryOp.h"
#include "ibtk/PhysicalBoundaryUtilities.h"
#include "ibtk/namespaces.h" // IWYU pragma: keep

#include "SAMRAI/hier/BoundaryBox.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/hier/ComponentSelector.h"
#include "SAMRAI/pdat/FaceData.h"
#include "SAMRAI/pdat/FaceIndex.h"
#include "SAMRAI/pdat/FaceIterator.h"
#include "SAMRAI/pdat/FaceVariable.h"
#include "SAMRAI/pdat/NodeData.h"
#include "SAMRAI/pdat/NodeIndex.h"
#include "SAMRAI/pdat/NodeIterator.h"
#include "SAMRAI/pdat/NodeVariable.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/PatchGeometry.h"
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/pdat/SideIndex.h"
#include "SAMRAI/pdat/SideIterator.h"
#include "SAMRAI/pdat/SideVariable.h"
#include "SAMRAI/hier/Variable.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/tbox/Array.h"

#include "SAMRAI/tbox/Utilities.h"

#include <array>
#include <ostream>
#include <set>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
static const int REFINE_OP_STENCIL_WIDTH = 0;

template <typename D, typename I>
inline double
compute_linear_extrap(D& patch_data, const I& i, const I& i_intr, const IntVector& i_shft, const int depth)
{
    double ret_val = patch_data(i_intr, depth);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        if (i_shft(d) != 0)
        {
            const I& i_intr0 = i_intr;
            I i_intr1 = i_intr;
            i_intr1(d) += i_shft(d);

            const double& f0 = patch_data(i_intr0, depth);
            const double& f1 = patch_data(i_intr1, depth);

            const double du = f0 - f1;
            const double delta = std::abs(i(d) - i_intr(d));

            ret_val += du * delta;
        }
    }
    return ret_val;
} // compute_linear_extrap

template <typename D, typename I>
inline double
compute_quadratic_extrap(D& patch_data,
                         const I& i,
                         const I& i_intr,
                         const IntVector& i_shft,
                         const int depth,
                         const int codim)
{
    if (codim == 1)
    {
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            if (i_shft(d) != 0)
            {
                const I& i_intr0 = i_intr;

                I i_intr1 = i_intr;
                i_intr1(d) += i_shft(d);

                I i_intr2 = i_intr1;
                i_intr2(d) += i_shft(d);

                const double& f0 = patch_data(i_intr0, depth);
                const double& f1 = patch_data(i_intr1, depth);
                const double& f2 = patch_data(i_intr2, depth);

                const double x = std::abs(i(d) - i_intr(d));

                return (1.0 / 2.0 * f2 - f1 + 1.0 / 2.0 * f0) * x * x +
                       (1.0 / 2.0 * f2 - 2.0 * f1 + 3.0 / 2.0 * f0) * x + f0;
            }
        }
    }
    else
    {
        return compute_linear_extrap(patch_data, i, i_intr, i_shft, depth);
    }
    return 0.0; // this statement should not be reached
} // compute_quadratic_extrap
} // namespace

/////////////////////////////// PUBLIC ///////////////////////////////////////

CartExtrapPhysBdryOp::CartExtrapPhysBdryOp(const int patch_data_index, const std::string& extrap_type)
{
    setPatchDataIndex(patch_data_index);
    setExtrapolationType(extrap_type);
    return;
} // CartExtrapPhysBdryOp

CartExtrapPhysBdryOp::CartExtrapPhysBdryOp(const std::set<int>& patch_data_indices, const std::string& extrap_type)
{
    setPatchDataIndices(patch_data_indices);
    setExtrapolationType(extrap_type);
    return;
} // CartExtrapPhysBdryOp

CartExtrapPhysBdryOp::CartExtrapPhysBdryOp(const ComponentSelector& patch_data_indices, const std::string& extrap_type)
{
    setPatchDataIndices(patch_data_indices);
    setExtrapolationType(extrap_type);
    return;
} // CartExtrapPhysBdryOp

void
CartExtrapPhysBdryOp::setPatchDataIndex(const int patch_data_index)
{
    std::set<int> patch_data_indices;
    patch_data_indices.insert(patch_data_index);
    setPatchDataIndices(patch_data_indices);
    return;
} // setPatchDataIndex

void
CartExtrapPhysBdryOp::setPatchDataIndices(const std::set<int>& patch_data_indices)
{
    d_patch_data_indices.clear();
    d_patch_data_indices = patch_data_indices;
    return;
} // setPatchDataIndices

void
CartExtrapPhysBdryOp::setPatchDataIndices(const ComponentSelector& patch_data_indices)
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
CartExtrapPhysBdryOp::setExtrapolationType(const std::string& extrap_type)
{
    // Ensure that the extrapolation type is supported by this class.
    if (extrap_type != "CONSTANT" && extrap_type != "LINEAR" && extrap_type != "QUADRATIC")
    {
        TBOX_ERROR("CartExtrapPhysBdryOp::setExtrapolationType():\n"
                   << "  unknown extrapolation type: " << extrap_type << "\n"
                   << "  valid selections are: CONSTANT, LINEAR, or QUADRATIC" << std::endl);
    }
    d_extrap_type = extrap_type;
    return;
} // setExtrapolationType

void
CartExtrapPhysBdryOp::setPhysicalBoundaryConditions(Patch& patch,
                                                    const double /*fill_time*/,
                                                    const IntVector& ghost_width_to_fill)
{
    if (ghost_width_to_fill == IntVector(Dimension(NDIM), 0)) return;

    std::shared_ptr<PatchGeometry > pgeom = std::static_pointer_cast<PatchGeometry >(patch.getPatchGeometry());
    const Box& patch_box = patch.getBox();

    std::vector<std::pair<Box, std::pair<int, int> > > bdry_fill_boxes;

#if (NDIM > 1)
#if (NDIM > 2)
    // Compute the co-dimension three boundary fill boxes.
    const std::vector<BoundaryBox > physical_codim3_boxes =
        PhysicalBoundaryUtilities::getPhysicalBoundaryCodim3Boxes(patch);
    const int n_physical_codim3_boxes = physical_codim3_boxes.size();
    for (int n = 0; n < n_physical_codim3_boxes; ++n)
    {
        const BoundaryBox& bdry_box = physical_codim3_boxes[n];
        const Box bdry_fill_box = pgeom->getBoundaryFillBox(bdry_box, patch_box, ghost_width_to_fill);
        const unsigned int location_index = bdry_box.getLocationIndex();
        const int codim = 3;
        bdry_fill_boxes.push_back(std::make_pair(bdry_fill_box, std::make_pair(location_index, codim)));
    }
#endif
    // Compute the co-dimension two boundary fill boxes.
    const std::vector<BoundaryBox > physical_codim2_boxes =
        PhysicalBoundaryUtilities::getPhysicalBoundaryCodim2Boxes(patch);
    const int n_physical_codim2_boxes = physical_codim2_boxes.size();
    for (int n = 0; n < n_physical_codim2_boxes; ++n)
    {
        const BoundaryBox& bdry_box = physical_codim2_boxes[n];
        const Box bdry_fill_box = pgeom->getBoundaryFillBox(bdry_box, patch_box, ghost_width_to_fill);
        const unsigned int location_index = bdry_box.getLocationIndex();
        const int codim = 2;
        bdry_fill_boxes.push_back(std::make_pair(bdry_fill_box, std::make_pair(location_index, codim)));
    }
#endif
    // Compute the co-dimension one boundary fill boxes.
    const std::vector<BoundaryBox > physical_codim1_boxes =
        PhysicalBoundaryUtilities::getPhysicalBoundaryCodim1Boxes(patch);
    const int n_physical_codim1_boxes = physical_codim1_boxes.size();
    for (int n = 0; n < n_physical_codim1_boxes; ++n)
    {
        const BoundaryBox& bdry_box = physical_codim1_boxes[n];
        const Box bdry_fill_box = pgeom->getBoundaryFillBox(bdry_box, patch_box, ghost_width_to_fill);
        const unsigned int location_index = bdry_box.getLocationIndex();
        const int codim = 1;
        bdry_fill_boxes.push_back(std::make_pair(bdry_fill_box, std::make_pair(location_index, codim)));
    }

    // Set the boundary values.
    setPhysicalBoundaryConditions_cell(patch, bdry_fill_boxes);
    setPhysicalBoundaryConditions_face(patch, bdry_fill_boxes);
    setPhysicalBoundaryConditions_node(patch, bdry_fill_boxes);
    setPhysicalBoundaryConditions_side(patch, bdry_fill_boxes);
    return;
} // setPhysicalBoundaryConditions

IntVector
CartExtrapPhysBdryOp::getRefineOpStencilWidth(const SAMRAI::tbox::Dimension& dim) const
{
    return IntVector(dim, REFINE_OP_STENCIL_WIDTH);
} // getRefineOpStencilWidth

void
CartExtrapPhysBdryOp::preprocessRefine(Patch& /*fine*/,
                                       const Patch& /*coarse*/,
                                       const Box& /*fine_box*/,
                                       const IntVector& /*ratio*/)
{
    // intentionally blank
    return;
} // preprocessRefine

void
CartExtrapPhysBdryOp::postprocessRefine(Patch& /*fine*/,
                                        const Patch& /*coarse*/,
                                        const Box& /*fine_box*/,
                                        const IntVector& /*ratio*/)
{
    // intentionally blank
    return;
} // postprocessRefine

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

void
CartExtrapPhysBdryOp::setPhysicalBoundaryConditions_cell(
    Patch& patch,
    const std::vector<std::pair<Box, std::pair<int, int> > >& bdry_fill_boxes)
{
    const Box& patch_box = patch.getBox();
    const hier::Index& patch_lower = patch_box.lower();
    const hier::Index& patch_upper = patch_box.upper();

    const int extrap_type =
        (d_extrap_type == "CONSTANT" ? 0 : (d_extrap_type == "LINEAR" ? 1 : (d_extrap_type == "QUADRATIC" ? 2 : -1)));
    if (extrap_type != 0 && extrap_type != 1 && extrap_type != 2)
    {
        TBOX_ERROR("CartExtrapPhysBdryOp::setPhysicalBoundaryConditions():\n"
                   << "  unknown extrapolation type: " << d_extrap_type << "\n"
                   << "  valid selections are: CONSTANT, LINEAR, or QUADRATIC" << std::endl);
    }

    // Set the physical boundary conditions for the specified patch data
    // indices.
    for (const auto& patch_data_idx : d_patch_data_indices)
    {
        VariableDatabase* var_db = VariableDatabase::getDatabase();
        std::shared_ptr<Variable > var;
        var_db->mapIndexToVariable(patch_data_idx, var);
        auto cc_var = std::static_pointer_cast<CellVariable<double>>(var);
        if (!cc_var) continue;

        auto patch_data = std::static_pointer_cast<CellData<double> >(patch.getPatchData(patch_data_idx));
        const Box& ghost_box = patch_data->getGhostBox();

        // Loop over the boundary fill boxes and extrapolate the data.
        for (const auto& bdry_fill_box_pair : bdry_fill_boxes)
        {
            const Box& bdry_fill_box = bdry_fill_box_pair.first;
            const unsigned int location_index = bdry_fill_box_pair.second.first;
            const int codim = bdry_fill_box_pair.second.second;
#if (NDIM == 2)
            const std::array<bool, NDIM> is_lower = { { PhysicalBoundaryUtilities::isLower(location_index, codim, 0),
                                                        PhysicalBoundaryUtilities::isLower(
                                                            location_index, codim, 1) } };
            const std::array<bool, NDIM> is_upper = { { PhysicalBoundaryUtilities::isUpper(location_index, codim, 0),
                                                        PhysicalBoundaryUtilities::isUpper(
                                                            location_index, codim, 1) } };
#endif
#if (NDIM == 3)
            const std::array<bool, NDIM> is_lower = { { PhysicalBoundaryUtilities::isLower(location_index, codim, 0),
                                                        PhysicalBoundaryUtilities::isLower(location_index, codim, 1),
                                                        PhysicalBoundaryUtilities::isLower(
                                                            location_index, codim, 2) } };
            const std::array<bool, NDIM> is_upper = { { PhysicalBoundaryUtilities::isUpper(location_index, codim, 0),
                                                        PhysicalBoundaryUtilities::isUpper(location_index, codim, 1),
                                                        PhysicalBoundaryUtilities::isUpper(
                                                            location_index, codim, 2) } };
#endif
            // Loop over the boundary box indices and compute the nearest
            // interior index.
            for (int depth = 0; depth < patch_data->getDepth(); ++depth)
            {
                Box intersect = bdry_fill_box * ghost_box;
                CellIterator end(CellGeometry::end(intersect));
                for (CellIterator b = CellGeometry::begin(intersect); b != end; b++)
                {
                    const CellIndex& i = *b;
                    CellIndex i_intr = i;
                    IntVector i_shft = IntVector::getZero(Dimension(NDIM));
                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        if (is_lower[d])
                        {
                            i_intr(d) = patch_lower(d);
                            i_shft(d) = +1; // use interior data for extrapolation
                        }
                        else if (is_upper[d])
                        {
                            i_intr(d) = patch_upper(d);
                            i_shft(d) = -1; // use interior data for extrapolation
                        }
                    }

                    // Perform constant, linear, or quadratic extrapolation.
                    switch (extrap_type)
                    {
                    case 0:
                        (*patch_data)(i, depth) = (*patch_data)(i_intr, depth);
                        break;
                    case 1:
                        (*patch_data)(i, depth) = compute_linear_extrap(*patch_data, i, i_intr, i_shft, depth);
                        break;
                    case 2:
                        (*patch_data)(i, depth) =
                            compute_quadratic_extrap(*patch_data, i, i_intr, i_shft, depth, codim);
                        break;
                    }
                }
            }
        }
    }
    return;
} // setPhysicalBoundaryConditions_cell

void
CartExtrapPhysBdryOp::setPhysicalBoundaryConditions_face(
    Patch& patch,
    const std::vector<std::pair<Box, std::pair<int, int> > >& bdry_fill_boxes)
{
    const Box& patch_box = patch.getBox();
    const hier::Index& patch_lower = patch_box.lower();
    const hier::Index& patch_upper = patch_box.upper();

    const int extrap_type =
        (d_extrap_type == "CONSTANT" ? 0 : (d_extrap_type == "LINEAR" ? 1 : (d_extrap_type == "QUADRATIC" ? 2 : -1)));
    if (extrap_type != 0 && extrap_type != 1 && extrap_type != 2)
    {
        TBOX_ERROR("CartExtrapPhysBdryOp::setPhysicalBoundaryConditions():\n"
                   << "  unknown extrapolation type: " << d_extrap_type << "\n"
                   << "  valid selections are: CONSTANT, LINEAR, or QUADRATIC" << std::endl);
    }

    // Set the physical boundary conditions for the specified patch data
    // indices.
    for (const auto& patch_data_idx : d_patch_data_indices)
    {
        VariableDatabase* var_db = VariableDatabase::getDatabase();
        std::shared_ptr<Variable > var;
        var_db->mapIndexToVariable(patch_data_idx, var);
        auto fc_var = std::static_pointer_cast<FaceVariable<double>>(var);
        if (!fc_var) continue;
        auto patch_data = std::static_pointer_cast<FaceData<double> >(patch.getPatchData(patch_data_idx));
        const Box& ghost_box = patch_data->getGhostBox();

        // Loop over the boundary fill boxes and extrapolate the data.
        for (const auto& bdry_fill_box_pair : bdry_fill_boxes)
        {
            const Box& bdry_fill_box = bdry_fill_box_pair.first;
            const unsigned int location_index = bdry_fill_box_pair.second.first;
            const int codim = bdry_fill_box_pair.second.second;
#if (NDIM == 2)
            const std::array<bool, NDIM> is_lower = { { PhysicalBoundaryUtilities::isLower(location_index, codim, 0),
                                                        PhysicalBoundaryUtilities::isLower(
                                                            location_index, codim, 1) } };
            const std::array<bool, NDIM> is_upper = { { PhysicalBoundaryUtilities::isUpper(location_index, codim, 0),
                                                        PhysicalBoundaryUtilities::isUpper(
                                                            location_index, codim, 1) } };
#endif
#if (NDIM == 3)
            const std::array<bool, NDIM> is_lower = { { PhysicalBoundaryUtilities::isLower(location_index, codim, 0),
                                                        PhysicalBoundaryUtilities::isLower(location_index, codim, 1),
                                                        PhysicalBoundaryUtilities::isLower(
                                                            location_index, codim, 2) } };
            const std::array<bool, NDIM> is_upper = { { PhysicalBoundaryUtilities::isUpper(location_index, codim, 0),
                                                        PhysicalBoundaryUtilities::isUpper(location_index, codim, 1),
                                                        PhysicalBoundaryUtilities::isUpper(
                                                            location_index, codim, 2) } };
#endif
            for (int depth = 0; depth < patch_data->getDepth(); ++depth)
            {
                for (unsigned int axis = 0; axis < NDIM; ++axis)
                {
                    Box intersect = bdry_fill_box * ghost_box;
                    FaceIterator end(FaceGeometry::end(intersect, axis));
                    for (FaceIterator b = FaceGeometry::begin(intersect, axis); b != end; b++)
                    {
                        const FaceIndex i = *b;
                        FaceIndex i_bdry = i;
                        IntVector i_shft = IntVector::getZero(Dimension(NDIM));
                        for (unsigned int d = 0; d < NDIM; ++d)
                        {
                            if (is_lower[d])
                            {
                                i_bdry((NDIM - axis + d) % NDIM) = patch_lower(d);
                                i_shft((NDIM - axis + d) % NDIM) = +1; // use interior data for extrapolation
                            }
                            else if (is_upper[d])
                            {
                                if (axis != d)
                                {
                                    i_bdry((NDIM - axis + d) % NDIM) = patch_upper(d);
                                }
                                else
                                {
                                    i_bdry((NDIM - axis + d) % NDIM) = patch_upper(d) + 1;
                                }
                                i_shft((NDIM - axis + d) % NDIM) = -1; // use interior data for extrapolation
                            }
                        }

                        // Perform constant, linear, or quadratic extrapolation.
                        switch (extrap_type)
                        {
                        case 0:
                            (*patch_data)(i, depth) = (*patch_data)(i_bdry, depth);
                            break;
                        case 1:
                            (*patch_data)(i, depth) = compute_linear_extrap(*patch_data, i, i_bdry, i_shft, depth);
                            break;
                        case 2:
                            (*patch_data)(i, depth) =
                                compute_quadratic_extrap(*patch_data, i, i_bdry, i_shft, depth, codim);
                            break;
                        }
                    }
                }
            }
        }
    }
    return;
} // setPhysicalBoundaryConditions_face

void
CartExtrapPhysBdryOp::setPhysicalBoundaryConditions_node(
    Patch& patch,
    const std::vector<std::pair<Box, std::pair<int, int> > >& bdry_fill_boxes)
{
    const Box& patch_box = patch.getBox();
    const hier::Index& patch_lower = patch_box.lower();
    const hier::Index& patch_upper = patch_box.upper();

    const int extrap_type =
        (d_extrap_type == "CONSTANT" ? 0 : (d_extrap_type == "LINEAR" ? 1 : (d_extrap_type == "QUADRATIC" ? 2 : -1)));
    if (extrap_type != 0 && extrap_type != 1 && extrap_type != 2)
    {
        TBOX_ERROR("CartExtrapPhysBdryOp::setPhysicalBoundaryConditions():\n"
                   << "  unknown extrapolation type: " << d_extrap_type << "\n"
                   << "  valid selections are: CONSTANT, LINEAR, or QUADRATIC" << std::endl);
    }

    // Set the physical boundary conditions for the specified patch data
    // indices.
    for (int patch_data_idx : d_patch_data_indices)
    {
        VariableDatabase* var_db = VariableDatabase::getDatabase();
        std::shared_ptr<Variable > var;
        var_db->mapIndexToVariable(patch_data_idx, var);
        auto nc_var = std::static_pointer_cast<NodeVariable<double>>(var);
        if (!nc_var) continue;
        auto patch_data = std::static_pointer_cast<NodeData<double> >(patch.getPatchData(patch_data_idx));
        const Box& ghost_box = patch_data->getGhostBox();

        // Loop over the boundary fill boxes and extrapolate the data.
        for (const auto& bdry_fill_box_pair : bdry_fill_boxes)
        {
            const Box& bdry_fill_box = bdry_fill_box_pair.first;
            const unsigned int location_index = bdry_fill_box_pair.second.first;
            const int codim = bdry_fill_box_pair.second.second;
#if (NDIM == 2)
            const std::array<bool, NDIM> is_lower = { { PhysicalBoundaryUtilities::isLower(location_index, codim, 0),
                                                        PhysicalBoundaryUtilities::isLower(
                                                            location_index, codim, 1) } };
            const std::array<bool, NDIM> is_upper = { { PhysicalBoundaryUtilities::isUpper(location_index, codim, 0),
                                                        PhysicalBoundaryUtilities::isUpper(
                                                            location_index, codim, 1) } };
#endif
#if (NDIM == 3)
            const std::array<bool, NDIM> is_lower = { { PhysicalBoundaryUtilities::isLower(location_index, codim, 0),
                                                        PhysicalBoundaryUtilities::isLower(location_index, codim, 1),
                                                        PhysicalBoundaryUtilities::isLower(
                                                            location_index, codim, 2) } };
            const std::array<bool, NDIM> is_upper = { { PhysicalBoundaryUtilities::isUpper(location_index, codim, 0),
                                                        PhysicalBoundaryUtilities::isUpper(location_index, codim, 1),
                                                        PhysicalBoundaryUtilities::isUpper(
                                                            location_index, codim, 2) } };
#endif
            // Loop over the boundary box indices and compute the
            // nearest interior index.
            for (int depth = 0; depth < patch_data->getDepth(); ++depth)
            {
                NodeIterator end(NodeGeometry::end(bdry_fill_box * ghost_box));
                for (NodeIterator b = NodeGeometry::begin(bdry_fill_box * ghost_box); b != end; b++)
                {
                    const NodeIndex& i = *b;
                    NodeIndex i_bdry = i;
                    IntVector i_shft = IntVector::getZero(Dimension(NDIM));
                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        if (is_lower[d])
                        {
                            i_bdry(d) = patch_lower(d);
                            i_shft(d) = +1; // use interior data for extrapolation
                        }
                        else if (is_upper[d])
                        {
                            i_bdry(d) = patch_upper(d) + 1;
                            i_shft(d) = -1; // use interior data for extrapolation
                        }
                    }

                    // Perform constant, linear, or quadratic extrapolation.
                    switch (extrap_type)
                    {
                    case 0:
                        (*patch_data)(i, depth) = (*patch_data)(i_bdry, depth);
                        break;
                    case 1:
                        (*patch_data)(i, depth) = compute_linear_extrap(*patch_data, i, i_bdry, i_shft, depth);
                        break;
                    case 2:
                        (*patch_data)(i, depth) =
                            compute_quadratic_extrap(*patch_data, i, i_bdry, i_shft, depth, codim);
                        break;
                    }
                }
            }
        }
    }
    return;
} // setPhysicalBoundaryConditions_node

void
CartExtrapPhysBdryOp::setPhysicalBoundaryConditions_side(
    Patch& patch,
    const std::vector<std::pair<Box, std::pair<int, int> > >& bdry_fill_boxes)
{
    const Box& patch_box = patch.getBox();
    const hier::Index& patch_lower = patch_box.lower();
    const hier::Index& patch_upper = patch_box.upper();

    const int extrap_type =
        (d_extrap_type == "CONSTANT" ? 0 : (d_extrap_type == "LINEAR" ? 1 : (d_extrap_type == "QUADRATIC" ? 2 : -1)));
    if (extrap_type != 0 && extrap_type != 1 && extrap_type != 2)
    {
        TBOX_ERROR("CartExtrapPhysBdryOp::setPhysicalBoundaryConditions():\n"
                   << "  unknown extrapolation type: " << d_extrap_type << "\n"
                   << "  valid selections are: CONSTANT, LINEAR, or QUADRATIC" << std::endl);
    }

    // Set the physical boundary conditions for the specified patch data
    // indices.
    for (const auto& patch_data_idx : d_patch_data_indices)
    {
        VariableDatabase* var_db = VariableDatabase::getDatabase();
        std::shared_ptr<Variable > var;
        var_db->mapIndexToVariable(patch_data_idx, var);
        auto sc_var = std::static_pointer_cast<SideVariable<double>>(var);
        if (!sc_var) continue;
        auto patch_data = std::static_pointer_cast<SideData<double> >(patch.getPatchData(patch_data_idx));
        const Box& ghost_box = patch_data->getGhostBox();

        // Loop over the boundary fill boxes and extrapolate the data.
        for (const auto& bdry_fill_box_map : bdry_fill_boxes)
        {
            const Box& bdry_fill_box = bdry_fill_box_map.first;
            const unsigned int location_index = bdry_fill_box_map.second.first;
            const int codim = bdry_fill_box_map.second.second;
#if (NDIM == 2)
            const std::array<bool, NDIM> is_lower = { { PhysicalBoundaryUtilities::isLower(location_index, codim, 0),
                                                        PhysicalBoundaryUtilities::isLower(
                                                            location_index, codim, 1) } };
            const std::array<bool, NDIM> is_upper = { { PhysicalBoundaryUtilities::isUpper(location_index, codim, 0),
                                                        PhysicalBoundaryUtilities::isUpper(
                                                            location_index, codim, 1) } };
#endif
#if (NDIM == 3)
            const std::array<bool, NDIM> is_lower = { { PhysicalBoundaryUtilities::isLower(location_index, codim, 0),
                                                        PhysicalBoundaryUtilities::isLower(location_index, codim, 1),
                                                        PhysicalBoundaryUtilities::isLower(
                                                            location_index, codim, 2) } };
            const std::array<bool, NDIM> is_upper = { { PhysicalBoundaryUtilities::isUpper(location_index, codim, 0),
                                                        PhysicalBoundaryUtilities::isUpper(location_index, codim, 1),
                                                        PhysicalBoundaryUtilities::isUpper(
                                                            location_index, codim, 2) } };
#endif
            for (int depth = 0; depth < patch_data->getDepth(); ++depth)
            {
                for (unsigned int axis = 0; axis < NDIM; ++axis)
                {
                    SideIterator end(SideGeometry::end(bdry_fill_box * ghost_box, axis));
                    for (SideIterator b = SideGeometry::begin(bdry_fill_box * ghost_box, axis); b != end; b++)
                    {
                        const SideIndex i = *b;
                        SideIndex i_bdry = i;
                        IntVector i_shft = IntVector::getZero(Dimension(NDIM));
                        for (unsigned int d = 0; d < NDIM; ++d)
                        {
                            if (is_lower[d])
                            {
                                i_bdry(d) = patch_lower(d);
                                i_shft(d) = +1; // use interior data for extrapolation
                            }
                            else if (is_upper[d])
                            {
                                if (axis != d)
                                {
                                    i_bdry(d) = patch_upper(d);
                                }
                                else
                                {
                                    i_bdry(d) = patch_upper(d) + 1;
                                }
                                i_shft(d) = -1; // use interior data for extrapolation
                            }
                        }

                        // Perform constant, linear, or quadratic extrapolation.
                        switch (extrap_type)
                        {
                        case 0:
                            (*patch_data)(i, depth) = (*patch_data)(i_bdry, depth);
                            break;
                        case 1:
                            (*patch_data)(i, depth) = compute_linear_extrap(*patch_data, i, i_bdry, i_shft, depth);
                            break;
                        case 2:
                            (*patch_data)(i, depth) =
                                compute_quadratic_extrap(*patch_data, i, i_bdry, i_shft, depth, codim);
                            break;
                        }
                    }
                }
            }
        }
    }
    return;
} // setPhysicalBoundaryConditions_side

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
