// Filename: CartFaceRobinPhysBdryOp.cpp
// Created on 21 May 2008 by Boyce Griffith
//
// Copyright (c) 2002-2017, Boyce Griffith
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//    * Redistributions of source code must retain the above copyright notice,
//      this list of conditions and the following disclaimer.
//
//    * Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions and the following disclaimer in the
//      documentation and/or other materials provided with the distribution.
//
//    * Neither the name of The University of North Carolina nor the names of
//      its contributors may be used to endorse or promote products derived from
//      this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ostream>
#include <set>
#include <vector>

#include "ArrayData.h"
#include "BoundaryBox.h"
#include "Box.h"
#include "CartesianPatchGeometry.h"
#include "ComponentSelector.h"
#include "IBTK_config.h"
#include "Index.h"
#include "IntVector.h"
#include "Patch.h"
#include "PatchGeometry.h"
#include "RobinBcCoefStrategy.h"
#include "SideData.h"
#include "SideGeometry.h"
#include "Variable.h"
#include "VariableDatabase.h"
#include "boost/array.hpp"
#include "ibtk/CartFaceRobinPhysBdryOp.h"
#include "ibtk/ExtendedRobinBcCoefStrategy.h"
#include "ibtk/PhysicalBoundaryUtilities.h"
#include "ibtk/RobinPhysBdryPatchStrategy.h"
#include "ibtk/namespaces.h" // IWYU pragma: keep
#include "tbox/Array.h"
#include "tbox/Pointer.h"
#include "tbox/Utilities.h"

// FORTRAN ROUTINES
#if (NDIM == 2)
#define CC_ROBIN_PHYS_BDRY_OP_1_X_FC IBTK_FC_FUNC(ccrobinphysbdryop1x2d, CCROBINPHYSBDRYOP1X2D)
#define CC_ROBIN_PHYS_BDRY_OP_1_Y_FC IBTK_FC_FUNC(ccrobinphysbdryop1y2d, CCROBINPHYSBDRYOP1Y2D)
#define FC_ROBIN_PHYS_BDRY_OP_1_X_FC IBTK_FC_FUNC(fcrobinphysbdryop1x2d, FCROBINPHYSBDRYOP1X2D)
#define FC_ROBIN_PHYS_BDRY_OP_1_Y_FC IBTK_FC_FUNC(fcrobinphysbdryop1y2d, FCROBINPHYSBDRYOP1Y2D)
#define FC_ROBIN_PHYS_BDRY_OP_2_FC IBTK_FC_FUNC(fcrobinphysbdryop22d, FCROBINPHYSBDRYOP22D)
#endif // if (NDIM == 2)

#if (NDIM == 3)
#define CC_ROBIN_PHYS_BDRY_OP_1_X_FC IBTK_FC_FUNC(ccrobinphysbdryop1x3d, CCROBINPHYSBDRYOP1X3D)
#define CC_ROBIN_PHYS_BDRY_OP_1_Y_FC IBTK_FC_FUNC(ccrobinphysbdryop1y3d, CCROBINPHYSBDRYOP1Y3D)
#define CC_ROBIN_PHYS_BDRY_OP_1_Z_FC IBTK_FC_FUNC(ccrobinphysbdryop1z3d, CCROBINPHYSBDRYOP1Z3D)
#define FC_ROBIN_PHYS_BDRY_OP_1_X_FC IBTK_FC_FUNC(fcrobinphysbdryop1x3d, FCROBINPHYSBDRYOP1X3D)
#define FC_ROBIN_PHYS_BDRY_OP_1_Y_FC IBTK_FC_FUNC(fcrobinphysbdryop1y3d, FCROBINPHYSBDRYOP1Y3D)
#define FC_ROBIN_PHYS_BDRY_OP_1_Z_FC IBTK_FC_FUNC(fcrobinphysbdryop1z3d, FCROBINPHYSBDRYOP1Z3D)
#define FC_ROBIN_PHYS_BDRY_OP_2_FC IBTK_FC_FUNC(fcrobinphysbdryop23d, FCROBINPHYSBDRYOP22D)
#define CC_ROBIN_PHYS_BDRY_OP_2_FC IBTK_FC_FUNC(ccrobinphysbdryop23d, CCROBINPHYSBDRYOP23D)
#define FC_ROBIN_PHYS_BDRY_OP_3_FC IBTK_FC_FUNC(fcrobinphysbdryop33d, FCROBINPHYSBDRYOP32D)
#endif // if (NDIM == 3)

extern "C"
{
    void CC_ROBIN_PHYS_BDRY_OP_1_X_FC(double* u,
                                      const int& u_gcw,
                                      const double* acoef,
                                      const double* bcoef,
                                      const double* gcoef,
                                      const int& location_index,
                                      const int& ilower0,
                                      const int& iupper0,
                                      const int& ilower1,
                                      const int& iupper1,
#if (NDIM == 3)
                                      const int& ilower2,
                                      const int& iupper2,
#endif
                                      const int& blower1,
                                      const int& bupper1,
#if (NDIM == 3)
                                      const int& blower2,
                                      const int& bupper2,
#endif
                                      const double* dx,
                                      const int& adjoint_op);

    void CC_ROBIN_PHYS_BDRY_OP_1_Y_FC(double* u,
                                      const int& u_gcw,
                                      const double* acoef,
                                      const double* bcoef,
                                      const double* gcoef,
                                      const int& location_index,
                                      const int& ilower0,
                                      const int& iupper0,
                                      const int& ilower1,
                                      const int& iupper1,
#if (NDIM == 3)
                                      const int& ilower2,
                                      const int& iupper2,
#endif
                                      const int& blower0,
                                      const int& bupper0,
#if (NDIM == 3)
                                      const int& blower2,
                                      const int& bupper2,
#endif
                                      const double* dx,
                                      const int& adjoint_op);

#if (NDIM == 3)
    void CC_ROBIN_PHYS_BDRY_OP_1_Z_FC(double* U,
                                      const int& U_gcw,
                                      const double* acoef,
                                      const double* bcoef,
                                      const double* gcoef,
                                      const int& location_index,
                                      const int& ilower0,
                                      const int& iupper0,
                                      const int& ilower1,
                                      const int& iupper1,
                                      const int& ilower2,
                                      const int& iupper2,
                                      const int& blower0,
                                      const int& bupper0,
                                      const int& blower1,
                                      const int& bupper1,
                                      const double* dx,
                                      const int& adjoint_op);
#endif

    void FC_ROBIN_PHYS_BDRY_OP_1_X_FC(double* u0,
                                      const int& u_gcw,
                                      const double* acoef,
                                      const double* bcoef,
                                      const double* gcoef,
                                      const int& location_index,
                                      const int& ilower0,
                                      const int& iupper0,
                                      const int& ilower1,
                                      const int& iupper1,
#if (NDIM == 3)
                                      const int& ilower2,
                                      const int& iupper2,
#endif
                                      const int& blower1,
                                      const int& bupper1,
#if (NDIM == 3)
                                      const int& blower2,
                                      const int& bupper2,
#endif
                                      const double* dx,
                                      const int& adjoint_op);

    void FC_ROBIN_PHYS_BDRY_OP_1_Y_FC(double* u1,
                                      const int& u_gcw,
                                      const double* acoef,
                                      const double* bcoef,
                                      const double* gcoef,
                                      const int& location_index,
                                      const int& ilower0,
                                      const int& iupper0,
                                      const int& ilower1,
                                      const int& iupper1,
#if (NDIM == 3)
                                      const int& ilower2,
                                      const int& iupper2,
#endif
                                      const int& blower0,
                                      const int& bupper0,
#if (NDIM == 3)
                                      const int& blower2,
                                      const int& bupper2,
#endif
                                      const double* dx,
                                      const int& adjoint_op);

#if (NDIM == 3)
    void FC_ROBIN_PHYS_BDRY_OP_1_Z_FC(double* u2,
                                      const int& u_gcw,
                                      const double* acoef,
                                      const double* bcoef,
                                      const double* gcoef,
                                      const int& location_index,
                                      const int& ilower0,
                                      const int& iupper0,
                                      const int& ilower1,
                                      const int& iupper1,
                                      const int& ilower2,
                                      const int& iupper2,
                                      const int& blower0,
                                      const int& bupper0,
                                      const int& blower1,
                                      const int& bupper1,
                                      const double* dx,
                                      const int& adjoint_op);
#endif

    void FC_ROBIN_PHYS_BDRY_OP_2_FC(double* u0,
                                    double* u1,
#if (NDIM == 3)
                                    double* u2,
#endif
                                    const int& u_gcw,
                                    const int& location_index,
                                    const int& ilower0,
                                    const int& iupper0,
                                    const int& ilower1,
                                    const int& iupper1,
#if (NDIM == 3)
                                    const int& ilower2,
                                    const int& iupper2,
#endif
                                    const int& blower0,
                                    const int& bupper0,
                                    const int& blower1,
                                    const int& bupper1,
#if (NDIM == 3)
                                    const int& blower2,
                                    const int& bupper2,
#endif
                                    const int& adjoint_op);

#if (NDIM == 3)
    void CC_ROBIN_PHYS_BDRY_OP_2_FC(double* U,
                                    const int& U_gcw,
                                    const int& location_index,
                                    const int& ilower0,
                                    const int& iupper0,
                                    const int& ilower1,
                                    const int& iupper1,
                                    const int& ilower2,
                                    const int& iupper2,
                                    const int& blower0,
                                    const int& bupper0,
                                    const int& blower1,
                                    const int& bupper1,
                                    const int& blower2,
                                    const int& bupper2,
                                    const int& adjoint_op);

    void FC_ROBIN_PHYS_BDRY_OP_3_FC(double* u0,
                                    double* u1,
                                    double* u2,
                                    const int& u_gcw,
                                    const int& location_index,
                                    const int& ilower0,
                                    const int& iupper0,
                                    const int& ilower1,
                                    const int& iupper1,
                                    const int& ilower2,
                                    const int& iupper2,
                                    const int& blower0,
                                    const int& bupper0,
                                    const int& blower1,
                                    const int& bupper1,
                                    const int& blower2,
                                    const int& bupper2,
                                    const int& adjoint_op);
#endif
}

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
static const int REFINE_OP_STENCIL_WIDTH = 1;

inline Box<NDIM>
compute_tangential_extension(const Box<NDIM>& box, const int data_axis)
{
    Box<NDIM> extended_box = box;
    extended_box.upper()(data_axis) += 1;
    return extended_box;
} // compute_tangential_extension
} // namespace

/////////////////////////////// PUBLIC ///////////////////////////////////////

CartFaceRobinPhysBdryOp::CartFaceRobinPhysBdryOp() : RobinPhysBdryPatchStrategy()
{
    // intentionally blank
    return;
} // CartFaceRobinPhysBdryOp

CartFaceRobinPhysBdryOp::CartFaceRobinPhysBdryOp(const int patch_data_index,
                                                 const std::vector<RobinBcCoefStrategy<NDIM>*>& bc_coefs,
                                                 const bool homogeneous_bc)
    : RobinPhysBdryPatchStrategy()
{
#if !defined(NDEBUG)
    TBOX_ASSERT(bc_coefs.size() == NDIM);
#endif
    setPatchDataIndex(patch_data_index);
    setPhysicalBcCoefs(bc_coefs);
    setHomogeneousBc(homogeneous_bc);
    return;
} // CartFaceRobinPhysBdryOp

CartFaceRobinPhysBdryOp::CartFaceRobinPhysBdryOp(const std::set<int>& patch_data_indices,
                                                 const std::vector<RobinBcCoefStrategy<NDIM>*>& bc_coefs,
                                                 const bool homogeneous_bc)
    : RobinPhysBdryPatchStrategy()
{
#if !defined(NDEBUG)
    TBOX_ASSERT(bc_coefs.size() == NDIM);
#endif
    setPatchDataIndices(patch_data_indices);
    setPhysicalBcCoefs(bc_coefs);
    setHomogeneousBc(homogeneous_bc);
    return;
} // CartFaceRobinPhysBdryOp

CartFaceRobinPhysBdryOp::CartFaceRobinPhysBdryOp(const ComponentSelector& patch_data_indices,
                                                 const std::vector<RobinBcCoefStrategy<NDIM>*>& bc_coefs,
                                                 const bool homogeneous_bc)
    : RobinPhysBdryPatchStrategy()
{
#if !defined(NDEBUG)
    TBOX_ASSERT(bc_coefs.size() == NDIM);
#endif
    setPatchDataIndices(patch_data_indices);
    setPhysicalBcCoefs(bc_coefs);
    setHomogeneousBc(homogeneous_bc);
    return;
} // CartFaceRobinPhysBdryOp

CartFaceRobinPhysBdryOp::~CartFaceRobinPhysBdryOp()
{
    // intentionally blank
    return;
} // ~CartFaceRobinPhysBdryOp

void
CartFaceRobinPhysBdryOp::setPhysicalBoundaryConditions(Patch<NDIM>& patch,
                                                       const double fill_time,
                                                       const IntVector<NDIM>& ghost_width_to_fill)
{
    if (ghost_width_to_fill == IntVector<NDIM>(0)) return;

    // Ensure the target patch data corresponds to a side centered variable and
    // that the proper number of boundary condition objects have been provided.
    for (std::set<int>::const_iterator cit = d_patch_data_indices.begin(); cit != d_patch_data_indices.end(); ++cit)
    {
        const int patch_data_idx = (*cit);
        Pointer<FaceData<NDIM, double> > patch_data = patch.getPatchData(patch_data_idx);
        if (!patch_data)
        {
            TBOX_ERROR("CartFaceRobinPhysBdryOp::setPhysicalBoundaryConditions():\n"
                       << "  patch data index " << patch_data_idx
                       << " does not correspond to a face-centered double precision variable." << std::endl);
        }
        if (NDIM * patch_data->getDepth() != static_cast<int>(d_bc_coefs.size()))
        {
            TBOX_ERROR("CartFaceRobinPhysBdryOp::setPhysicalBoundaryConditions():\n"
                       << "  data depth for patch data index " << patch_data_idx << " is " << patch_data->getDepth()
                       << "\n"
                       << "  but " << d_bc_coefs.size()
                       << " boundary condition coefficient objects were provided to the class "
                          "constructor."
                       << std::endl);
        }
    }

    // Set the boundary conditions along the co-dimension one boundary boxes,
    // then extrapolate those values to the co-dimension two and three boundary
    // boxes.
    static const bool adjoint_op = false;
    const Array<BoundaryBox<NDIM> > physical_codim1_boxes =
        PhysicalBoundaryUtilities::getPhysicalBoundaryCodim1Boxes(patch);
    for (std::set<int>::const_iterator cit = d_patch_data_indices.begin(); cit != d_patch_data_indices.end(); ++cit)
    {
        const int patch_data_idx = (*cit);
        fillGhostCellValuesCodim1Normal(
            patch_data_idx, physical_codim1_boxes, fill_time, ghost_width_to_fill, patch, adjoint_op);
    }
    for (std::set<int>::const_iterator cit = d_patch_data_indices.begin(); cit != d_patch_data_indices.end(); ++cit)
    {
        const int patch_data_idx = (*cit);
        fillGhostCellValuesCodim1Transverse(
            patch_data_idx, physical_codim1_boxes, fill_time, ghost_width_to_fill, patch, adjoint_op);
    }
    const Array<BoundaryBox<NDIM> > physical_codim2_boxes =
        PhysicalBoundaryUtilities::getPhysicalBoundaryCodim2Boxes(patch);
    for (std::set<int>::const_iterator cit = d_patch_data_indices.begin(); cit != d_patch_data_indices.end(); ++cit)
    {
        const int patch_data_idx = (*cit);
        fillGhostCellValuesCodim2(patch_data_idx, physical_codim2_boxes, ghost_width_to_fill, patch, adjoint_op);
    }
#if (NDIM > 2)
    const Array<BoundaryBox<NDIM> > physical_codim3_boxes =
        PhysicalBoundaryUtilities::getPhysicalBoundaryCodim3Boxes(patch);
    for (std::set<int>::const_iterator cit = d_patch_data_indices.begin(); cit != d_patch_data_indices.end(); ++cit)
    {
        const int patch_data_idx = (*cit);
        fillGhostCellValuesCodim3(patch_data_idx, physical_codim3_boxes, ghost_width_to_fill, patch, adjoint_op);
    }
#endif
    return;
} // setPhysicalBoundaryConditions

IntVector<NDIM>
CartFaceRobinPhysBdryOp::getRefineOpStencilWidth() const
{
    return REFINE_OP_STENCIL_WIDTH;
} // getRefineOpStencilWidth

void
CartFaceRobinPhysBdryOp::accumulateFromPhysicalBoundaryData(Patch<NDIM>& patch,
                                                            const double fill_time,
                                                            const IntVector<NDIM>& ghost_width_to_fill)
{
    if (ghost_width_to_fill == IntVector<NDIM>(0)) return;

    // Ensure the target patch data corresponds to a side centered variable and
    // that the proper number of boundary condition objects have been provided.
    for (std::set<int>::const_iterator cit = d_patch_data_indices.begin(); cit != d_patch_data_indices.end(); ++cit)
    {
        const int patch_data_idx = (*cit);
        Pointer<FaceData<NDIM, double> > patch_data = patch.getPatchData(patch_data_idx);
        if (!patch_data)
        {
            TBOX_ERROR("CartFaceRobinPhysBdryOp::accumulateFromPhysicalBoundaryData():\n"
                       << "  patch data index " << patch_data_idx
                       << " does not correspond to a face-centered double precision variable." << std::endl);
        }
        if (NDIM * patch_data->getDepth() != static_cast<int>(d_bc_coefs.size()))
        {
            TBOX_ERROR("CartFaceRobinPhysBdryOp::accumulateFromPhysicalBoundaryData():\n"
                       << "  data depth for patch data index " << patch_data_idx << " is " << patch_data->getDepth()
                       << "\n"
                       << "  but " << d_bc_coefs.size()
                       << " boundary condition coefficient objects were provided to the class "
                          "constructor."
                       << std::endl);
        }
    }

    // Set the boundary conditions along the co-dimension one boundary boxes,
    // then extrapolate those values to the co-dimension two and three boundary
    // boxes.
    static const bool adjoint_op = true;
#if (NDIM > 2)
    const Array<BoundaryBox<NDIM> > physical_codim3_boxes =
        PhysicalBoundaryUtilities::getPhysicalBoundaryCodim3Boxes(patch);
    for (std::set<int>::const_iterator cit = d_patch_data_indices.begin(); cit != d_patch_data_indices.end(); ++cit)
    {
        const int patch_data_idx = (*cit);
        fillGhostCellValuesCodim3(patch_data_idx, physical_codim3_boxes, ghost_width_to_fill, patch, adjoint_op);
    }
#endif
    const Array<BoundaryBox<NDIM> > physical_codim2_boxes =
        PhysicalBoundaryUtilities::getPhysicalBoundaryCodim2Boxes(patch);
    for (std::set<int>::const_iterator cit = d_patch_data_indices.begin(); cit != d_patch_data_indices.end(); ++cit)
    {
        const int patch_data_idx = (*cit);
        fillGhostCellValuesCodim2(patch_data_idx, physical_codim2_boxes, ghost_width_to_fill, patch, adjoint_op);
    }
    const Array<BoundaryBox<NDIM> > physical_codim1_boxes =
        PhysicalBoundaryUtilities::getPhysicalBoundaryCodim1Boxes(patch);
    for (std::set<int>::const_iterator cit = d_patch_data_indices.begin(); cit != d_patch_data_indices.end(); ++cit)
    {
        const int patch_data_idx = (*cit);
        fillGhostCellValuesCodim1Transverse(
            patch_data_idx, physical_codim1_boxes, fill_time, ghost_width_to_fill, patch, adjoint_op);
    }
    for (std::set<int>::const_iterator cit = d_patch_data_indices.begin(); cit != d_patch_data_indices.end(); ++cit)
    {
        const int patch_data_idx = (*cit);
        fillGhostCellValuesCodim1Normal(
            patch_data_idx, physical_codim1_boxes, fill_time, ghost_width_to_fill, patch, adjoint_op);
    }
    return;
} // accumulateFromPhysicalBoundaryData

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

void
CartFaceRobinPhysBdryOp::fillGhostCellValuesCodim1Normal(const int patch_data_idx,
                                                         const Array<BoundaryBox<NDIM> >& physical_codim1_boxes,
                                                         const double fill_time,
                                                         const IntVector<NDIM>& ghost_width_to_fill,
                                                         Patch<NDIM>& patch,
                                                         const bool adjoint_op)
{
    const int n_physical_codim1_boxes = physical_codim1_boxes.size();
    if (n_physical_codim1_boxes == 0) return;

    const Box<NDIM>& patch_box = patch.getBox();
    Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch.getPatchGeometry();
    const double* const dx = pgeom->getDx();
    Pointer<FaceData<NDIM, double> > patch_data = patch.getPatchData(patch_data_idx);
    const int patch_data_depth = patch_data->getDepth();
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<Variable<NDIM> > var;
    var_db->mapIndexToVariable(patch_data_idx, var);
    const int patch_data_gcw = (patch_data->getGhostCellWidth()).max();
#if !defined(NDEBUG)
    if (patch_data_gcw != (patch_data->getGhostCellWidth()).min())
    {
        TBOX_ERROR(
            "CartFaceRobinPhysBdryOp::fillGhostCellValuesCodim1Normal():\n"
            "  patch data for patch data index "
            << patch_data_idx << " does not have uniform ghost cell widths." << std::endl);
    }
#endif
    const IntVector<NDIM> gcw_to_fill = IntVector<NDIM>::min(patch_data->getGhostCellWidth(), ghost_width_to_fill);

    // Set the boundary condition coefficients and then set the ghost cell
    // values.
    for (int n = 0; n < n_physical_codim1_boxes; ++n)
    {
        const BoundaryBox<NDIM>& bdry_box = physical_codim1_boxes[n];
        const unsigned int location_index = bdry_box.getLocationIndex();
        const unsigned int bdry_normal_axis = location_index / 2;
        const Box<NDIM> bc_fill_box = pgeom->getBoundaryFillBox(bdry_box, patch_box, gcw_to_fill);
        const BoundaryBox<NDIM> trimmed_bdry_box(
            bdry_box.getBox() * bc_fill_box, bdry_box.getBoundaryType(), location_index);
        const Box<NDIM> bc_coef_box = PhysicalBoundaryUtilities::makeSideBoundaryCodim1Box(trimmed_bdry_box);
        Pointer<ArrayData<NDIM, double> > acoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);
        Pointer<ArrayData<NDIM, double> > bcoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);
        Pointer<ArrayData<NDIM, double> > gcoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);
        for (int d = 0; d < patch_data_depth; ++d)
        {
            RobinBcCoefStrategy<NDIM>* bc_coef = d_bc_coefs[NDIM * d + bdry_normal_axis];
            ExtendedRobinBcCoefStrategy* const extended_bc_coef = dynamic_cast<ExtendedRobinBcCoefStrategy*>(bc_coef);
            if (extended_bc_coef)
            {
                extended_bc_coef->setTargetPatchDataIndex(patch_data_idx);
                extended_bc_coef->setHomogeneousBc(d_homogeneous_bc);
            }
            bc_coef->setBcCoefs(acoef_data, bcoef_data, gcoef_data, var, patch, trimmed_bdry_box, fill_time);
            if (d_homogeneous_bc && !extended_bc_coef) gcoef_data->fillAll(0.0);
            if (extended_bc_coef) extended_bc_coef->clearTargetPatchDataIndex();
            if (location_index == 0 || location_index == 1)
            {
                FC_ROBIN_PHYS_BDRY_OP_1_X_FC(patch_data->getPointer(bdry_normal_axis, d),
                                             patch_data_gcw,
                                             acoef_data->getPointer(),
                                             bcoef_data->getPointer(),
                                             gcoef_data->getPointer(),
                                             location_index,
                                             patch_box.lower(0),
                                             patch_box.upper(0),
                                             patch_box.lower(1),
                                             patch_box.upper(1),
#if (NDIM == 3)
                                             patch_box.lower(2),
                                             patch_box.upper(2),
#endif
                                             bc_coef_box.lower(1),
                                             bc_coef_box.upper(1),
#if (NDIM == 3)
                                             bc_coef_box.lower(2),
                                             bc_coef_box.upper(2),
#endif
                                             dx,
                                             adjoint_op ? 1 : 0);
            }
            else if (location_index == 2 || location_index == 3)
            {
                SC_ROBIN_PHYS_BDRY_OP_1_Y_FC(patch_data->getPointer(bdry_normal_axis, d),
                                             patch_data_gcw,
                                             acoef_data->getPointer(),
                                             bcoef_data->getPointer(),
                                             gcoef_data->getPointer(),
                                             location_index,
                                             patch_box.lower(0),
                                             patch_box.upper(0),
                                             patch_box.lower(1),
                                             patch_box.upper(1),
#if (NDIM == 3)
                                             patch_box.lower(2),
                                             patch_box.upper(2),
#endif
                                             bc_coef_box.lower(0),
                                             bc_coef_box.upper(0),
#if (NDIM == 3)
                                             bc_coef_box.lower(2),
                                             bc_coef_box.upper(2),
#endif
                                             dx,
                                             adjoint_op ? 1 : 0);
            }
#if (NDIM == 3)
            else if (location_index == 4 || location_index == 5)
            {
                SC_ROBIN_PHYS_BDRY_OP_1_Z_FC(patch_data->getPointer(bdry_normal_axis, d),
                                             patch_data_gcw,
                                             acoef_data->getPointer(),
                                             bcoef_data->getPointer(),
                                             gcoef_data->getPointer(),
                                             location_index,
                                             patch_box.lower(0),
                                             patch_box.upper(0),
                                             patch_box.lower(1),
                                             patch_box.upper(1),
                                             patch_box.lower(2),
                                             patch_box.upper(2),
                                             bc_coef_box.lower(0),
                                             bc_coef_box.upper(0),
                                             bc_coef_box.lower(1),
                                             bc_coef_box.upper(1),
                                             dx,
                                             adjoint_op ? 1 : 0);
            }
#endif
        }
    }
    return;
} // fillGhostCellValuesCodim1Normal

void
CartFaceRobinPhysBdryOp::fillGhostCellValuesCodim1Transverse(const int patch_data_idx,
                                                             const Array<BoundaryBox<NDIM> >& physical_codim1_boxes,
                                                             const double fill_time,
                                                             const IntVector<NDIM>& ghost_width_to_fill,
                                                             Patch<NDIM>& patch,
                                                             const bool adjoint_op)
{
    const int n_physical_codim1_boxes = physical_codim1_boxes.size();
    if (n_physical_codim1_boxes == 0) return;

    const Box<NDIM>& patch_box = patch.getBox();
    Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch.getPatchGeometry();
    const double* const dx = pgeom->getDx();
    const double* const patch_x_lower = pgeom->getXLower();
    const double* const patch_x_upper = pgeom->getXUpper();
    Pointer<SideData<NDIM, double> > patch_data = patch.getPatchData(patch_data_idx);
    const int patch_data_depth = patch_data->getDepth();
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<Variable<NDIM> > var;
    var_db->mapIndexToVariable(patch_data_idx, var);
    const int patch_data_gcw = (patch_data->getGhostCellWidth()).max();
#if !defined(NDEBUG)
    if (patch_data_gcw != (patch_data->getGhostCellWidth()).min())
    {
        TBOX_ERROR(
            "CartFaceRobinPhysBdryOp::fillGhostCellValuesCodim1Transverse():\n"
            "  patch data for patch data index "
            << patch_data_idx << " does not have uniform ghost cell widths." << std::endl);
    }
#endif
    const IntVector<NDIM> gcw_to_fill = IntVector<NDIM>::min(patch_data->getGhostCellWidth(), ghost_width_to_fill);

    boost::array<Box<NDIM>, NDIM> side_box;
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        side_box[axis] = SideGeometry<NDIM>::toSideBox(patch_box, axis);
    }

    const IntVector<NDIM>& ratio_to_level_zero = pgeom->getRatio();
    Array<Array<bool> > touches_regular_bdry(NDIM), touches_periodic_bdry(NDIM);
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        touches_regular_bdry[axis].resizeArray(2);
        touches_periodic_bdry[axis].resizeArray(2);
        for (int upperlower = 0; upperlower < 2; ++upperlower)
        {
            touches_regular_bdry[axis][upperlower] = pgeom->getTouchesRegularBoundary(axis, upperlower);
            touches_periodic_bdry[axis][upperlower] = pgeom->getTouchesPeriodicBoundary(axis, upperlower);
        }
    }

    // Set the boundary condition coefficients and then set the ghost cell
    // values.
    for (int n = 0; n < n_physical_codim1_boxes; ++n)
    {
        const BoundaryBox<NDIM>& bdry_box = physical_codim1_boxes[n];
        const unsigned int location_index = bdry_box.getLocationIndex();
        const unsigned int bdry_normal_axis = location_index / 2;
        const Box<NDIM> bc_fill_box = pgeom->getBoundaryFillBox(bdry_box, patch_box, gcw_to_fill);
        const BoundaryBox<NDIM> trimmed_bdry_box(
            bdry_box.getBox() * bc_fill_box, bdry_box.getBoundaryType(), location_index);
        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            if (axis != bdry_normal_axis)
            {
                const Box<NDIM> bc_coef_box = compute_tangential_extension(
                    PhysicalBoundaryUtilities::makeSideBoundaryCodim1Box(trimmed_bdry_box), axis);
                Pointer<ArrayData<NDIM, double> > acoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);
                Pointer<ArrayData<NDIM, double> > bcoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);
                Pointer<ArrayData<NDIM, double> > gcoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);

                // Temporarily reset the patch geometry object associated with
                // the patch so that boundary conditions are set at the correct
                // spatial locations.
                boost::array<double, NDIM> shifted_patch_x_lower, shifted_patch_x_upper;
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    shifted_patch_x_lower[d] = patch_x_lower[d];
                    shifted_patch_x_upper[d] = patch_x_upper[d];
                }
                shifted_patch_x_lower[axis] -= 0.5 * dx[axis];
                shifted_patch_x_upper[axis] -= 0.5 * dx[axis];
                patch.setPatchGeometry(new CartesianPatchGeometry<NDIM>(ratio_to_level_zero,
                                                                        touches_regular_bdry,
                                                                        touches_periodic_bdry,
                                                                        dx,
                                                                        shifted_patch_x_lower.data(),
                                                                        shifted_patch_x_upper.data()));

                // Set the boundary condition coefficients.
                for (int d = 0; d < patch_data_depth; ++d)
                {
                    RobinBcCoefStrategy<NDIM>* bc_coef = d_bc_coefs[NDIM * d + axis];
                    ExtendedRobinBcCoefStrategy* const extended_bc_coef =
                        dynamic_cast<ExtendedRobinBcCoefStrategy*>(bc_coef);
                    if (extended_bc_coef)
                    {
                        extended_bc_coef->setTargetPatchDataIndex(patch_data_idx);
                        extended_bc_coef->setHomogeneousBc(d_homogeneous_bc);
                    }
                    bc_coef->setBcCoefs(acoef_data, bcoef_data, gcoef_data, var, patch, trimmed_bdry_box, fill_time);
                    if (d_homogeneous_bc && !extended_bc_coef) gcoef_data->fillAll(0.0);
                    if (extended_bc_coef) extended_bc_coef->clearTargetPatchDataIndex();

                    // Restore the original patch geometry object.
                    patch.setPatchGeometry(pgeom);

                    // Set the boundary values.
                    if (location_index == 0 || location_index == 1)
                    {
                        CC_ROBIN_PHYS_BDRY_OP_1_X_FC(patch_data->getPointer(axis, d),
                                                     patch_data_gcw,
                                                     acoef_data->getPointer(),
                                                     bcoef_data->getPointer(),
                                                     gcoef_data->getPointer(),
                                                     location_index,
                                                     side_box[axis].lower(0),
                                                     side_box[axis].upper(0),
                                                     side_box[axis].lower(1),
                                                     side_box[axis].upper(1),
#if (NDIM == 3)
                                                     side_box[axis].lower(2),
                                                     side_box[axis].upper(2),
#endif
                                                     bc_coef_box.lower(1),
                                                     bc_coef_box.upper(1),
#if (NDIM == 3)
                                                     bc_coef_box.lower(2),
                                                     bc_coef_box.upper(2),
#endif
                                                     dx,
                                                     adjoint_op ? 1 : 0);
                    }
                    else if (location_index == 2 || location_index == 3)
                    {
                        CC_ROBIN_PHYS_BDRY_OP_1_Y_FC(patch_data->getPointer(axis, d),
                                                     patch_data_gcw,
                                                     acoef_data->getPointer(),
                                                     bcoef_data->getPointer(),
                                                     gcoef_data->getPointer(),
                                                     location_index,
                                                     side_box[axis].lower(0),
                                                     side_box[axis].upper(0),
                                                     side_box[axis].lower(1),
                                                     side_box[axis].upper(1),
#if (NDIM == 3)
                                                     side_box[axis].lower(2),
                                                     side_box[axis].upper(2),
#endif
                                                     bc_coef_box.lower(0),
                                                     bc_coef_box.upper(0),
#if (NDIM == 3)
                                                     bc_coef_box.lower(2),
                                                     bc_coef_box.upper(2),
#endif
                                                     dx,
                                                     adjoint_op ? 1 : 0);
                    }
#if (NDIM == 3)
                    else if (location_index == 4 || location_index == 5)
                    {
                        CC_ROBIN_PHYS_BDRY_OP_1_Z_FC(patch_data->getPointer(axis, d),
                                                     patch_data_gcw,
                                                     acoef_data->getPointer(),
                                                     bcoef_data->getPointer(),
                                                     gcoef_data->getPointer(),
                                                     location_index,
                                                     side_box[axis].lower(0),
                                                     side_box[axis].upper(0),
                                                     side_box[axis].lower(1),
                                                     side_box[axis].upper(1),
                                                     side_box[axis].lower(2),
                                                     side_box[axis].upper(2),
                                                     bc_coef_box.lower(0),
                                                     bc_coef_box.upper(0),
                                                     bc_coef_box.lower(1),
                                                     bc_coef_box.upper(1),
                                                     dx,
                                                     adjoint_op ? 1 : 0);
                    }
#endif
                }
            }
        }
    }
    return;
} // fillGhostCellValuesCodim1Transverse

void
CartFaceRobinPhysBdryOp::fillGhostCellValuesCodim2(const int patch_data_idx,
                                                   const Array<BoundaryBox<NDIM> >& physical_codim2_boxes,
                                                   const IntVector<NDIM>& ghost_width_to_fill,
                                                   const Patch<NDIM>& patch,
                                                   const bool adjoint_op)
{
    const int n_physical_codim2_boxes = physical_codim2_boxes.size();
    if (n_physical_codim2_boxes == 0) return;

    const Box<NDIM>& patch_box = patch.getBox();
    Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch.getPatchGeometry();
    Pointer<SideData<NDIM, double> > patch_data = patch.getPatchData(patch_data_idx);
    const int patch_data_depth = patch_data->getDepth();
    const int patch_data_gcw = (patch_data->getGhostCellWidth()).max();
#if !defined(NDEBUG)
    if (patch_data_gcw != (patch_data->getGhostCellWidth()).min())
    {
        TBOX_ERROR(
            "CartFaceRobinPhysBdryOp::fillGhostCellValuesCodim2():\n"
            "  patch data for patch data index "
            << patch_data_idx << " does not have uniform ghost cell widths." << std::endl);
    }
#endif
    const IntVector<NDIM> gcw_to_fill = IntVector<NDIM>::min(patch_data->getGhostCellWidth(), ghost_width_to_fill);

#if (NDIM == 3)
    boost::array<Box<NDIM>, NDIM> side_box;
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        side_box[axis] = SideGeometry<NDIM>::toSideBox(patch_box, axis);
    }
#endif

    for (int n = 0; n < n_physical_codim2_boxes; ++n)
    {
        const BoundaryBox<NDIM>& bdry_box = physical_codim2_boxes[n];
        const unsigned int location_index = bdry_box.getLocationIndex();
        const Box<NDIM> bc_fill_box = pgeom->getBoundaryFillBox(bdry_box, patch_box, gcw_to_fill);
        for (int d = 0; d < patch_data_depth; ++d)
        {
            SC_ROBIN_PHYS_BDRY_OP_2_FC(patch_data->getPointer(0, d),
                                       patch_data->getPointer(1, d),
#if (NDIM == 3)
                                       patch_data->getPointer(2, d),
#endif
                                       patch_data_gcw,
                                       location_index,
                                       patch_box.lower(0),
                                       patch_box.upper(0),
                                       patch_box.lower(1),
                                       patch_box.upper(1),
#if (NDIM == 3)
                                       patch_box.lower(2),
                                       patch_box.upper(2),
#endif
                                       bc_fill_box.lower(0),
                                       bc_fill_box.upper(0),
                                       bc_fill_box.lower(1),
                                       bc_fill_box.upper(1),
#if (NDIM == 3)
                                       bc_fill_box.lower(2),
                                       bc_fill_box.upper(2),
#endif
                                       adjoint_op ? 1 : 0);
#if (NDIM == 3)
            if (location_index < 4)
            {
                static const unsigned int axis = 0;
                Box<NDIM> side_bc_fill_box = SideGeometry<NDIM>::toSideBox(bc_fill_box, axis);
                CC_ROBIN_PHYS_BDRY_OP_2_FC(patch_data->getPointer(axis, d),
                                           patch_data_gcw,
                                           location_index,
                                           side_box[axis].lower(0),
                                           side_box[axis].upper(0),
                                           side_box[axis].lower(1),
                                           side_box[axis].upper(1),
                                           side_box[axis].lower(2),
                                           side_box[axis].upper(2),
                                           side_bc_fill_box.lower(0),
                                           side_bc_fill_box.upper(0),
                                           side_bc_fill_box.lower(1),
                                           side_bc_fill_box.upper(1),
                                           side_bc_fill_box.lower(2),
                                           side_bc_fill_box.upper(2),
                                           adjoint_op ? 1 : 0);
            }
            else if (location_index >= 4 && location_index < 8)
            {
                static const unsigned int axis = 1;
                Box<NDIM> side_bc_fill_box = SideGeometry<NDIM>::toSideBox(bc_fill_box, axis);
                CC_ROBIN_PHYS_BDRY_OP_2_FC(patch_data->getPointer(axis, d),
                                           patch_data_gcw,
                                           location_index,
                                           side_box[axis].lower(0),
                                           side_box[axis].upper(0),
                                           side_box[axis].lower(1),
                                           side_box[axis].upper(1),
                                           side_box[axis].lower(2),
                                           side_box[axis].upper(2),
                                           side_bc_fill_box.lower(0),
                                           side_bc_fill_box.upper(0),
                                           side_bc_fill_box.lower(1),
                                           side_bc_fill_box.upper(1),
                                           side_bc_fill_box.lower(2),
                                           side_bc_fill_box.upper(2),
                                           adjoint_op ? 1 : 0);
            }
            else if (location_index >= 8 && location_index < 12)
            {
                static const unsigned int axis = 2;
                Box<NDIM> side_bc_fill_box = SideGeometry<NDIM>::toSideBox(bc_fill_box, axis);
                CC_ROBIN_PHYS_BDRY_OP_2_FC(patch_data->getPointer(axis, d),
                                           patch_data_gcw,
                                           location_index,
                                           side_box[axis].lower(0),
                                           side_box[axis].upper(0),
                                           side_box[axis].lower(1),
                                           side_box[axis].upper(1),
                                           side_box[axis].lower(2),
                                           side_box[axis].upper(2),
                                           side_bc_fill_box.lower(0),
                                           side_bc_fill_box.upper(0),
                                           side_bc_fill_box.lower(1),
                                           side_bc_fill_box.upper(1),
                                           side_bc_fill_box.lower(2),
                                           side_bc_fill_box.upper(2),
                                           adjoint_op ? 1 : 0);
            }
#endif
        }
    }
    return;
} // fillGhostCellValuesCodim2

#if (NDIM > 2)
void
CartFaceRobinPhysBdryOp::fillGhostCellValuesCodim3(const int patch_data_idx,
                                                   const Array<BoundaryBox<NDIM> >& physical_codim3_boxes,
                                                   const IntVector<NDIM>& ghost_width_to_fill,
                                                   const Patch<NDIM>& patch,
                                                   const bool adjoint_op)
{
    const int n_physical_codim3_boxes = physical_codim3_boxes.size();
    if (n_physical_codim3_boxes == 0) return;

    const Box<NDIM>& patch_box = patch.getBox();
    Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch.getPatchGeometry();
    Pointer<SideData<NDIM, double> > patch_data = patch.getPatchData(patch_data_idx);
    const int patch_data_depth = patch_data->getDepth();
    const int patch_data_gcw = (patch_data->getGhostCellWidth()).max();
#if !defined(NDEBUG)
    if (patch_data_gcw != (patch_data->getGhostCellWidth()).min())
    {
        TBOX_ERROR(
            "CartFaceRobinPhysBdryOp::fillGhostCellValuesCodim3():\n"
            "  patch data for patch data index "
            << patch_data_idx << " does not have uniform ghost cell widths." << std::endl);
    }
#endif
    const IntVector<NDIM> gcw_to_fill = IntVector<NDIM>::min(patch_data->getGhostCellWidth(), ghost_width_to_fill);

    for (int n = 0; n < n_physical_codim3_boxes; ++n)
    {
        const BoundaryBox<NDIM>& bdry_box = physical_codim3_boxes[n];
        const unsigned int location_index = bdry_box.getLocationIndex();
        const Box<NDIM> bc_fill_box = pgeom->getBoundaryFillBox(bdry_box, patch_box, gcw_to_fill);
        for (int d = 0; d < patch_data_depth; ++d)
        {
            SC_ROBIN_PHYS_BDRY_OP_3_FC(patch_data->getPointer(0, d),
                                       patch_data->getPointer(1, d),
                                       patch_data->getPointer(2, d),
                                       patch_data_gcw,
                                       location_index,
                                       patch_box.lower(0),
                                       patch_box.upper(0),
                                       patch_box.lower(1),
                                       patch_box.upper(1),
                                       patch_box.lower(2),
                                       patch_box.upper(2),
                                       bc_fill_box.lower(0),
                                       bc_fill_box.upper(0),
                                       bc_fill_box.lower(1),
                                       bc_fill_box.upper(1),
                                       bc_fill_box.lower(2),
                                       bc_fill_box.upper(2),
                                       adjoint_op ? 1 : 0);
        }
    }
    return;
} // fillGhostCellValuesCodim3
#endif

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
