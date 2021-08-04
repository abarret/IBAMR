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

#include "ibtk/CartCellDoubleQuadraticRefine.h"
#include "ibtk/ibtk_utilities.h"
#include "ibtk/namespaces.h" // IWYU pragma: keep

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/hier/Patch.h"

#include <array>
#include <string>
#include <vector>

namespace SAMRAI
{
namespace hier
{

class Variable;
} // namespace hier
} // namespace SAMRAI

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

const std::string CartCellDoubleQuadraticRefine::s_op_name = "QUADRATIC_REFINE";

namespace
{
static const int REFINE_OP_PRIORITY = 0;
static const int REFINE_OP_STENCIL_WIDTH = 1;

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
} // namespace

/////////////////////////////// PUBLIC ///////////////////////////////////////

bool
CartCellDoubleQuadraticRefine::findRefineOperator(const std::shared_ptr<Variable >& var, const std::string& op_name) const
{
    const std::shared_ptr<CellVariable<double> > cc_var = var;
    return (cc_var && op_name == s_op_name);
} // findRefineOperator

const std::string&
CartCellDoubleQuadraticRefine::getOperatorName() const
{
    return s_op_name;
} // getOperatorName

int
CartCellDoubleQuadraticRefine::getOperatorPriority() const
{
    return REFINE_OP_PRIORITY;
} // getOperatorPriority

IntVector
CartCellDoubleQuadraticRefine::getStencilWidth() const
{
    return REFINE_OP_STENCIL_WIDTH;
} // getStencilWidth

void
CartCellDoubleQuadraticRefine::refine(Patch& fine,
                                      const Patch& coarse,
                                      const int dst_component,
                                      const int src_component,
                                      const Box& fine_box,
                                      const IntVector& ratio) const
{
    // Get the patch data.
    std::shared_ptr<CellData<double> > fdata = std::static_pointer_cast<CellData<double> >(fine.getPatchData(dst_component));
    std::shared_ptr<CellData<double> > cdata = std::static_pointer_cast<CellData<double> >(coarse.getPatchData(src_component));
#if !defined(NDEBUG)
    TBOX_ASSERT(fdata);
    TBOX_ASSERT(cdata);
    TBOX_ASSERT(fdata->getDepth() == cdata->getDepth());
#endif
    const int data_depth = fdata->getDepth();

    const Box& patch_box_fine = fine.getBox();
    const hier::Index& patch_lower_fine = patch_box_fine.lower();
    std::shared_ptr<CartesianPatchGeometry > pgeom_fine = std::static_pointer_cast<CartesianPatchGeometry >(fine.getPatchGeometry());
    const double* const XLower_fine = pgeom_fine->getXLower();
    const double* const dx_fine = pgeom_fine->getDx();

    const Box& patch_box_crse = coarse.getBox();
    const hier::Index& patch_lower_crse = patch_box_crse.lower();
    std::shared_ptr<CartesianPatchGeometry > pgeom_crse = std::static_pointer_cast<CartesianPatchGeometry >(coarse.getPatchGeometry());
    const double* const XLower_crse = pgeom_crse->getXLower();
    const double* const dx_crse = pgeom_crse->getDx();

    // Set all values in the fine box via quadratic interpolation from the
    // overlying coarse grid data.
    for (Box::Iterator b(fine_box); b; b++)
    {
        const hier::Index& i_fine = b();
        const hier::Index i_crse = coarsen(i_fine, ratio);

        // Determine the interpolation stencil in the coarse index space.
        Box stencil_box_crse(i_crse, i_crse);
        stencil_box_crse.grow(IntVector(1));

        // Determine the interpolation weights.
        static const int degree = 2;
        std::array<std::array<double, degree + 1>, NDIM> wgts(array_constant<std::array<double, degree + 1>, NDIM>(
            std::array<double, degree + 1>(array_constant<double, degree + 1>(0.0))));
        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            const double X =
                XLower_fine[axis] + dx_fine[axis] * (static_cast<double>(i_fine(axis) - patch_lower_fine(axis)) + 0.5);
            std::vector<double> X_crse(degree + 1, 0.0);
            for (int i_crse = stencil_box_crse.lower()(axis), k = 0; i_crse <= stencil_box_crse.upper()(axis);
                 ++i_crse, ++k)
            {
                X_crse[k] =
                    XLower_crse[axis] + dx_crse[axis] * (static_cast<double>(i_crse - patch_lower_crse(axis)) + 0.5);
            }
            wgts[axis][0] = ((X - X_crse[1]) * (X - X_crse[2])) / ((X_crse[0] - X_crse[1]) * (X_crse[0] - X_crse[2]));
            wgts[axis][1] = ((X - X_crse[0]) * (X - X_crse[2])) / ((X_crse[1] - X_crse[0]) * (X_crse[1] - X_crse[2]));
            wgts[axis][2] = ((X - X_crse[0]) * (X - X_crse[1])) / ((X_crse[2] - X_crse[0]) * (X_crse[2] - X_crse[1]));
        }

        // Interpolate from the coarse grid to the fine grid.
        hier::Index i_intrp(Dimension(NDIM));
        for (int d = 0; d < data_depth; ++d)
        {
            (*fdata)(i_fine, d) = 0.0;
#if (NDIM > 2)
            for (int i2 = 0; i2 <= degree; ++i2)
            {
                const double& wgt2 = wgts[2][i2];
                i_intrp(2) = stencil_box_crse.lower()(2) + i2;
#else
            const double wgt2 = 1.0;
#endif
#if (NDIM > 1)
                for (int i1 = 0; i1 <= degree; ++i1)
                {
                    const double& wgt1 = wgts[1][i1];
                    i_intrp(1) = stencil_box_crse.lower()(1) + i1;
#else
            const double wgt1 = 1.0;
#endif
                    for (int i0 = 0; i0 <= degree; ++i0)
                    {
                        const double& wgt0 = wgts[0][i0];
                        i_intrp(0) = stencil_box_crse.lower()(0) + i0;

                        (*fdata)(i_fine, d) += wgt0 * wgt1 * wgt2 * (*cdata)(i_intrp, d);
                    }
#if (NDIM > 1)
                }
#endif
#if (NDIM > 2)
            }
#endif
        }
    }
    return;
} // refine

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
