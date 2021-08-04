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

#include "ibtk/CartCellDoubleCubicCoarsen.h"
#include "ibtk/ibtk_utilities.h"
#include "ibtk/namespaces.h" // IWYU pragma: keep

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/geom/CartesianCellDoubleWeightedAverage.h"
#include "SAMRAI/pdat/CellData.h"

#include "SAMRAI/tbox/Utilities.h"

#include <ostream>
#include <string>

namespace SAMRAI
{
namespace hier
{

class Variable;
} // namespace hier
} // namespace SAMRAI

// FORTRAN ROUTINES
#if (NDIM == 2)
#define CC_CUBIC_COARSEN_FC IBTK_FC_FUNC(cccubiccoarsen2d, CCCUBICCOARSEN2D)
#endif
#if (NDIM == 3)
#define CC_CUBIC_COARSEN_FC IBTK_FC_FUNC(cccubiccoarsen3d, CCCUBICCOARSEN3D)
#endif

// Function interfaces
extern "C"
{
    void CC_CUBIC_COARSEN_FC(double* U_coarse,
                             const int& U_crse_gcw,
                             const double* U_fine0,
                             const int& U_fine_gcw,
                             const int& ilowerc0,
                             const int& iupperc0,
                             const int& ilowerc1,
                             const int& iupperc1,
#if (NDIM == 3)
                             const int& ilowerc2,
                             const int& iupperc2,
#endif
                             const int& ilowerf0,
                             const int& iupperf0,
                             const int& ilowerf1,
                             const int& iupperf1,
#if (NDIM == 3)
                             const int& ilowerf2,
                             const int& iupperf2,
#endif
                             const int* ratio_to_coarser,
                             const int* fblower,
                             const int* fbupper);
}

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

const std::string CartCellDoubleCubicCoarsen::s_op_name = "CUBIC_COARSEN";

namespace
{
static const int COARSEN_OP_PRIORITY = 0;
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

CartCellDoubleCubicCoarsen::CartCellDoubleCubicCoarsen() : CoarsenOperator(s_op_name)
{
    // intentionally blank
}

int
CartCellDoubleCubicCoarsen::getOperatorPriority() const
{
    return COARSEN_OP_PRIORITY;
} // getOperatorPriority

IntVector
CartCellDoubleCubicCoarsen::getStencilWidth(const SAMRAI::tbox::Dimension& dim) const
{
    return d_weighted_average_coarsen_op.getStencilWidth(dim);
} // getStencilWidth

void
CartCellDoubleCubicCoarsen::coarsen(Patch& coarse,
                                    const Patch& fine,
                                    const int dst_component,
                                    const int src_component,
                                    const Box& coarse_box,
                                    const IntVector& ratio) const
{
    if (ratio.min() < 4)
    {
        IBTK_DO_ONCE(TBOX_WARNING("CartCellDoubleCubicCoarsen::coarsen():\n"
                                  << "  cubic coarsening requires a refinement ratio of 4 or larger.\n"
                                  << "  reverting to weighted averaging." << std::endl););
        d_weighted_average_coarsen_op.coarsen(coarse, fine, dst_component, src_component, coarse_box, ratio);
        return;
    }
    std::shared_ptr<CellData<double> > cdata = std::static_pointer_cast<CellData<double> >(coarse.getPatchData(dst_component));
    std::shared_ptr<CellData<double> > fdata = std::static_pointer_cast<CellData<double> >(fine.getPatchData(src_component));
    const int U_fine_ghosts = (fdata->getGhostCellWidth()).max();
    const int U_crse_ghosts = (cdata->getGhostCellWidth()).max();
#if !defined(NDEBUG)
    if (U_fine_ghosts != (fdata->getGhostCellWidth()).min())
    {
        TBOX_ERROR("CartCellDoubleCubicCoarsen::coarsen():\n"
                   << "   fine patch data does not have uniform ghost cell widths" << std::endl);
    }
    if (U_crse_ghosts != (cdata->getGhostCellWidth()).min())
    {
        TBOX_ERROR("CartCellDoubleCubicCoarsen::coarsen():\n"
                   << "   coarse patch data does not have uniform ghost cell widths" << std::endl);
    }
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        if (ratio(d) % 2 == 1)
        {
            TBOX_ERROR("CartCellDoubleCubicCoarsen::coarsen():\n"
                       << "   refinement ratio between coarse and fine index spaces is odd" << std::endl);
        }
    }
#endif
    const int data_depth = cdata->getDepth();
#if !defined(NDEBUG)
    TBOX_ASSERT(data_depth == fdata->getDepth());
#endif
    const Box& patch_box_fine = fine.getBox();
    const Box& patch_box_crse = coarse.getBox();
    for (int depth = 0; depth < data_depth; ++depth)
    {
        double* const U_crse = cdata->getPointer(depth);
        const double* const U_fine = fdata->getPointer(depth);
        CC_CUBIC_COARSEN_FC(U_crse,
                            U_crse_ghosts,
                            U_fine,
                            U_fine_ghosts,
                            patch_box_crse.lower(0),
                            patch_box_crse.upper(0),
                            patch_box_crse.lower(1),
                            patch_box_crse.upper(1),
#if (NDIM == 3)
                            patch_box_crse.lower(2),
                            patch_box_crse.upper(2),
#endif
                            patch_box_fine.lower(0),
                            patch_box_fine.upper(0),
                            patch_box_fine.lower(1),
                            patch_box_fine.upper(1),
#if (NDIM == 3)
                            patch_box_fine.lower(2),
                            patch_box_fine.upper(2),
#endif
                            &(ratio(0)),
                            &(coarse_box.lower()(0)),
                            &(coarse_box.upper()(0)));
    }
    return;
} // coarsen

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
