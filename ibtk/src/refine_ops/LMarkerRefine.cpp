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

#include "ibtk/IndexUtilities.h"
#include "ibtk/LMarker.h"
#include "ibtk/LMarkerRefine.h"
#include "ibtk/LMarkerSet.h"
#include "ibtk/LMarkerSetData.h"
#include "ibtk/LMarkerSetVariable.h"
#include "ibtk/LSetData.h"
#include "ibtk/ibtk_utilities.h"
#include "ibtk/namespaces.h" // IWYU pragma: keep

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/tbox/MathUtilities.h"


#include <algorithm>
#include <array>
#include <string>

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

const std::string LMarkerRefine::s_op_name = "LMARKER_REFINE";

namespace
{
static const int REFINE_OP_PRIORITY = 0;
static const int REFINE_OP_STENCIL_WIDTH = 0;
} // namespace

/////////////////////////////// PUBLIC ///////////////////////////////////////

bool
LMarkerRefine::findRefineOperator(const std::shared_ptr<Variable >& var, const std::string& op_name) const
{
    std::shared_ptr<LMarkerSetVariable> mark_var = var;
    return (mark_var && op_name == s_op_name);
} // findRefineOperator

const std::string&
LMarkerRefine::getOperatorName() const
{
    return s_op_name;
} // getOperatorName

int
LMarkerRefine::getOperatorPriority() const
{
    return REFINE_OP_PRIORITY;
} // getOperatorPriority

IntVector
LMarkerRefine::getStencilWidth() const
{
    return REFINE_OP_STENCIL_WIDTH;
} // getStencilWidth

void
LMarkerRefine::refine(Patch& fine,
                      const Patch& coarse,
                      const int dst_component,
                      const int src_component,
                      const Box& fine_box,
                      const IntVector& ratio) const
{
    std::shared_ptr<LMarkerSetData> dst_mark_data = std::static_pointer_cast<LMarkerSetData>(fine.getPatchData(dst_component));
    std::shared_ptr<LMarkerSetData> src_mark_data = std::static_pointer_cast<LMarkerSetData>(coarse.getPatchData(src_component));

    const Box& fine_patch_box = fine.getBox();
    const std::shared_ptr<CartesianPatchGeometry > fine_patch_geom = std::static_pointer_cast<CartesianPatchGeometry >(fine.getPatchGeometry());
    const hier::Index& fine_patch_lower = fine_patch_box.lower();
    const hier::Index& fine_patch_upper = fine_patch_box.upper();
    const double* const fine_patchXLower = fine_patch_geom->getXLower();
    const double* const fine_patchXUpper = fine_patch_geom->getXUpper();

    const std::shared_ptr<CartesianPatchGeometry > coarse_patch_geom = std::static_pointer_cast<CartesianPatchGeometry >(coarse.getPatchGeometry());
    const double* const coarse_patchDx = coarse_patch_geom->getDx();

    const Box coarse_box = Box::coarsen(fine_box, ratio);
    for (LMarkerSetData::SetIterator it(*src_mark_data); it; it++)
    {
        const hier::Index& coarse_i = it.getIndex();
        if (coarse_box.contains(coarse_i))
        {
            const LMarkerSet& coarse_mark_set = it();
            for (const auto& coarse_mark : coarse_mark_set)
            {
                const Point& X = coarse_mark->getPosition();
                const IntVector& offset = coarse_mark->getPeriodicOffset();
                std::array<double, NDIM> X_shifted;
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    X_shifted[d] = X[d] + static_cast<double>(offset(d)) * coarse_patchDx[d];
                }
                hier::Index fine_i = IndexUtilities::getCellIndex(X_shifted, fine_patch_geom, fine_patch_box);

                // Catch edge cases in which roundoff error can cause problems.
                //
                // NOTE: This can permit markers to "escape" out the "top" of the domain if
                // the marker position is equal to the upper domain extent.  The marker
                // advection code needs to keep markers from hitting the domain boundaries.
                // (Note that bad things also happen if IB points hit the domain boundaries,
                // so this is not an issue that is unique to markers.)
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    if (MathUtilities<double>::equalEps(X_shifted[d], fine_patchXLower[d]))
                    {
                        X_shifted[d] = fine_patchXLower[d];
                        fine_i(d) = fine_patch_lower(d);
                    }
                    else if (MathUtilities<double>::equalEps(X_shifted[d], fine_patchXUpper[d]))
                    {
                        X_shifted[d] = fine_patchXUpper[d];
                        fine_i(d) = fine_patch_upper(d) + 1;
                    }
                }

                if (fine_box.contains(fine_i))
                {
                    if (!dst_mark_data->isElement(fine_i))
                    {
                        dst_mark_data->appendItemPointer(fine_i, std::make_shared<LMarkerSet>());
                    }
                    LMarkerSet& fine_mark_set = *(dst_mark_data->getItem(fine_i));
                    fine_mark_set.push_back(coarse_mark);
                }
            }
        }
    }
    return;
} // refine

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
