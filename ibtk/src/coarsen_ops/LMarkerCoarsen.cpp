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

#include "ibtk/LMarkerCoarsen.h"
#include "ibtk/LMarkerSet.h"
#include "ibtk/LMarkerSetData.h"
#include "ibtk/LMarkerSetVariable.h"
#include "ibtk/namespaces.h" // IWYU pragma: keep

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/Patch.h"


#include <algorithm>
#include <string>
#include <vector>
#include <memory>

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

const std::string LMarkerCoarsen::s_op_name = "LMARKER_COARSEN";

namespace
{
static const int COARSEN_OP_PRIORITY = 0;
static const int COARSEN_OP_STENCIL_WIDTH = 0;

inline int
coarsen(const int index, const int ratio)
{
    return (index < 0 ? (index + 1) / ratio - 1 : index / ratio);
} // coarsen

inline hier::Index
coarsen_index(const hier::Index& i, const IntVector& ratio)
{
    hier::Index coarse_i(Dimension(NDIM));
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        coarse_i(d) = coarsen(i(d), ratio(d));
    }
    return coarse_i;
} // coarsen_index
} // namespace

/////////////////////////////// PUBLIC ///////////////////////////////////////

LMarkerCoarsen::LMarkerCoarsen() : CoarsenOperator(s_op_name)
{
    // intentionally blank
}

int
LMarkerCoarsen::getOperatorPriority() const
{
    return COARSEN_OP_PRIORITY;
} // getOperatorPriority

IntVector
LMarkerCoarsen::getStencilWidth(const SAMRAI::tbox::Dimension& /*dim*/) const
{
    return IntVector(Dimension(NDIM), COARSEN_OP_STENCIL_WIDTH);
} // getStencilWidth

void
LMarkerCoarsen::coarsen(Patch& coarse,
                        const Patch& fine,
                        const int dst_component,
                        const int src_component,
                        const Box& coarse_box,
                        const IntVector& ratio) const
{
    auto dst_mark_data = std::static_pointer_cast<LMarkerSetData>(coarse.getPatchData(dst_component));
    auto src_mark_data = std::static_pointer_cast<LMarkerSetData>(fine.getPatchData(src_component));

    const Box fine_box = Box::refine(coarse_box, ratio);
    LMarkerSetData::SetIterator it_end(*src_mark_data, false);
    for (LMarkerSetData::SetIterator it(*src_mark_data, true); it != it_end; it++)
    {
        const hier::Index& fine_i = it.getIndex();
        const hier::Index coarse_i = coarsen_index(fine_i, ratio);
        if (fine_box.contains(fine_i) && coarse_box.contains(coarse_i))
        {
            const LMarkerSet& fine_mark_set = *it;
            if (!dst_mark_data->isElement(coarse_i))
            {
                dst_mark_data->appendItemPointer(coarse_i, new LMarkerSet());
            }
            LMarkerSet& coarse_mark_set = *(dst_mark_data->getItem(coarse_i));
            coarse_mark_set.insert(coarse_mark_set.end(), fine_mark_set.begin(), fine_mark_set.end());
        }
    }
    return;
} // coarsen

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
