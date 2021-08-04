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

#include "ibtk/NodeSynchCopyFillPattern.h"
#include "ibtk/namespaces.h" // IWYU pragma: keep

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/BoxGeometry.h"
#include "SAMRAI/pdat/NodeGeometry.h"
#include "SAMRAI/pdat/NodeOverlap.h"


#include <string>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
static const std::string PATTERN_NAME = "NODE_SYNCH_COPY_FILL_PATTERN";
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

NodeSynchCopyFillPattern::NodeSynchCopyFillPattern(const unsigned int axis) : d_axis(axis)
{
    // intentionally blank
    return;
} // NodeSynchCopyFillPattern

std::shared_ptr<BoxOverlap >
NodeSynchCopyFillPattern::calculateOverlap(const BoxGeometry& dst_geometry,
                                           const BoxGeometry& src_geometry,
                                           const Box& /*dst_patch_box*/,
                                           const Box& src_mask,
                                           const bool overwrite_interior,
                                           const IntVector& src_offset) const
{
    std::shared_ptr<NodeOverlap > box_geom_overlap =
        dst_geometry.calculateOverlap(src_geometry, src_mask, overwrite_interior, src_offset);
#if !defined(NDEBUG)
    TBOX_ASSERT(box_geom_overlap);
#endif
    if (box_geom_overlap->isOverlapEmpty()) return box_geom_overlap;

    auto const t_dst_geometry = dynamic_cast<const NodeGeometry*>(&dst_geometry);
#if !defined(NDEBUG)
    TBOX_ASSERT(t_dst_geometry);
#endif
    BoxList dst_boxes;
    bool skip = false;
    for (unsigned int d = 0; d < NDIM && !skip; ++d)
    {
        if (d != d_axis)
        {
            skip = skip || (src_offset(d) != 0);
        }
    }
    if (!skip)
    {
        // Determine the stencil box.
        const Box& dst_box = t_dst_geometry->getBox();
        Box stencil_box = NodeGeometry::toNodeBox(dst_box);
        stencil_box.lower(d_axis) = stencil_box.upper(d_axis);

        // Intersect the original overlap boxes with the stencil box.
        const BoxList& box_geom_overlap_boxes = box_geom_overlap->getDestinationBoxList();
        for (BoxList::Iterator it(box_geom_overlap_boxes); it; it++)
        {
            const Box overlap_box = stencil_box * it();
            if (!overlap_box.empty()) dst_boxes.appendItem(overlap_box);
        }
    }
    return std::make_shared<NodeOverlap>(dst_boxes, src_offset);
} // calculateOverlap

IntVector&
NodeSynchCopyFillPattern::getStencilWidth()
{
    return d_stencil_width;
} // getStencilWidth

const std::string&
NodeSynchCopyFillPattern::getPatternName() const
{
    return PATTERN_NAME;
} // getPatternName

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
