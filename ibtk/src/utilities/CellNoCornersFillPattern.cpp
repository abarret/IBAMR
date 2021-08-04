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

#include "ibtk/CellNoCornersFillPattern.h"
#include "ibtk/namespaces.h" // IWYU pragma: keep

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/BoxGeometry.h"
#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/hier/BoxOverlap.h"
#include "SAMRAI/pdat/CellGeometry.h"
#include "SAMRAI/pdat/CellOverlap.h"
#include "SAMRAI/hier/IntVector.h"


#include <ostream>
#include <string>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
static const std::string PATTERN_NAME = "CELL_NO_CORNERS_FILL_PATTERN";
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

CellNoCornersFillPattern::CellNoCornersFillPattern(const int stencil_width,
                                                   const bool include_dst_patch_box,
                                                   const bool include_edges_on_dst_level,
                                                   const bool include_edges_on_src_level)
    : d_stencil_width(stencil_width),
      d_include_dst_patch_box(include_dst_patch_box),
      d_include_edges_on_dst_level(include_edges_on_dst_level),
      d_include_edges_on_src_level(include_edges_on_src_level)
{
    // intentionally blank
    return;
} // CellNoCornersFillPattern

std::shared_ptr<BoxOverlap >
CellNoCornersFillPattern::calculateOverlap(const BoxGeometry& dst_geometry,
                                           const BoxGeometry& src_geometry,
                                           const Box& /*dst_patch_box*/,
                                           const Box& src_mask,
                                           const bool overwrite_interior,
                                           const IntVector& src_offset) const
{
    std::shared_ptr<CellOverlap > box_geom_overlap =
        dst_geometry.calculateOverlap(src_geometry, src_mask, overwrite_interior, src_offset);
#if !defined(NDEBUG)
    TBOX_ASSERT(box_geom_overlap);
#endif
    auto const t_dst_geometry = dynamic_cast<const CellGeometry*>(&dst_geometry);
#if !defined(NDEBUG)
    TBOX_ASSERT(t_dst_geometry);
#endif
    BoxList dst_boxes;
    if (!box_geom_overlap->isOverlapEmpty())
    {
        const Box& dst_box = t_dst_geometry->getBox();
        const BoxList& box_geom_overlap_boxes = box_geom_overlap->getDestinationBoxList();

        // Determine the stencil boxes with the specified ghost cell width.
        BoxList stencil_boxes;
        if (NDIM == 2 || (!d_include_edges_on_src_level && !d_include_edges_on_dst_level))
        {
            for (unsigned int i = 0; i < NDIM; ++i)
            {
                Box box = dst_box;
                box.lower()(i) -= d_stencil_width(i);
                box.upper()(i) += d_stencil_width(i);
                stencil_boxes.appendItem(box);
            }
        }
        else
        {
            for (unsigned int j = 0; j < NDIM; ++j)
            {
                for (unsigned int i = 0; i < NDIM; ++i)
                {
                    if (i == j) continue;
                    Box box = dst_box;
                    box.lower()(i) -= d_stencil_width(i);
                    box.upper()(i) += d_stencil_width(i);
                    box.lower()(j) -= d_stencil_width(j);
                    box.upper()(j) += d_stencil_width(j);
                    stencil_boxes.appendItem(box);
                }
            }
        }

        // Intersect the overlap boxes with the stencil boxes.
        for (BoxList::Iterator it1(box_geom_overlap_boxes); it1; it1++)
        {
            BoxList overlap_boxes(stencil_boxes);
            overlap_boxes.intersectBoxes(it1());
            for (BoxList::Iterator it2(overlap_boxes); it2; it2++)
            {
                const Box& overlap_box = it2();
                if (!overlap_box.empty()) dst_boxes.appendItem(overlap_box);
            }
        }
    }
    return std::make_shared<CellOverlap>(dst_boxes, src_offset);
} // calculateOverlap

std::shared_ptr<BoxOverlap >
CellNoCornersFillPattern::calculateOverlapOnLevel(const BoxGeometry& dst_geometry,
                                                  const BoxGeometry& src_geometry,
                                                  const Box& dst_patch_box,
                                                  const Box& src_mask,
                                                  const bool overwrite_interior,
                                                  const IntVector& src_offset,
                                                  const int dst_level_num,
                                                  const int /*src_level_num*/) const
{
    std::shared_ptr<CellOverlap > box_geom_overlap =
        dst_geometry.calculateOverlap(src_geometry, src_mask, overwrite_interior, src_offset);
#if !defined(NDEBUG)
    TBOX_ASSERT(box_geom_overlap);
#endif
    auto const t_dst_geometry = dynamic_cast<const CellGeometry*>(&dst_geometry);
#if !defined(NDEBUG)
    TBOX_ASSERT(t_dst_geometry);
#endif
    BoxList dst_boxes;
    if (!box_geom_overlap->isOverlapEmpty())
    {
        const Box& dst_box = t_dst_geometry->getBox();
        const BoxList& box_geom_overlap_boxes = box_geom_overlap->getDestinationBoxList();

        // Determine the stencil boxes with the specified ghost cell width.
        BoxList stencil_boxes;
        if (NDIM == 2 || (!d_include_edges_on_dst_level && dst_level_num == d_target_level_num) ||
            (!d_include_edges_on_src_level && dst_level_num != d_target_level_num))
        {
            for (unsigned int i = 0; i < NDIM; ++i)
            {
                Box box = dst_box;
                box.lower()(i) -= d_stencil_width(i);
                box.upper()(i) += d_stencil_width(i);
                stencil_boxes.appendItem(box);
            }
        }
        else
        {
            for (unsigned int j = 0; j < NDIM; ++j)
            {
                for (unsigned int i = 0; i < NDIM; ++i)
                {
                    if (i == j) continue;
                    Box box = dst_box;
                    box.lower()(i) -= d_stencil_width(i);
                    box.upper()(i) += d_stencil_width(i);
                    box.lower()(j) -= d_stencil_width(j);
                    box.upper()(j) += d_stencil_width(j);
                    stencil_boxes.appendItem(box);
                }
            }
        }

        // Intersect the overlap boxes with the stencil boxes.
        for (BoxList::Iterator it1(box_geom_overlap_boxes); it1; it1++)
        {
            BoxList overlap_boxes(stencil_boxes);
            overlap_boxes.intersectBoxes(it1());
            if (dst_level_num == d_target_level_num && !d_include_dst_patch_box)
            {
                overlap_boxes.removeIntersections(dst_patch_box);
            }
            for (BoxList::Iterator it2(overlap_boxes); it2; it2++)
            {
                const Box& overlap_box = it2();
                if (!overlap_box.empty()) dst_boxes.appendItem(overlap_box);
            }
        }
    }
    return std::make_shared<CellOverlap>(dst_boxes, src_offset);
} // calculateOverlapOnLevel

void
CellNoCornersFillPattern::setTargetPatchLevelNumber(const int level_num)
{
    d_target_level_num = level_num;
    return;
} // setTargetPatchLevelNumber

IntVector&
CellNoCornersFillPattern::getStencilWidth()
{
    return d_stencil_width;
} // getStencilWidth

const std::string&
CellNoCornersFillPattern::getPatternName() const
{
    return PATTERN_NAME;
} // getPatternName

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
