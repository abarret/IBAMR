// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2018 by the IBAMR developers
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

#include "ibtk/RefinePatchStrategySet.h"
#include "ibtk/namespaces.h" // IWYU pragma: keep

namespace SAMRAI
{
namespace hier
{

class Patch;
} // namespace hier
} // namespace SAMRAI

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

RefinePatchStrategySet::~RefinePatchStrategySet()
{
    if (d_managed)
    {
        for (const auto& strategy : d_strategy_set)
        {
            delete strategy;
        }
    }
    return;
} // ~RefinePatchStrategySet

void
RefinePatchStrategySet::setPhysicalBoundaryConditions(Patch& patch,
                                                      const double fill_time,
                                                      const IntVector& ghost_width_to_fill)
{
    for (const auto& strategy : d_strategy_set)
    {
        strategy->setPhysicalBoundaryConditions(patch, fill_time, ghost_width_to_fill);
    }
    return;
} // setPhysicalBoundaryConditions

IntVector
RefinePatchStrategySet::getRefineOpStencilWidth(const SAMRAI::tbox::Dimension& dim) const
{
    auto width = IntVector::getZero(dim);
    for (const auto& strategy : d_strategy_set)
    {
        width = IntVector::max(width, strategy->getRefineOpStencilWidth(dim));
    }
    return width;
} // getRefineOpStencilWidth()

void
RefinePatchStrategySet::preprocessRefine(Patch& fine,
                                         const Patch& coarse,
                                         const Box& fine_box,
                                         const IntVector& ratio)
{
    for (const auto& strategy : d_strategy_set)
    {
        strategy->preprocessRefine(fine, coarse, fine_box, ratio);
    }
    return;
} // preprocessRefine

void
RefinePatchStrategySet::postprocessRefine(Patch& fine,
                                          const Patch& coarse,
                                          const Box& fine_box,
                                          const IntVector& ratio)
{
    for (const auto& strategy : d_strategy_set)
    {
        strategy->postprocessRefine(fine, coarse, fine_box, ratio);
    }
    return;
} // postprocessRefine

void
RefinePatchStrategySet::preprocessRefineBoxes(Patch& fine,
                                              const Patch& coarse,
                                              const BoxList& fine_boxes,
                                              const IntVector& ratio)
{
    for (const auto& strategy : d_strategy_set)
    {
        strategy->preprocessRefineBoxes(fine, coarse, fine_boxes, ratio);
    }
    return;
} // preprocessRefineBoxes

void
RefinePatchStrategySet::postprocessRefineBoxes(Patch& fine,
                                               const Patch& coarse,
                                               const BoxList& fine_boxes,
                                               const IntVector& ratio)
{
    for (const auto& strategy : d_strategy_set)
    {
        strategy->postprocessRefineBoxes(fine, coarse, fine_boxes, ratio);
    }
    return;
} // postprocessRefineBoxes

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
