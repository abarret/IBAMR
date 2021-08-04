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

#include "ibtk/CartGridFunction.h"
#include "ibtk/namespaces.h" // IWYU pragma: keep

#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/Variable.h"

#include <string>
#include <utility>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

CartGridFunction::CartGridFunction(std::string object_name) : d_object_name(std::move(object_name))
{
    // intentionally blank
    return;
} // CartGridFunction

void
CartGridFunction::setDataOnPatchHierarchy(const int data_idx,
                                          std::shared_ptr<Variable > var,
                                          std::shared_ptr<PatchHierarchy > hierarchy,
                                          const double data_time,
                                          const bool initial_time,
                                          const int coarsest_ln_in,
                                          const int finest_ln_in)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(hierarchy);
#endif
    const int coarsest_ln = (coarsest_ln_in == -1 ? 0 : coarsest_ln_in);
    const int finest_ln = (finest_ln_in == -1 ? hierarchy->getFinestLevelNumber() : finest_ln_in);
    for (int level_num = coarsest_ln; level_num <= finest_ln; ++level_num)
    {
        setDataOnPatchLevel(data_idx, var, hierarchy->getPatchLevel(level_num), data_time, initial_time);
    }
    return;
} // setDataOnPatchHierarchy

void
CartGridFunction::setDataOnPatchLevel(const int data_idx,
                                      std::shared_ptr<Variable > var,
                                      std::shared_ptr<PatchLevel > level,
                                      const double data_time,
                                      const bool initial_time)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(level);
#endif
    for (PatchLevel::Iterator p = level->begin(); p != level->end(); p++)
    {
        setDataOnPatch(data_idx, var, *p, data_time, initial_time, level);
    }
    return;
} // setDataOnPatchLevel

//////////////////////////////////////////////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
