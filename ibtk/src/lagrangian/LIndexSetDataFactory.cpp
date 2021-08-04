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

#include "ibtk/LIndexSetData.h"
#include "ibtk/LIndexSetDataFactory.h"
#include "ibtk/LNode.h"
#include "ibtk/LNodeIndex.h"
#include "ibtk/LSetDataFactory.h"
#include "ibtk/namespaces.h" // IWYU pragma: keep

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/PatchData.h"
#include "SAMRAI/hier/PatchDataFactory.h"
#include "SAMRAI/tbox/Arena.h"
#include "SAMRAI/tbox/ArenaManager.h"


#include <algorithm>
#include <utility>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

template <class T>
LIndexSetDataFactory<T>::LIndexSetDataFactory(IntVector ghosts) : LSetDataFactory<T>(std::move(ghosts))
{
    // intentionally blank
    return;
} // LIndexSetDataFactory

template <class T>
std::shared_ptr<PatchDataFactory >
LIndexSetDataFactory<T>::cloneFactory(const IntVector& ghosts)
{
    return std::make_shared<LIndexSetDataFactory<T>>(ghosts);
} // cloneFactory

template <class T>
std::shared_ptr<PatchData >
LIndexSetDataFactory<T>::allocate(const Box& box, std::shared_ptr<Arena> pool) const
{
    if (!pool)
    {
        pool = ArenaManager::getManager()->getStandardAllocator();
    }
    PatchData* pd = std::make_shared<>(pool) LIndexSetData<T>(box, LSetDataFactory<T>::getGhostCellWidth());
    return std::shared_ptr<PatchData >(pd, pool);
} // allocate

template <class T>
std::shared_ptr<PatchData >
LIndexSetDataFactory<T>::allocate(const Patch& patch, std::shared_ptr<Arena> pool) const
{
    return allocate(patch.getBox(), pool);
} // allocate

template <class T>
size_t
LIndexSetDataFactory<T>::getSizeOfMemory(const Box& /*box*/) const
{
    return Arena::align(sizeof(LIndexSetData<T>));
} // getSizeOfMemory

template <class T>
bool
LIndexSetDataFactory<T>::validCopyTo(const std::shared_ptr<PatchDataFactory >& dst_pdf) const
{
    const std::shared_ptr<LIndexSetDataFactory<T> > lnidf = dst_pdf;
    return lnidf;
} // validCopyTo

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

template class IBTK::LIndexSetDataFactory<IBTK::LNode>;
template class IBTK::LIndexSetDataFactory<IBTK::LNodeIndex>;

//////////////////////////////////////////////////////////////////////////////
