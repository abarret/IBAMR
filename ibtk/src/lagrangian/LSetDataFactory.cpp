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

#include "ibtk/LMarker.h"
#include "ibtk/LNode.h"
#include "ibtk/LNodeIndex.h"
#include "ibtk/LSet.h"
#include "ibtk/LSetData.h"
#include "ibtk/LSetDataFactory.h"
#include "ibtk/namespaces.h" // IWYU pragma: keep

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/pdat/CellGeometry.h"
#include "SAMRAI/pdat/IndexDataFactory.h"
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
LSetDataFactory<T>::LSetDataFactory(IntVector ghosts)
    : IndexDataFactory<LSet<T>, CellGeometry >(std::move(ghosts))
{
    // intentionally blank
    return;
} // LSetDataFactory

template <class T>
std::shared_ptr<PatchDataFactory >
LSetDataFactory<T>::cloneFactory(const IntVector& ghosts)
{
    return std::make_shared<LSetDataFactory<T>>(ghosts);
} // cloneFactory

template <class T>
std::shared_ptr<PatchData >
LSetDataFactory<T>::allocate(const Box& box, std::shared_ptr<Arena> pool) const
{
    if (!pool)
    {
        pool = ArenaManager::getManager()->getStandardAllocator();
    }
    PatchData* pd =
        std::make_shared<>(pool) LSetData<T>(box, IndexDataFactory<LSet<T>, CellGeometry >::getGhostCellWidth());
    return std::shared_ptr<PatchData >(pd, pool);
} // allocate

template <class T>
std::shared_ptr<PatchData >
LSetDataFactory<T>::allocate(const Patch& patch, std::shared_ptr<Arena> pool) const
{
    return allocate(patch.getBox(), pool);
} // allocate

template <class T>
size_t
LSetDataFactory<T>::getSizeOfMemory(const Box& /*box*/) const
{
    return Arena::align(sizeof(LSetData<T>));
} // getSizeOfMemory

template <class T>
bool
LSetDataFactory<T>::validCopyTo(const std::shared_ptr<PatchDataFactory >& dst_pdf) const
{
    std::shared_ptr<LSetDataFactory<T> > lnidf = dst_pdf;
    return lnidf;
} // validCopyTo

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

template class IBTK::LSetDataFactory<IBTK::LMarker>;
template class IBTK::LSetDataFactory<IBTK::LNode>;
template class IBTK::LSetDataFactory<IBTK::LNodeIndex>;

//////////////////////////////////////////////////////////////////////////////
