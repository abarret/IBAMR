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
#include "ibtk/LSetData.h"
#include "ibtk/namespaces.h" // IWYU pragma: keep

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/pdat/CellGeometry.h"
#include "SAMRAI/pdat/IndexData.h"
#include "SAMRAI/pdat/IndexDataFactory.h"
#include "SAMRAI/pdat/IndexVariable.h"
#include "SAMRAI/hier/IntVector.h"


#include <algorithm>
#include <ostream>
#include <utility>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

template <class T>
LSetData<T>::LSetData(Box box, IntVector ghosts)
    : IndexData<LSet<T>, CellGeometry >(std::move(box), std::move(ghosts))
{
    // intentionally blank
    return;
} // LSetData

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

#include "ibtk/LMarkerSet.h"

template class SAMRAI::pdat::IndexData<IBTK::LMarkerSet, CellGeometry >;
template class SAMRAI::pdat::IndexDataFactory<IBTK::LMarkerSet, CellGeometry >;
template class SAMRAI::pdat::IndexDataNode<IBTK::LMarkerSet, CellGeometry >;
template class SAMRAI::pdat::IndexIterator<IBTK::LMarkerSet, CellGeometry >;
template class SAMRAI::pdat::IndexVariable<IBTK::LMarkerSet, CellGeometry >;
template class IBTK::LSetData<IBTK::LMarker>;

#include "ibtk/LNodeSet.h"

template class SAMRAI::pdat::IndexData<IBTK::LNodeSet, CellGeometry >;
template class SAMRAI::pdat::IndexDataFactory<IBTK::LNodeSet, CellGeometry >;
template class SAMRAI::pdat::IndexDataNode<IBTK::LNodeSet, CellGeometry >;
template class SAMRAI::pdat::IndexIterator<IBTK::LNodeSet, CellGeometry >;
template class SAMRAI::pdat::IndexVariable<IBTK::LNodeSet, CellGeometry >;
template class IBTK::LSetData<IBTK::LNode>;

#include "ibtk/LNodeIndexSet.h"

template class SAMRAI::pdat::IndexData<IBTK::LNodeIndexSet, CellGeometry >;
template class SAMRAI::pdat::IndexDataFactory<IBTK::LNodeIndexSet, CellGeometry >;
template class SAMRAI::pdat::IndexDataNode<IBTK::LNodeIndexSet, CellGeometry >;
template class SAMRAI::pdat::IndexIterator<IBTK::LNodeIndexSet, CellGeometry >;
template class SAMRAI::pdat::IndexVariable<IBTK::LNodeIndexSet, CellGeometry >;
template class IBTK::LSetData<IBTK::LNodeIndex>;

//////////////////////////////////////////////////////////////////////////////
