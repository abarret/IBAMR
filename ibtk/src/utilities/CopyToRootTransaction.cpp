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

#include "ibtk/CopyToRootTransaction.h"
#include "ibtk/namespaces.h" // IWYU pragma: keep

#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/hier/BoxGeometry.h"
#include "SAMRAI/hier/BoxOverlap.h"
#include "SAMRAI/geom/GridGeometry.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/PatchData.h"
#include "SAMRAI/hier/PatchDataFactory.h"
#include "SAMRAI/hier/PatchDescriptor.h"
#include "SAMRAI/hier/PatchLevel.h"


#include <ostream>

namespace SAMRAI
{
namespace hier
{

class Box;
} // namespace hier
} // namespace SAMRAI

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

CopyToRootTransaction::CopyToRootTransaction(const int src_proc,
                                             const int dst_proc,
                                             std::shared_ptr<PatchLevel > patch_level,
                                             const int src_patch_data_idx,
                                             std::shared_ptr<PatchData > dst_patch_data)
    : d_src_proc(src_proc),
      d_dst_proc(dst_proc),
      d_patch_level(patch_level),
      d_src_patch_data_idx(src_patch_data_idx),
      d_dst_patch_data(dst_patch_data)
{
    // intentionally blank
    return;
} // CopyToRootTransaction

std::shared_ptr<PatchData >
CopyToRootTransaction::getRootPatchData() const
{
    return d_dst_patch_data;
} // getRootPatchData

bool
CopyToRootTransaction::canEstimateIncomingMessageSize()
{
    return false;
} // canEstimateIncomingMessageSize

int
CopyToRootTransaction::computeIncomingMessageSize()
{
    return 0;
} // computeIncomingMessageSize

int
CopyToRootTransaction::computeOutgoingMessageSize()
{
    std::shared_ptr<PatchDataFactory > pdat_factory =
        d_patch_level->getPatchDescriptor()->getPatchDataFactory(d_src_patch_data_idx);

    std::shared_ptr<GridGeometry > grid_geom = d_patch_level->getGridGeometry();
#if !defined(NDEBUG)
    TBOX_ASSERT(grid_geom->getDomainIsSingleBox());
#endif
    const Box& dst_box = grid_geom->getPhysicalDomain()[0];
    std::shared_ptr<BoxGeometry > dst_box_geometry = pdat_factory->getBoxGeometry(dst_box);

    int size = MessageStream::sizeofInt();
    for (PatchLevel::Iterator p(d_patch_level); p != level->end(); p++)
    {
        const int src_patch_num = p();
        size += MessageStream::sizeofInt();
        std::shared_ptr<Patch > patch = d_patch_level->getPatch(src_patch_num);
        const Box& src_box = patch->getBox();
        std::shared_ptr<BoxGeometry > src_box_geometry = pdat_factory->getBoxGeometry(src_box);
        const Box& src_mask = dst_box;
        const bool overwrite_interior = true;
        const IntVector src_shift = 0;
        std::shared_ptr<BoxOverlap > box_overlap =
            dst_box_geometry->calculateOverlap(*src_box_geometry, src_mask, overwrite_interior, src_shift);
        size += std::static_pointer_cast<>(patch->getPatchData(d_src_patch_data_idx))->getDataStreamSize(*box_overlap);
    }
    return size;
} // computeOutgoingMessageSize

int
CopyToRootTransaction::getSourceProcessor()
{
    return d_src_proc;
} // getSourceProcessor

int
CopyToRootTransaction::getDestinationProcessor()
{
    return d_dst_proc;
} // getDestinationProcessor

void
CopyToRootTransaction::packStream(MessageStream& stream)
{
    std::shared_ptr<PatchDataFactory > pdat_factory =
        d_patch_level->getPatchDescriptor()->getPatchDataFactory(d_src_patch_data_idx);

    std::shared_ptr<GridGeometry > grid_geom = d_patch_level->getGridGeometry();
#if !defined(NDEBUG)
    TBOX_ASSERT(grid_geom->getDomainIsSingleBox());
#endif
    const Box& dst_box = grid_geom->getPhysicalDomain()[0];
    std::shared_ptr<BoxGeometry > dst_box_geometry = pdat_factory->getBoxGeometry(dst_box);

    int src_patch_count = 0;
    for (PatchLevel::Iterator p(d_patch_level); p != level->end(); p++)
    {
        ++src_patch_count;
    }
    stream << src_patch_count;

    for (PatchLevel::Iterator p(d_patch_level); p != level->end(); p++)
    {
        const int src_patch_num = p();
        stream << src_patch_num;
        std::shared_ptr<Patch > patch = d_patch_level->getPatch(src_patch_num);
        const Box& src_box = patch->getBox();
        std::shared_ptr<BoxGeometry > src_box_geometry = pdat_factory->getBoxGeometry(src_box);
        const Box& src_mask = dst_box;
        const bool overwrite_interior = true;
        const IntVector src_shift = 0;
        std::shared_ptr<BoxOverlap > box_overlap =
            dst_box_geometry->calculateOverlap(*src_box_geometry, src_mask, overwrite_interior, src_shift);
        std::static_pointer_cast<>(patch->getPatchData(d_src_patch_data_idx))->packStream(stream, *box_overlap);
    }
    return;
} // packStream

void
CopyToRootTransaction::unpackStream(MessageStream& stream)
{
    std::shared_ptr<PatchDataFactory > pdat_factory =
        d_patch_level->getPatchDescriptor()->getPatchDataFactory(d_src_patch_data_idx);

    std::shared_ptr<GridGeometry > grid_geom = d_patch_level->getGridGeometry();
#if !defined(NDEBUG)
    TBOX_ASSERT(grid_geom->getDomainIsSingleBox());
#endif
    const Box& dst_box = grid_geom->getPhysicalDomain()[0];
    std::shared_ptr<BoxGeometry > dst_box_geometry = pdat_factory->getBoxGeometry(dst_box);

    int src_patch_count;
    stream >> src_patch_count;
    for (int p = 0; p < src_patch_count; ++p)
    {
        int src_patch_num;
        stream >> src_patch_num;
        const Box& src_box = d_patch_level->getBoxes()[src_patch_num];
        std::shared_ptr<BoxGeometry > src_box_geometry = pdat_factory->getBoxGeometry(src_box);
        const Box& src_mask = dst_box;
        const bool overwrite_interior = true;
        const IntVector src_shift = 0;
        std::shared_ptr<BoxOverlap > box_overlap =
            dst_box_geometry->calculateOverlap(*src_box_geometry, src_mask, overwrite_interior, src_shift);
        d_dst_patch_data->unpackStream(stream, *box_overlap);
    }
    return;
} // unpackStream

void
CopyToRootTransaction::copyLocalData()
{
    std::shared_ptr<PatchDataFactory > pdat_factory =
        d_patch_level->getPatchDescriptor()->getPatchDataFactory(d_src_patch_data_idx);

    std::shared_ptr<GridGeometry > grid_geom = d_patch_level->getGridGeometry();
#if !defined(NDEBUG)
    TBOX_ASSERT(grid_geom->getDomainIsSingleBox());
#endif
    const Box& dst_box = grid_geom->getPhysicalDomain()[0];
    std::shared_ptr<BoxGeometry > dst_box_geometry = pdat_factory->getBoxGeometry(dst_box);

    for (PatchLevel::Iterator p(d_patch_level); p != level->end(); p++)
    {
        int src_patch_num = p();
        std::shared_ptr<Patch > patch = d_patch_level->getPatch(src_patch_num);
        const Box& src_box = patch->getBox();
        std::shared_ptr<BoxGeometry > src_box_geometry = pdat_factory->getBoxGeometry(src_box);
        const Box& src_mask = dst_box;
        const bool overwrite_interior = true;
        const IntVector src_shift = 0;
        std::shared_ptr<BoxOverlap > box_overlap =
            dst_box_geometry->calculateOverlap(*src_box_geometry, src_mask, overwrite_interior, src_shift);
        d_dst_patch_data->copy(*std::static_pointer_cast<>(patch->getPatchData(d_src_patch_data_idx)), *box_overlap);
    }
    return;
} // copyLocalData

void
CopyToRootTransaction::printClassData(std::ostream& stream) const
{
    stream << "CopyToRootTransaction::printClassData() is not implemented\n";
    return;
} // printClassData

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
