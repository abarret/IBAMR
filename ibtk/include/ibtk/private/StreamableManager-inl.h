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

#ifndef included_IBTK_StreamableManager_inl_h
#define included_IBTK_StreamableManager_inl_h

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "ibtk/Streamable.h"
#include "ibtk/StreamableManager.h"

#include "SAMRAI/tbox/MessageStream.h"
#include "SAMRAI/tbox/Utilities.h"

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// PUBLIC ///////////////////////////////////////

inline size_t
StreamableManager::getDataStreamSize(const std::shared_ptr<Streamable> data_item) const
{
    return SAMRAI::tbox::MessageStream::getSizeof<int>() + data_item->getDataStreamSize();
} // getDataStreamSize

inline size_t
StreamableManager::getDataStreamSize(const std::vector<std::shared_ptr<Streamable> >& data_items) const
{
    size_t size = SAMRAI::tbox::MessageStream::getSizeof<int>();
    for (const auto& data_item : data_items)
    {
        size += getDataStreamSize(data_item);
    }
    return size;
} // getDataStreamSize

inline void
StreamableManager::packStream(SAMRAI::tbox::MessageStream& stream, std::shared_ptr<Streamable> data_item)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(data_item);
#endif
    const int streamable_id = data_item->getStreamableClassID();
    stream.pack(&streamable_id, 1);
    data_item->packStream(stream);
    return;
} // packStream

inline void
StreamableManager::packStream(SAMRAI::tbox::MessageStream& stream,
                              std::vector<std::shared_ptr<Streamable> >& data_items)
{
    const int num_data = static_cast<int>(data_items.size());
    stream.pack(&num_data, 1);
    for (auto& data_item : data_items)
    {
        packStream(stream, data_item);
    }
    return;
} // packStream

inline std::shared_ptr<Streamable>
StreamableManager::unpackStream(SAMRAI::tbox::MessageStream& stream, const SAMRAI::hier::IntVector& offset)
{
    int streamable_id;
    stream.unpack(&streamable_id, 1);
#if !defined(NDEBUG)
    TBOX_ASSERT(d_factory_map.count(streamable_id) == 1);
#endif
    return d_factory_map[streamable_id]->unpackStream(stream, offset);
} // unpackStream

inline void
StreamableManager::unpackStream(SAMRAI::tbox::MessageStream& stream,
                                const SAMRAI::hier::IntVector& offset,
                                std::vector<std::shared_ptr<Streamable> >& data_items)
{
    int num_data;
    stream.unpack(&num_data, 1);
    data_items.resize(num_data);
    for (auto& data_item : data_items)
    {
        data_item = unpackStream(stream, offset);
    }
    std::vector<std::shared_ptr<Streamable> >(data_items).swap(data_items); // trim-to-fit
    return;
} // unpackStream

//////////////////////////////////////////////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBTK_StreamableManager_inl_h
