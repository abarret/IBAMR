// ---------------------------------------------------------------------
//
// Copyright (c) 2026 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. See the COPYRIGHT file at the top level of the distribution.
//
// ---------------------------------------------------------------------

#ifndef included_IBAMR_IBAngularSpringForceSpec_inl
#define included_IBAMR_IBAngularSpringForceSpec_inl

#include <ibamr/config.h>

#include <ibamr/IBAngularSpringForceSpec.h>

#include <ibtk/StreamableManager.h>

#include <tbox/Utilities.h>

namespace IBAMR
{
inline bool
IBAngularSpringForceSpec::getIsRegisteredWithStreamableManager()
{
    return STREAMABLE_CLASS_ID != IBTK::StreamableManager::getUnregisteredID();
}

inline IBAngularSpringForceSpec::IBAngularSpringForceSpec(const unsigned int num_springs)
    : d_neighbor_idxs(num_springs), d_stiffnesses(num_springs), d_rest_angles(num_springs)
{
#if !defined(NDEBUG)
    if (!getIsRegisteredWithStreamableManager())
    {
        TBOX_ERROR("IBAngularSpringForceSpec::IBAngularSpringForceSpec():\n"
                   << "  must call IBAngularSpringForceSpec::registerWithStreamableManager() before\n"
                   << "  creating any IBAngularSpringForceSpec objects.\n");
    }
#endif
}

inline IBAngularSpringForceSpec::IBAngularSpringForceSpec(const int master_idx,
                                                          const std::vector<NeighborIdxs>& neighbor_idxs,
                                                          const std::vector<double>& stiffnesses,
                                                          const std::vector<double>& rest_angles)
    : d_master_idx(master_idx), d_neighbor_idxs(neighbor_idxs), d_stiffnesses(stiffnesses), d_rest_angles(rest_angles)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_neighbor_idxs.size() == d_stiffnesses.size());
    TBOX_ASSERT(d_neighbor_idxs.size() == d_rest_angles.size());
    if (!getIsRegisteredWithStreamableManager())
    {
        TBOX_ERROR("IBAngularSpringForceSpec::IBAngularSpringForceSpec():\n"
                   << "  must call IBAngularSpringForceSpec::registerWithStreamableManager() before\n"
                   << "  creating any IBAngularSpringForceSpec objects.\n");
    }
#endif
}

inline IBAngularSpringForceSpec::~IBAngularSpringForceSpec() = default;

inline unsigned int
IBAngularSpringForceSpec::getNumberOfAngularSprings() const
{
    return static_cast<unsigned int>(d_neighbor_idxs.size());
}

inline const int&
IBAngularSpringForceSpec::getMasterNodeIndex() const
{
    return d_master_idx;
}

inline int&
IBAngularSpringForceSpec::getMasterNodeIndex()
{
    return d_master_idx;
}

inline const std::vector<IBAngularSpringForceSpec::NeighborIdxs>&
IBAngularSpringForceSpec::getNeighborNodeIndices() const
{
    return d_neighbor_idxs;
}

inline std::vector<IBAngularSpringForceSpec::NeighborIdxs>&
IBAngularSpringForceSpec::getNeighborNodeIndices()
{
    return d_neighbor_idxs;
}

inline const std::vector<double>&
IBAngularSpringForceSpec::getStiffnesses() const
{
    return d_stiffnesses;
}

inline std::vector<double>&
IBAngularSpringForceSpec::getStiffnesses()
{
    return d_stiffnesses;
}

inline const std::vector<double>&
IBAngularSpringForceSpec::getRestAngles() const
{
    return d_rest_angles;
}

inline std::vector<double>&
IBAngularSpringForceSpec::getRestAngles()
{
    return d_rest_angles;
}

inline int
IBAngularSpringForceSpec::getStreamableClassID() const
{
    return STREAMABLE_CLASS_ID;
}

inline size_t
IBAngularSpringForceSpec::getDataStreamSize() const
{
    const size_t n = d_neighbor_idxs.size();
    return (2 + 2 * n) * SAMRAI::tbox::AbstractStream::sizeofInt() +
           2 * n * SAMRAI::tbox::AbstractStream::sizeofDouble();
}

inline void
IBAngularSpringForceSpec::packStream(SAMRAI::tbox::AbstractStream& stream)
{
    const auto n = static_cast<unsigned int>(d_neighbor_idxs.size());
    std::vector<int> neighbors(2 * n);
    for (unsigned int k = 0; k < n; ++k)
    {
        neighbors[2 * k] = d_neighbor_idxs[k].first;
        neighbors[2 * k + 1] = d_neighbor_idxs[k].second;
    }
    stream << static_cast<int>(n);
    stream.pack(&d_master_idx, 1);
    if (n > 0)
    {
        stream.pack(neighbors.data(), 2 * n);
        stream.pack(d_stiffnesses.data(), n);
        stream.pack(d_rest_angles.data(), n);
    }
}
} // namespace IBAMR

#endif // #ifndef included_IBAMR_IBAngularSpringForceSpec_inl
