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

#include <ibamr/IBAngularSpringForceSpec.h>

#include <ibtk/Streamable.h>
#include <ibtk/StreamableManager.h>

#include <tbox/AbstractStream.h>
#include <tbox/Pointer.h>

#include <vector>

#include <ibamr/namespaces.h> // IWYU pragma: keep

namespace IBAMR
{
IBAngularSpringForceSpec::Factory::Factory()
{
    setStreamableClassID(StreamableManager::getUnregisteredID());
}

int
IBAngularSpringForceSpec::Factory::getStreamableClassID() const
{
    return STREAMABLE_CLASS_ID;
}

void
IBAngularSpringForceSpec::Factory::setStreamableClassID(const int class_id)
{
    STREAMABLE_CLASS_ID = class_id;
}

Pointer<Streamable>
IBAngularSpringForceSpec::Factory::unpackStream(AbstractStream& stream, const IntVector<NDIM>& /*offset*/)
{
    int n = 0;
    stream.unpack(&n, 1);
    Pointer<IBAngularSpringForceSpec> result = new IBAngularSpringForceSpec(n);
    stream.unpack(&result->d_master_idx, 1);
    if (n > 0)
    {
        std::vector<int> neighbors(2 * n);
        stream.unpack(neighbors.data(), 2 * n);
        for (int k = 0; k < n; ++k)
        {
            result->d_neighbor_idxs[k] = { neighbors[2 * k], neighbors[2 * k + 1] };
        }
        stream.unpack(result->d_stiffnesses.data(), n);
        stream.unpack(result->d_rest_angles.data(), n);
    }
    return result;
}
} // namespace IBAMR
