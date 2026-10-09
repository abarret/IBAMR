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

#include <ibtk/IBTK_MPI.h>
#include <ibtk/StreamableFactory.h>
#include <ibtk/StreamableManager.h>

#include <ibamr/namespaces.h> // IWYU pragma: keep

namespace IBAMR
{
int IBAngularSpringForceSpec::STREAMABLE_CLASS_ID = StreamableManager::getUnregisteredID();

void
IBAngularSpringForceSpec::registerWithStreamableManager()
{
    IBTK_MPI::barrier();
    if (!getIsRegisteredWithStreamableManager())
    {
#if !defined(NDEBUG)
        TBOX_ASSERT(STREAMABLE_CLASS_ID == StreamableManager::getUnregisteredID());
#endif
        STREAMABLE_CLASS_ID = StreamableManager::getManager()->registerFactory(new IBAngularSpringForceSpecFactory());
    }
    IBTK_MPI::barrier();
}
} // namespace IBAMR
