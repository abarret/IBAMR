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

#include "ibamr/config.h"

#include "ibamr/cut_cells/LSCartGridFunction.h"

#include <SAMRAI_config.h>

#include "ibamr/app_namespaces.h"

namespace LS
{
/////////////////////////////// PUBLIC ///////////////////////////////////////

LSCartGridFunction::LSCartGridFunction(const string& object_name) : CartGridFunction(object_name)
{
    // intentionally blank
    return;
} // LSCartGridFunction

} // namespace LS