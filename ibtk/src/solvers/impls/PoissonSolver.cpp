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

#include "ibtk/PoissonSolver.h"
#include "ibtk/namespaces.h" // IWYU pragma: keep

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/solv/LocationIndexRobinBcCoefs.h"
#include "SAMRAI/solv/PoissonSpecifications.h"
#include "SAMRAI/solv/RobinBcCoefStrategy.h"
#include "SAMRAI/tbox/Database.h"


#include <string>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

void
PoissonSolver::setPoissonSpecifications(const PoissonSpecifications& poisson_spec)
{
    d_poisson_spec = poisson_spec;
    return;
} // setPoissonSpecifications

void
PoissonSolver::setPhysicalBcCoef(RobinBcCoefStrategy* const bc_coef)
{
    setPhysicalBcCoefs(std::vector<RobinBcCoefStrategy*>(1, bc_coef));
    return;
} // setPhysicalBcCoef

void
PoissonSolver::setPhysicalBcCoefs(const std::vector<RobinBcCoefStrategy*>& bc_coefs)
{
    d_bc_coefs.resize(bc_coefs.size());
    for (unsigned int l = 0; l < bc_coefs.size(); ++l)
    {
        if (bc_coefs[l])
        {
            d_bc_coefs[l] = bc_coefs[l];
        }
        else
        {
            d_bc_coefs[l] = d_default_bc_coef.get();
        }
    }
    return;
} // setPhysicalBcCoefs

void
PoissonSolver::initSpecialized(const std::string& object_name, const bool /*homogeneous_bc*/)
{
    // Initialize the Poisson specifications.
    PoissonSpecifications poisson_spec(object_name + "::poisson_spec");
    poisson_spec.setCZero();
    poisson_spec.setDConstant(-1.0);
    setPoissonSpecifications(poisson_spec);

    // Initialize the boundary conditions.
    d_default_bc_coef.reset(
        std::make_shared<LocationIndexRobinBcCoefs>(object_name + "::default_bc_coef", std::shared_ptr<Database>(nullptr)));
    auto p_default_bc_coef = dynamic_cast<LocationIndexRobinBcCoefs*>(d_default_bc_coef.get());
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        p_default_bc_coef->setBoundaryValue(2 * d, 0.0);
        p_default_bc_coef->setBoundaryValue(2 * d + 1, 0.0);
    }
    setPhysicalBcCoef(d_default_bc_coef.get());
} // initSpecialized

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
