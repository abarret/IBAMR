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

#include "ibtk/IBTK_CHKERRQ.h"
#include "ibtk/IBTK_MPI.h"
#include "ibtk/JacobianOperator.h"
#include "ibtk/PETScSAMRAIVectorReal.h"
#include "ibtk/PETScSNESJacobianJOWrapper.h"
#include "ibtk/namespaces.h" // IWYU pragma: keep

#include "SAMRAI/hier/IntVector.h"
#include "MultiblockDataTranslator.h"
#include "SAMRAI/solv/SAMRAIVectorReal.h"


#include "petscsnes.h"
#include "petscvec.h"
#include <petsclog.h>

#include <mpi.h>

#include <string>
#include <utility>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

PETScSNESJacobianJOWrapper::PETScSNESJacobianJOWrapper(
    std::string object_name,
    SNES petsc_snes,
    PetscErrorCode (*const petsc_snes_form_jac)(SNES, Vec, Mat, Mat, void*),
    void* const petsc_snes_jac_ctx)
    : JacobianOperator(std::move(object_name)),
      d_petsc_snes(std::move(petsc_snes)),
      d_petsc_snes_form_jac(petsc_snes_form_jac),
      d_petsc_snes_jac_ctx(petsc_snes_jac_ctx)
{
    // intentionally blank
}

PETScSNESJacobianJOWrapper::~PETScSNESJacobianJOWrapper()
{
    if (d_is_initialized) deallocateOperatorState();
}

const SNES&
PETScSNESJacobianJOWrapper::getPETScSNES() const
{
    return d_petsc_snes;
}

PetscErrorCode (*PETScSNESJacobianJOWrapper::getPETScSNESFormJacobian())(SNES, Vec, Mat, Mat, void*)
{
    return d_petsc_snes_form_jac;
}

void*
PETScSNESJacobianJOWrapper::getPETScSNESJacobianContext() const
{
    return d_petsc_snes_jac_ctx;
}

void
PETScSNESJacobianJOWrapper::formJacobian(SAMRAIVectorReal<double>& x)
{
    // Create the PETSc Vec wrappers.
    Vec petsc_x = PETScSAMRAIVectorReal::createPETScVector(std::shared_ptr<SAMRAIVectorReal<double> >(&x, false));

    // Setup the Jacobian matrix.
    int ierr = d_petsc_snes_form_jac(d_petsc_snes, petsc_x, d_petsc_snes_jac, nullptr, d_petsc_snes_jac_ctx);
    IBTK_CHKERRQ(ierr);

    // Destroy the PETSc Vec wrappers.
    PETScSAMRAIVectorReal::destroyPETScVector(petsc_x);
    petsc_x = nullptr;
}

std::shared_ptr<SAMRAIVectorReal<double> >
PETScSNESJacobianJOWrapper::getBaseVector() const
{
    Vec petsc_x;
    int ierr = SNESGetSolution(d_petsc_snes, &petsc_x);
    IBTK_CHKERRQ(ierr);
    std::shared_ptr<SAMRAIVectorReal<double> > samrai_x;
    PETScSAMRAIVectorReal::getSAMRAIVectorRead(petsc_x, &samrai_x);
    std::shared_ptr<SAMRAIVectorReal<double> > samrai_x_ptr = samrai_x;
    PETScSAMRAIVectorReal::restoreSAMRAIVectorRead(petsc_x, &samrai_x);
    return samrai_x_ptr;
}

void
PETScSNESJacobianJOWrapper::apply(SAMRAIVectorReal<double>& x, SAMRAIVectorReal<NDIM, double>& y)
{
    if (!d_is_initialized) initializeOperatorState(x, y);

    // Update the PETSc Vec wrappers.
    PETScSAMRAIVectorReal::replaceSAMRAIVector(d_petsc_x, std::shared_ptr<SAMRAIVectorReal<double> >(&x, false));
    PETScSAMRAIVectorReal::replaceSAMRAIVector(d_petsc_y, std::shared_ptr<SAMRAIVectorReal<double> >(&y, false));

    // Apply the operator.
    int ierr = MatMult(d_petsc_snes_jac, d_petsc_x, d_petsc_y);
    IBTK_CHKERRQ(ierr);
    return;
}

void
PETScSNESJacobianJOWrapper::applyAdd(SAMRAIVectorReal<double>& x,
                                     SAMRAIVectorReal<double>& y,
                                     SAMRAIVectorReal<double>& z)
{
    if (!d_is_initialized) initializeOperatorState(x, y);

    // Update the PETSc Vec wrappers.
    PETScSAMRAIVectorReal::replaceSAMRAIVector(d_petsc_x, std::shared_ptr<SAMRAIVectorReal<double> >(&x, false));
    PETScSAMRAIVectorReal::replaceSAMRAIVector(d_petsc_y, std::shared_ptr<SAMRAIVectorReal<double> >(&y, false));
    PETScSAMRAIVectorReal::replaceSAMRAIVector(d_petsc_z, std::shared_ptr<SAMRAIVectorReal<double> >(&z, false));

    // Apply the operator.
    int ierr = MatMultAdd(d_petsc_snes_jac, d_petsc_x, d_petsc_y, d_petsc_z);
    IBTK_CHKERRQ(ierr);
}

void
PETScSNESJacobianJOWrapper::initializeOperatorState(const SAMRAIVectorReal<double>& in,
                                                    const SAMRAIVectorReal<double>& out)
{
    if (d_is_initialized) deallocateOperatorState();
    d_x = in.cloneVector("");
    d_y = out.cloneVector("");
    d_z = out.cloneVector("");
    MPI_Comm comm;
    int ierr = PetscObjectGetComm(reinterpret_cast<PetscObject>(d_petsc_snes_jac), &comm);
    IBTK_CHKERRQ(ierr);
    d_petsc_x = PETScSAMRAIVectorReal::createPETScVector(d_x, comm);
    d_petsc_y = PETScSAMRAIVectorReal::createPETScVector(d_y, comm);
    d_petsc_z = PETScSAMRAIVectorReal::createPETScVector(d_z, comm);
    d_is_initialized = true;
}

void
PETScSNESJacobianJOWrapper::deallocateOperatorState()
{
    if (!d_is_initialized) return;
    PETScSAMRAIVectorReal::destroyPETScVector(d_petsc_x);
    PETScSAMRAIVectorReal::destroyPETScVector(d_petsc_y);
    PETScSAMRAIVectorReal::destroyPETScVector(d_petsc_z);
    d_x->freeVectorComponents();
    d_x.setNull();
    d_y->freeVectorComponents();
    d_y.setNull();
    d_z->freeVectorComponents();
    d_z.setNull();
    d_is_initialized = false;
}

//////////////////////////////////////////////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
