// ---------------------------------------------------------------------
//
// Copyright (c) 2011 - 2019 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#ifndef included_IBTK_KrylovLinearSolver
#define included_IBTK_KrylovLinearSolver

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "ibtk/LinearOperator.h"
#include "ibtk/LinearSolver.h"

#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/solv/SAMRAIVectorReal.h"


namespace IBTK
{
class HierarchyMathOps;
} // namespace IBTK

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class KrylovLinearSolver provides an abstract interface for the
 * implementation of Krylov subspace solvers for linear problems of the form
 * \f$Ax=b\f$.
 */
class KrylovLinearSolver : public LinearSolver
{
public:
    /*!
     * \brief Default constructor.
     */
    KrylovLinearSolver() = default;

    /*!
     * \brief Empty destructor.
     */
    ~KrylovLinearSolver() = default;

    /*!
     * \brief Set the HierarchyMathOps object used by the solver.
     */
    void setHierarchyMathOps(std::shared_ptr<HierarchyMathOps> hier_math_ops) override;

    /*!
     * \name General-purpose solver functionality.
     */
    //\{

    /*!
     * \brief Set whether the solver should use homogeneous boundary conditions.
     */
    void setHomogeneousBc(bool homogeneous_bc) override;

    /*!
     * \brief Set the time at which the solution is to be evaluated.
     */
    void setSolutionTime(double solution_time) override;

    /*!
     * \brief Set the current time interval.
     */
    void setTimeInterval(double current_time, double new_time) override;

    //\}

    /*!
     * \name Krylov solver functionality.
     */
    //\{

    /*!
     * \brief Set the linear operator used when solving \f$Ax=b\f$.
     */
    virtual void setOperator(std::shared_ptr<LinearOperator> A);

    /*!
     * \brief Retrieve the linear operator used when solving \f$Ax=b\f$.
     */
    virtual std::shared_ptr<LinearOperator> getOperator() const;

    /*!
     * \brief Set the preconditioner used by the Krylov subspace method when
     * solving \f$Ax=b\f$.
     *
     * \note If the preconditioner is NULL, no preconditioning is performed.
     */
    virtual void setPreconditioner(std::shared_ptr<LinearSolver> pc_solver = NULL);

    /*!
     * \brief Retrieve the preconditioner used by the Krylov subspace method
     * when solving \f$Ax=b\f$.
     */
    virtual std::shared_ptr<LinearSolver> getPreconditioner() const;

    //\}

protected:
    // Solver components.
    std::shared_ptr<LinearOperator> d_A;
    std::shared_ptr<LinearSolver> d_pc_solver;
    std::shared_ptr<SAMRAI::solv::SAMRAIVectorReal<double> > d_x, d_b;

private:
    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    KrylovLinearSolver(const KrylovLinearSolver& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    KrylovLinearSolver& operator=(const KrylovLinearSolver& that) = delete;
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBTK_KrylovLinearSolver
