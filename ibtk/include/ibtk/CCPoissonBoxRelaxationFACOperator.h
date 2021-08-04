// ---------------------------------------------------------------------
//
// Copyright (c) 2015 - 2019 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#ifndef included_IBTK_CCPoissonBoxRelaxationFACOperator
#define included_IBTK_CCPoissonBoxRelaxationFACOperator

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "ibtk/PoissonFACPreconditioner.h"
#include "ibtk/PoissonFACPreconditionerStrategy.h"
#include "ibtk/PoissonSolver.h"

#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/solv/PoissonSpecifications.h"
#include "SAMRAI/tbox/Database.h"


#include "petscksp.h"
#include "petscmat.h"
#include "petscvec.h"

#include <map>
#include <string>
#include <vector>

namespace SAMRAI
{
namespace hier
{

class Box;

class BoxList;

class Patch;
} // namespace hier
namespace pdat
{
template <class TYPE>
class CellData;
template <class TYPE>
class SideData;
} // namespace pdat
namespace solv
{
template <class TYPE>
class SAMRAIVectorReal;
} // namespace solv
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class CCPoissonBoxRelaxationFACOperator is a concrete
 * PoissonFACPreconditionerStrategy for solving elliptic equations of the form
 * \f$ \mbox{$L u$} = \mbox{$(C I + \nabla \cdot D \nabla) u$} = f \f$ using a
 * globally second-order accurate cell-centered finite-volume discretization,
 * with support for Robin and periodic boundary conditions.
 *
 * This class provides operators that are used by class FACPreconditioner to
 * solve scalar Poisson-type equations of the form \f[ (C I + \nabla \cdot D
 * \nabla) u = f \f] using a cell-centered, globally second-order accurate
 * finite-volume discretization, where
 *
 * - \f$ C \f$, \f$ D \f$ and \f$ f \f$ are independent of \f$ u \f$,
 * - \f$ C \f$ is a cell-centered scalar field,
 * - \f$ D \f$ is a side-centered scalar field of diffusion coefficients, and
 * - \f$ f \f$ is a cell-centered scalar function.
 *
 * Robin boundary conditions may be specified at physical boundaries; see class
 * SAMRAI::solv::RobinBcCoefStrategy.
 *
 * By default, the class is configured to solve the Poisson problem \f$
 * -\nabla^2 u = f \f$, subject to homogeneous Dirichlet boundary conditions.
 *
 * Sample parameters for initialization from database (and their default
 * values): \verbatim

 smoother_type = "PATCH_GAUSS_SEIDEL"         // see setSmootherType()
 prolongation_method = "LINEAR_REFINE"        // see setProlongationMethod()
 restriction_method = "CONSERVATIVE_COARSEN"  // see setRestrictionMethod()
 coarse_solver_type = "HYPRE_LEVEL_SOLVER"    // see setCoarseSolverType()
 coarse_solver_rel_residual_tol = 1.0e-5      // see setCoarseSolverRelativeTolerance()
 coarse_solver_abs_residual_tol = 1.0e-50     // see setCoarseSolverAbsoluteTolerance()
 coarse_solver_max_iterations = 1             // see setCoarseSolverMaxIterations()
 coarse_solver_db {                           // SAMRAI::tbox::Database for initializing coarse
 level solver
    solver_type = "PFMG"
    num_pre_relax_steps = 0
    num_post_relax_steps = 2
 }
 \endverbatim
*/
class CCPoissonBoxRelaxationFACOperator : public PoissonFACPreconditionerStrategy
{
public:
    /*!
     * \brief Constructor.
     */
    CCPoissonBoxRelaxationFACOperator(const std::string& object_name,
                                      std::shared_ptr<SAMRAI::tbox::Database> input_db,
                                      std::string default_options_prefix);

    /*!
     * \brief Destructor.
     */
    ~CCPoissonBoxRelaxationFACOperator();

    /*!
     * \brief Static function to construct a PoissonFACPreconditioner with a
     * CCPoissonBoxRelaxationFACOperator FAC strategy.
     */
    static std::shared_ptr<PoissonSolver> allocate_solver(const std::string& object_name,
                                                                std::shared_ptr<SAMRAI::tbox::Database> input_db,
                                                                const std::string& default_options_prefix)
    {
        std::shared_ptr<PoissonFACPreconditionerStrategy> fac_operator = std::make_shared<CCPoissonBoxRelaxationFACOperator>(
            object_name + "::CCPoissonBoxRelaxationFACOperator", input_db, default_options_prefix);
        return std::make_shared<PoissonFACPreconditioner>(object_name, fac_operator, input_db, default_options_prefix);
    } // allocate

    /*!
     * \name Functions for configuring the solver.
     */
    //\{

    /*!
     * \brief Specify the smoother type.
     *
     * Select from:
     * - \c "PATCH_GAUSS_SEIDEL"
     * - \c "PROCESSOR_GAUSS_SEIDEL"
     */
    void setSmootherType(const std::string& smoother_type) override;

    /*!
     * \brief Specify the coarse level solver.
     */
    void setCoarseSolverType(const std::string& coarse_solver_type) override;

    //\}

    /*!
     * \name Implementation of FACPreconditionerStrategy interface.
     */
    //\{

    /*!
     * \brief Perform a given number of relaxations on the error.
     *
     * \param error error vector
     * \param residual residual vector
     * \param level_num level number
     * \param num_sweeps number of sweeps to perform
     * \param performing_pre_sweeps boolean value that is true when pre-smoothing sweeps are
     *being
     *performed
     * \param performing_post_sweeps boolean value that is true when post-smoothing sweeps are
     *being
     *performed
     */
    void smoothError(SAMRAI::solv::SAMRAIVectorReal<double>& error,
                     const SAMRAI::solv::SAMRAIVectorReal<double>& residual,
                     int level_num,
                     int num_sweeps,
                     bool performing_pre_sweeps,
                     bool performing_post_sweeps) override;

    /*!
     * \brief Solve the residual equation Ae=r on the coarsest level of the
     * patch hierarchy.
     *
     * \param error error vector
     * \param residual residual vector
     * \param coarsest_ln coarsest level number
     */
    bool solveCoarsestLevel(SAMRAI::solv::SAMRAIVectorReal<double>& error,
                            const SAMRAI::solv::SAMRAIVectorReal<double>& residual,
                            int coarsest_ln) override;

    /*!
     * \brief Compute composite grid residual on a range of levels.
     *
     * \param residual residual vector
     * \param solution solution vector
     * \param rhs source (right hand side) vector
     * \param coarsest_level_num coarsest level number
     * \param finest_level_num finest level number
     */
    void computeResidual(SAMRAI::solv::SAMRAIVectorReal<double>& residual,
                         const SAMRAI::solv::SAMRAIVectorReal<double>& solution,
                         const SAMRAI::solv::SAMRAIVectorReal<double>& rhs,
                         int coarsest_level_num,
                         int finest_level_num) override;

    //\}

protected:
    /*!
     * \brief Compute implementation-specific hierarchy-dependent data.
     */
    void initializeOperatorStateSpecialized(const SAMRAI::solv::SAMRAIVectorReal<double>& solution,
                                            const SAMRAI::solv::SAMRAIVectorReal<double>& rhs,
                                            int coarsest_reset_ln,
                                            int finest_reset_ln) override;

    /*!
     * \brief Remove implementation-specific hierarchy-dependent data.
     */
    void deallocateOperatorStateSpecialized(int coarsest_reset_ln, int finest_reset_ln) override;

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    CCPoissonBoxRelaxationFACOperator() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    CCPoissonBoxRelaxationFACOperator(const CCPoissonBoxRelaxationFACOperator& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    CCPoissonBoxRelaxationFACOperator& operator=(const CCPoissonBoxRelaxationFACOperator& that) = delete;

    /*!
     * \brief Construct a matrix corresponding to a Laplace operator restricted
     * to a single patch.
     */
    static void buildPatchLaplaceOperator(Mat& A,
                                          const SAMRAI::solv::PoissonSpecifications& poisson_spec,
                                          std::shared_ptr<SAMRAI::hier::Patch > patch,
                                          const SAMRAI::hier::IntVector& ghost_cell_width);

    /*!
     * \brief Construct a matrix corresponding to a Laplace operator restricted
     * to a single patch with grid aligned anisotropy.
     */
    static void buildPatchLaplaceOperator_aligned(Mat& A,
                                                  std::shared_ptr<SAMRAI::pdat::CellData<double> > C_data,
                                                  std::shared_ptr<SAMRAI::pdat::SideData<double> > D_data,
                                                  std::shared_ptr<SAMRAI::hier::Patch > patch,
                                                  const SAMRAI::hier::IntVector& ghost_cell_width);

    /*!
     * \brief Construct a matrix corresponding to a Laplace operator restricted
     * to a single patch with non-grid aligned anisotropy.
     */
    static void
    buildPatchLaplaceOperator_nonaligned(Mat& A,
                                         std::shared_ptr<SAMRAI::pdat::CellData<double> > C_data,
                                         std::shared_ptr<SAMRAI::pdat::SideData<double> > D_data,
                                         std::shared_ptr<SAMRAI::hier::Patch > patch,
                                         const SAMRAI::hier::IntVector& ghost_cell_width);

    /*
     * Coarse level solvers and solver parameters.
     */
    std::shared_ptr<PoissonSolver> d_coarse_solver;
    std::shared_ptr<SAMRAI::tbox::Database> d_coarse_solver_db;

    /*
     * Mappings from patch indices to patch operators.
     */
    std::string d_petsc_options_prefix = "cc_poisson_fac_";
    std::vector<std::vector<Vec> > d_patch_vec_e, d_patch_vec_f;
    std::vector<std::vector<Mat> > d_patch_mat;
    std::vector<std::vector<KSP> > d_patch_ksp;

    /*
     * Patch overlap data.
     */
    std::vector<std::vector<SAMRAI::hier::BoxList > > d_patch_bc_box_overlap;
    std::vector<std::vector<std::map<int, SAMRAI::hier::Box > > > d_patch_neighbor_overlap;
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBTK_CCPoissonBoxRelaxationFACOperator
