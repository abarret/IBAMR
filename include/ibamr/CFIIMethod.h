// ---------------------------------------------------------------------
//
// Copyright (c) 2018 - 2021 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#ifndef included_IBAMR_CFIIMethod
#define included_IBAMR_CFIIMethod

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "ibamr/IBStrategy.h"
#include <ibamr/IIMethod.h>

#include "ibtk/FEDataManager.h"
#include "ibtk/ibtk_utilities.h"
#include "ibtk/libmesh_utilities.h"

#include "GriddingAlgorithm.h"
#include "IntVector.h"
#include "LoadBalancer.h"
#include "PatchHierarchy.h"
#include "SideIndex.h"
#include "tbox/Pointer.h"

#include "libmesh/enum_fe_family.h"
#include "libmesh/enum_order.h"
#include "libmesh/enum_quadrature_type.h"
#include "libmesh/vector_value.h"

#include <limits>
#include <memory>
#include <set>
#include <string>
#include <vector>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class CFIIMethod is an implementation of the abstract base
 * class IBStrategy that provides functionality required by the IB method with
 * a finite element representation of a surface mesh.
 *
 * Coupling schemes include both IB formulations (integral operations with
 * regularized delta function kernels) and an immersed interface method (IIM)
 * scheme (E. M. Kolahdouz, A. P. S. Bhalla, B. A. Craven, and B. E. Griffith.
 * An immersed interface method for discrete surfaces. J Comput Phys,
 * 400:108854 (37 pages), 2020).
 *
 * \note When using the IIM implementation, it is recommended that users set
 * all linear solvers to use tight relative tolerances (1e-10).
 */
class CFIIMethod : public IIMethod
{
public:
    static const std::string EXTRA_STRESS_JUMP_SYSTEM_NAME;
    static const std::string EXTRA_STRESS_IN_SYSTEM_NAME;
    static const std::string EXTRA_STRESS_OUT_SYSTEM_NAME;

    /*!
     * \brief Constructor.
     */
    CFIIMethod(const std::string& object_name,
               SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
               libMesh::MeshBase* mesh,
               int max_level_number,
               bool register_for_restart = true,
               const std::string& restart_read_dirname = "",
               unsigned int restart_restore_number = 0);

    /*!
     * \brief Constructor.
     */
    CFIIMethod(const std::string& object_name,
               SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
               const std::vector<libMesh::MeshBase*>& meshes,
               int max_level_number,
               bool register_for_restart = true,
               const std::string& restart_read_dirname = "",
               unsigned int restart_restore_number = 0);

    /*!
     * \brief Destructor.
     */
    virtual ~CFIIMethod() = default;

    /*!
     * Compute the jump in the extra stress tensor
     */
    virtual void computeStressJumps(int sig_idx, const int ls_idx, double data_time, unsigned int part = 0);

    /*!
     * Compute the Lagrangian force at the specified time within the current
     * time interval.
     */
    void computeLagrangianForce(double data_time) override;

    void
    spreadForce(int f_data_idx,
                IBTK::RobinPhysBdryPatchStrategy* f_phys_bdry_op,
                const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > >& f_prolongation_scheds,
                double data_time) override;

    /*!
     * Initialize the FE equation systems objects.  This method must be called
     * prior to calling initializeFEData().
     */
    void initializeFEEquationSystems() override;

    /*!
     * Initialize FE data.  This method must be called prior to calling
     * IBHierarchyIntegrator::initializePatchHierarchy().
     */
    void initializeFEData() override;

    /*!
     * \brief Register Eulerian variables with the parent IBHierarchyIntegrator.
     */
    void registerEulerianVariables() override;

protected:
    /*!
     * Impose the jump conditions.
     */
    void imposeJumpConditions(const int f_data_idx,
                              libMesh::PetscVector<double>& P_jump_ghost_vec,
                              std::array<libMesh::PetscVector<double>*, NDIM>& DU_jump_ghost_vec,
                              libMesh::PetscVector<double>& X_ghost_vec,
                              const double data_time,
                              const unsigned int part) override;

    std::vector<std::array<libMesh::System*, NDIM*(NDIM + 1) / 2> > d_sig_jump_systems, d_sig_in_systems,
        d_sig_out_systems;

    /*
     * Method parameters.
     */
    bool d_use_extra_stress_jump_conditions = false;
    libMesh::FEFamily d_sig_jump_fe_family = libMesh::LAGRANGE;
    double d_sig_calc_width = 0.5;
    bool d_use_bilinear_interp = false;

    /*
     * Eulerian data.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > d_sig_var;
    int d_sig_scratch_idx = IBTK::invalid_index;

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    CFIIMethod() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    CFIIMethod(const CFIIMethod& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    CFIIMethod& operator=(const CFIIMethod& that) = delete;

    /*!
     * Implementation of class constructor.
     */
    void commonConstructor(const std::string& object_name,
                           SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                           const std::vector<libMesh::MeshBase*>& meshes,
                           int max_level_number,
                           bool register_for_restart,
                           const std::string& restart_read_dirname,
                           unsigned int restart_restore_number);

    /*!
     * Read input values from a given database.
     */
    void getFromInput(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db, bool is_from_restart);

    /*!
     * Read object state from the restart file and initialize class data
     * members.
     */
    void getFromRestart();

    /*!
     * Given a location and a candidate list of stencils, compute weights for a polyharmonic spline appended with linear
     * polynomials that are in the same order as those in xpts. The optional scale argument scales the linear
     * polynomials. The default value is 1.0, but the suggested value scales with the grid size.
     */
    template <typename Vector>
    std::vector<double>
    findLinearInterpWeights(const Vector& eval_pt, const std::vector<Vector>& xpts, const double scale = 1.0);

    /*!
     * Find a stencil consisting of stencil_size number of points and store the physical locations in xpts and the cell
     * indices in idxs. The stencil is chosen so that it contains points that are on the same normal as the vector in
     * normal. A flooding algorithm is used to find the stencil.
     */
    template <typename Vector>
    void findInterpPts(const SAMRAI::pdat::CellIndex<NDIM>& eval_idx,
                       int ls_idx,
                       SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                       int stencil_size,
                       std::vector<Vector>& xpts,
                       std::vector<SAMRAI::pdat::CellIndex<NDIM> >& idxs);

    /*!
     * Evaluate the stencil given the weights and point locations. This computes sum(wgts * data(idxs, depth)). The
     * default value of depth is 0.
     */
    double evaluateInterpolant(const SAMRAI::pdat::CellData<NDIM, double>& data,
                               const std::vector<SAMRAI::pdat::CellIndex<NDIM> >& idxs,
                               const std::vector<double>& wgts,
                               int depth = 0);

    int d_stencil_size = 10;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBAMR_CFIIMethod
