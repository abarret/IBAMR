// ---------------------------------------------------------------------
//
// Copyright (c) 2015 - 2023 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

/////////////////////////////// INCLUDE GUARD ////////////////////////////////

#ifndef included_IBAMR_CIBMethod
#define included_IBAMR_CIBMethod

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include "ibamr/CIBStrategy.h"
#include "ibamr/IBMethod.h"
#include "ibamr/ibamr_enums.h"

#include "ibtk/LData.h"
#include "ibtk/LDataManager.h"
#include "ibtk/ibtk_utilities.h"

#include "RobinBcCoefStrategy.h"
#include "Variable.h"
#include "VisItDataWriter.h"
#include "tbox/Pointer.h"

#include "petscmat.h"
#include "petscvec.h"

#include <iosfwd>
#include <string>
#include <utility>
#include <vector>

namespace IBTK
{
class LData;
class RobinPhysBdryPatchStrategy;
} // namespace IBTK
namespace SAMRAI
{
namespace hier
{
template <int DIM>
class BasePatchHierarchy;
template <int DIM>
class BasePatchLevel;
template <int DIM>
class PatchHierarchy;
} // namespace hier
namespace mesh
{
template <int DIM>
class GriddingAlgorithm;
} // namespace mesh
namespace tbox
{
class Database;
} // namespace tbox
namespace xfer
{
template <int DIM>
class CoarsenSchedule;
template <int DIM>
class RefineSchedule;
} // namespace xfer
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class CIBMethod is a concrete CIBStrategy and IBMethod
 * class which implements the motion of rigid bodies using the constraint
 * formulation. The immersed structure is discretized using standard IB
 * markers.
 */

class CIBMethod : public IBAMR::IBMethod, public IBAMR::CIBStrategy
{
public:
    /*!
     * Explicitly use the other overload of this function from the base class.
     */
    using CIBStrategy::setRigidBodyVelocity;

    /*!
     * Explicitly use the other overloads of this function from the base class.
     */
    using CIBStrategy::computeNetRigidGeneralizedForce;

    /*!
     * \brief Constructor of the class.
     */
    CIBMethod(std::string object_name,
              SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
              const int no_structures = 1,
              bool register_for_restart = true);

    /*!
     * \brief Destructor of the class.
     */
    ~CIBMethod();

    /*!
     * \brief Typedef specifying interface for specifying constrained body velocities.
     */
    using ConstrainedNodalVelocityFcnPtr = void (*)(Vec U_k,
                                                    const IBTK::RigidDOFVector& U,
                                                    Vec X,
                                                    const Eigen::Vector3d& X_com,
                                                    const Eigen::Matrix3d& rotation_mat,
                                                    double data_time,
                                                    void* ctx);

    using ConstrainedCOMVelocityFcnPtr = void (*)(double data_time,
                                                  Eigen::Vector3d& U_com,
                                                  Eigen::Vector3d& W_com,
                                                  void* ctx);

    /*!
     * \brief Typedef specifying interface for specifying net external force
     * and torque on structures.
     */
    using ExternalForceTorqueFcnPtr = void (*)(double data_time, Eigen::Vector3d& F, Eigen::Vector3d& T, void* ctx);

    /*!
     * \brief Callbacks before INS is integrated.
     */
    using preprocessSolveFluidEqn_callbackfcn = void (*)(const double, const double, const int, void*);

    /*!
     * \brief Struct encapsulating constrained velocity functions data.
     */
    struct ConstrainedVelocityFcnsData
    {
        ConstrainedVelocityFcnsData(ConstrainedNodalVelocityFcnPtr nodalvelfcn = nullptr,
                                    ConstrainedCOMVelocityFcnPtr comvelfcn = nullptr,
                                    void* ctx = nullptr)
            : nodalvelfcn(nodalvelfcn), comvelfcn(comvelfcn), ctx(ctx)
        {
            // intentionally left blank
        }
        ConstrainedNodalVelocityFcnPtr nodalvelfcn;
        ConstrainedCOMVelocityFcnPtr comvelfcn;
        void* ctx;
    };

    /*!
     * \brief Struct encapsulating external force and torque function data.
     */
    struct ExternalForceTorqueFcnData
    {
        ExternalForceTorqueFcnData(ExternalForceTorqueFcnPtr forcetorquefcn = nullptr, void* ctx = nullptr)
            : forcetorquefcn(forcetorquefcn), ctx(ctx)
        {
            // intentionally left blank
        }
        ExternalForceTorqueFcnPtr forcetorquefcn;
        void* ctx;
    };

    /*!
     * \brief Register user defined constrained velocity functions.
     */
    void registerConstrainedVelocityFunction(ConstrainedNodalVelocityFcnPtr nodalvelfcn,
                                             ConstrainedCOMVelocityFcnPtr comvelfcn,
                                             void* ctx = nullptr,
                                             unsigned int part = 0);

    /*!
     * \brief Register user defined constrained velocity function data.
     */
    void registerConstrainedVelocityFunction(const ConstrainedVelocityFcnsData& data, unsigned int part = 0);

    /*!
     * \brief Register an external force and torque function.
     */
    void registerExternalForceTorqueFunction(ExternalForceTorqueFcnPtr forcetorquefcn,
                                             void* ctx = nullptr,
                                             unsigned int part = 0);
    /*!
     * \brief Register an external force and torque function data.
     */
    void registerExternalForceTorqueFunction(const ExternalForceTorqueFcnData& data, unsigned int part = 0);

    /*!
     * \brief Get the level on which the structures reside.
     */
    int getStructuresLevelNumber() const;

    /*!
     * \brief Get the structure handle to which this Lagrangian index belongs.
     */
    int getStructureHandle(const int lag_idx) const;

    /*!
     * \brief Register any preprocess fluid solve callback functions.
     */
    void registerPreProcessSolveFluidEquationsCallBackFcn(preprocessSolveFluidEqn_callbackfcn callback, void* ctx);

    /*!
     * \brief Calculate any body forces for INS solver over here.
     */
    void preprocessSolveFluidEquations(double current_time, double new_time, int cycle_num) override;

    /*!
     * \brief Register Eulerian variables with the parent IBHierarchyIntegrator.
     */
    void registerEulerianVariables() override;

    /*!
     * \brief Register Eulerian refinement or coarsening algorithms with the parent
     * IBHierarchyIntegrator.
     */
    void registerEulerianCommunicationAlgorithms() override;

    /*!
     * \brief Method to prepare to advance data from current_time to new_time.
     */
    void preprocessIntegrateData(double current_time, double new_time, int num_cycles) override;

    /*!
     * \brief Method to clean up data following call(s) to integrateHierarchy().
     */
    void postprocessIntegrateData(double current_time, double new_time, int num_cycles) override;

    /*!
     * \brief Initialize data on a new level after it is inserted into an AMR patch
     * hierarchy by the gridding algorithm.
     *
     * \see SAMRAI::mesh::StandardTagAndInitStrategy::initializeLevelData
     */
    void initializeLevelData(SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
                             int level_number,
                             double init_data_time,
                             bool can_be_refined,
                             bool initial_time,
                             SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchLevel<NDIM> > old_level,
                             bool allocate_data) override;

    /*!
     * \brief Initialize Lagrangian data corresponding to the given AMR patch hierarchy
     * at the start of a computation.  If the computation is begun from a
     * restart file, data may be read from the restart databases.
     *
     * A patch data descriptor is provided for the Eulerian velocity in case
     * initialization requires interpolating Eulerian data.  Ghost cells for
     * Eulerian data will be filled upon entry to this function.
     */
    void initializePatchHierarchy(
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
        SAMRAI::tbox::Pointer<SAMRAI::mesh::GriddingAlgorithm<NDIM> > gridding_alg,
        int u_data_idx,
        const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenSchedule<NDIM> > >& u_synch_scheds,
        const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > >& u_ghost_fill_scheds,
        int integrator_step,
        double init_data_time,
        bool initial_time) override;

    /*!
     * \brief Interpolate the Eulerian velocity to the curvilinear mesh at the
     * specified time within the current time interval.
     */
    void interpolateVelocity(
        int u_data_idx,
        const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenSchedule<NDIM> > >& u_synch_scheds,
        const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > >& u_ghost_fill_scheds,
        double data_time) override;

    /*!
     * \brief Spread the Lagrangian force to the Cartesian grid at the specified time
     * within the current time interval.
     */
    void
    spreadForce(int f_data_idx,
                IBTK::RobinPhysBdryPatchStrategy* f_phys_bdry_op,
                const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > >& f_prolongation_scheds,
                double data_time) override;

    /*!
     * \brief Advance the positions of the Lagrangian structure using the forward Euler
     * method.
     */
    void forwardEulerStep(double current_time, double new_time) override;

    /*!
     * \brief Advance the positions of the Lagrangian structure using the backward Euler
     * method.
     */
    void backwardEulerStep(double current_time, double new_time) override;

    /*!
     * \brief Advance the positions of the Lagrangian structure using the (explicit)
     * midpoint rule.
     */
    void midpointStep(double current_time, double new_time) override;

    /*!
     * \brief Advance the positions of the Lagrangian structure using the
     * trapezoidal rule.
     */
    void trapezoidalStep(double current_time, double new_time) override;

    /*!
     * \brief Register VisIt data writer to output data files that
     * may be postprocessed with the VisIt visualization tool.
     */
    void registerVisItDataWriter(SAMRAI::tbox::Pointer<SAMRAI::appu::VisItDataWriter<NDIM> > visit_writer);

    /*!
     * \brief Override the putToDatabase method of the base Serializable class.
     */
    void putToDatabase(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db) override;

    // \{
    // The following are the concrete implementation of CIBStrategy methods:
    //

    // \see CIBStrategy::setConstraintForce() method.
    /*!
     * \brief Set the constraint force in the internal data structures of the
     * class.
     */
    void setConstraintForce(Vec L, const double data_time, const double scale = 1.0) override;

    // \see CIBStrategy::getConstraintForce()
    /*!
     * \brief Get the constraint rigid body force at the specified time within
     * the current time interval.
     */
    void getConstraintForce(Vec* L, const double data_time) override;

    // \see CIBStrategy::getFreeRigidVelocities()
    /*!
     * \brief Get the free rigid velocities (DOFs) at the specified time within
     * the current time interval.
     */
    virtual void getFreeRigidVelocities(Vec* U, const double data_time) override;

    // \see CIBStrategy::getNetExternalForceTorque()
    /*!
     * \brief Get net external force and torque at the specified time within
     * the current time interval.
     */
    virtual void getNetExternalForceTorque(Vec* F, const double data_time) override;

    // \see CIBStrategy::subtractMeanConstraintForce()
    /*!
     * \brief Subtract the mean of constraint force from the background Eulerian
     * grid.
     */
    void subtractMeanConstraintForce(Vec L, int f_data_idx, const double scale = 1.0) override;

    // \see CIBStrategy::setInterpolatedVelocityVector() method
    /*!
     * \brief Prepare the CIBMethod class to get the interpolated fluid
     * velocity.
     */
    void setInterpolatedVelocityVector(Vec V, const double data_time) override;

    // \see CIBStrategy::setInterpolatedVelocityVector() method
    /*!
     * \brief Get interpolated velocity from the Eulerian grid.
     */
    void getInterpolatedVelocity(Vec V, const double data_time, const double scale = 1.0) override;

    // \see CIBStrategy::computeMobilityRegularization method
    /*!
     * \brief Compute regularization vector for the mobility problem.
     *
     */
    void computeMobilityRegularization(Vec D, Vec L, const double scale = 1.0) override;

    // \see CIBStrategy::getNumberOfNodes method
    /*!
     * \brief Get number of nodes for a particular structure registered with CIBMethod.
     */
    unsigned int getNumberOfNodes(const unsigned int part) const override;

    // \see CIBStrategy::setRigidBodyVelocity method
    /*!
     * \brief Set the rigid body velocity at the nodal points
     * contained in the Vec V.
     */
    void setRigidBodyVelocity(const unsigned int part, const IBTK::RigidDOFVector& U, Vec V) override;

    // \see CIBStrategy::computeNetRigidGeneralizedForce() method.
    /*!
     * \brief Compute total force and torque on the rigid structure(s).
     */
    void computeNetRigidGeneralizedForce(const unsigned int part, Vec L, IBTK::RigidDOFVector& F) override;

    // \see CIBStrategy::copyVecToArray() method.
    /*!
     * \brief Copy PETSc Vec to raw array for specified structures.
     */
    void copyVecToArray(Vec b,
                        double* array,
                        const std::vector<unsigned>& struct_ids,
                        const int data_depth,
                        const int array_rank) override;

    // \see CIBStrategy::copyVecToArray() method.
    /*!
     * \brief Copy raw array to PETSc Vec for specified structures.
     */
    void copyArrayToVec(Vec b,
                        double* array,
                        const std::vector<unsigned>& struct_ids,
                        const int data_depth,
                        const int array_rank) override;

    // \see CIBStrategy::constructMobilityMatrix() method.
    /*!
     * \brief Generate dense mobility matrix for the prototypical structures
     * identified by their indices.
     */
    void constructMobilityMatrix(const std::string& mat_name,
                                 MobilityMatrixType mat_type,
                                 Mat& mobility_mat,
                                 const std::vector<unsigned>& prototype_struct_ids,
                                 const double* grid_dx,
                                 const double* domain_extents,
                                 const bool initial_time,
                                 double rho,
                                 double mu,
                                 const std::pair<double, double>& scale,
                                 double f_periodic_corr,
                                 const int managing_rank) override;

    // \see CIBStrategy::constructGeometricMatrix() method.
    /*!
     * \brief Generate block-diagonal geometric matrix for the prototypical structures
     * identified by their indices.
     */
    void constructGeometricMatrix(const std::string& mat_name,
                                  Mat& geometric_mat,
                                  const std::vector<unsigned>& prototype_struct_ids,
                                  const bool initial_time,
                                  const int managing_rank) override;

    //\see CIBStrategy:: rotateArray() method.
    /*!
     * \brief Rotate vector using rotation matrix to/from the reference frame
     * of the structures.
     */
    void rotateArray(double* array,
                     const std::vector<unsigned>& struct_ids,
                     const bool use_transpose,
                     const int managing_rank,
                     const int depth) override;
    //\}

    /*!
     * Function to determine whether regridding should occur at the current time
     * step.
     */
    bool flagRegrid() const;

    /*
     * Set velocity physical boundary options
     */
    void setVelocityPhysBdryOp(IBTK::RobinPhysBdryPatchStrategy* u_phys_bdry_op);

    //////////////////////////////////////////////////////////////////////////////

protected:
    //////////////////////////////////////////////////////////////////////////////

    /*!
     * Functions to set constrained velocities of the structures.
     */
    std::vector<ConstrainedVelocityFcnsData> d_constrained_velocity_fcns_data;

    /*!
     * Functions to set net external force and torque on free moving structures.
     */
    std::vector<ExternalForceTorqueFcnData> d_ext_force_torque_fcn_data;

    /*!
     * Pre and post fluid solve call back functions and contexts.
     */
    std::vector<preprocessSolveFluidEqn_callbackfcn> d_prefluidsolve_callback_fcns;
    std::vector<void*> d_prefluidsolve_callback_fcns_ctx;

    /*!
     * Booleans to control spreading constraint force and interpolating
     * to Lagrangian velocities.
     */
    bool d_constraint_force_is_initialized = false, d_lag_velvec_is_initialized = false;

    /*!
     * Boolean to flag if time integrator needs regriding
     */
    bool d_time_integrator_needs_regrid = false;

    /*!
     * Eulerian variables.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > d_eul_lambda_var;
    int d_eul_lambda_idx = IBTK::invalid_index;

    /*!
     * Vector of Lagrnagian indices of all structures.
     */
    std::vector<std::pair<int, int> > d_struct_lag_idx_range;

    /*!
     * The object used to write out data for postprocessing by the visIt
     * visualization tool.
     */
    SAMRAI::tbox::Pointer<SAMRAI::appu::VisItDataWriter<NDIM> > d_visit_writer;

    /*!
     * Control printing of S[lambda]
     */
    bool d_output_eul_lambda = false;

    /*!
     *  Input/output.
     */
    int d_lambda_dump_interval = 0;
    std::ofstream d_lambda_stream;
    std::vector<std::string> d_reg_filename;
    std::vector<std::string> d_lambda_filename;

    // Velocity boundary operator.
    IBTK::RobinPhysBdryPatchStrategy* d_u_phys_bdry_op = nullptr;

    // If we are using steady Stokes solver.
    bool d_use_steady_stokes = false;

private:
    /*!
     * \brief Set additional values from input database.
     */
    void getFromInput(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db);

    /*!
     * \brief Get values from restart file.
     */
    void getFromRestart();

    /*!
     * \brief Compute center of mass of structures.
     */
    void computeCOMOfStructures(std::vector<Eigen::Vector3d>& center_of_mass,
                                std::vector<SAMRAI::tbox::Pointer<IBTK::LData> >& X_data);

    /*!
     * \brief Set regularization weight for Lagrangian markers.
     */
    void setRegularizationWeight(const int level_number);

    /*!
     * \brief Set initial Lambda for Lagrangian markers.
     */
    void setInitialLambda(const int level_number);

}; // CIBMethod
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBAMR_CIBMethod
