// ---------------------------------------------------------------------
//
// Copyright (c) 2019 - 2021 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#include "ibamr/AdvDiffConvectiveOperatorManager.h"
#include "ibamr/CFOneSidedUpperConvectiveOperator.h"
#include "ibamr/StaggeredStokesPhysicalBoundaryHelper.h"

#include "ibtk/HierarchyGhostCellInterpolation.h"

#include "Box.h"
#include "CartesianGridGeometry.h"
#include "CartesianPatchGeometry.h"
#include "CellData.h"
#include "FaceData.h"
#include "FaceIndex.h"
#include "Index.h"
#include "MultiblockDataTranslator.h"
#include "Patch.h"
#include "PatchLevel.h"
#include "SAMRAIVectorReal.h"
#include "SideData.h"
#include "SideIndex.h"
#include "SideIterator.h"
#include "Variable.h"
#include "VariableContext.h"
#include "VariableDatabase.h"
#include "tbox/Database.h"
#include "tbox/Utilities.h"

#include <Eigen/Dense>

#include "ibamr/app_namespaces.h" // IWYU pragma: keep

namespace SAMRAI
{
namespace solv
{
template <int DIM>
class RobinBcCoefStrategy;
} // namespace solv
} // namespace SAMRAI

#if (NDIM == 2)
#define UPPER_CONVECTIVE_OP_FC IBAMR_FC_FUNC_(upper_convective_op2d, UPPER_CONVECTIVE_OP2D)
#define SQRT_UPPER_CONVECTIVE_OP_FC IBAMR_FC_FUNC_(sqrt_upper_convective_op2d, SQRT_UPPER_CONVECTIVE_OP2D)
#define LOG_UPPER_CONVECTIVE_OP_FC IBAMR_FC_FUNC_(log_upper_convective_op2d, LOG_UPPER_CONVECTIVE_OP2D)
#endif
#if (NDIM == 3)
#define UPPER_CONVECTIVE_OP_FC IBAMR_FC_FUNC_(upper_convective_op3d, UPPER_CONVECTIVE_OP3D)
#define SQRT_UPPER_CONVECTIVE_OP_FC IBAMR_FC_FUNC_(sqrt_upper_convective_op3d, SQRT_UPPER_CONVECTIVE_OP3D)
#define LOG_UPPER_CONVECTIVE_OP_FC IBAMR_FC_FUNC_(log_upper_convective_op3d, LOG_UPPER_CONVECTIVE_OP3D)
#endif

extern "C"
{
#if (NDIM == 2)
    void UPPER_CONVECTIVE_OP_FC(const double*,
                                const double*,
                                const double*,
                                const int&,
                                const double*,
                                const int&,
                                const double*,
                                const int&,
                                const double*,
                                const int&,
                                const double*,
                                const int&,
                                const int&,
                                const int&,
                                const int&,
                                const int&);
    void SQRT_UPPER_CONVECTIVE_OP_FC(const double*,
                                     const double*,
                                     const double*,
                                     const int&,
                                     const double*,
                                     const int&,
                                     const double*,
                                     const int&,
                                     const double*,
                                     const int&,
                                     const double*,
                                     const int&,
                                     const int&,
                                     const int&,
                                     const int&,
                                     const int&);
    void LOG_UPPER_CONVECTIVE_OP_FC(const double*,
                                    const double*,
                                    const double*,
                                    const int&,
                                    const double*,
                                    const int&,
                                    const double*,
                                    const int&,
                                    const double*,
                                    const int&,
                                    const double*,
                                    const int&,
                                    const int&,
                                    const int&,
                                    const int&,
                                    const int&);
#endif
#if (NDIM == 3)
    void UPPER_CONVECTIVE_OP_FC(const double*,
                                const double*,
                                const double*,
                                const double*,
                                const int&,
                                const double*,
                                const int&,
                                const double*,
                                const int&,
                                const double*,
                                const int&,
                                const double*,
                                const int&,
                                const int&,
                                const int&,
                                const int&,
                                const int&,
                                const int&,
                                const int&);
    void SQRT_UPPER_CONVECTIVE_OP_FC(const double*,
                                     const double*,
                                     const double*,
                                     const double*,
                                     const int&,
                                     const double*,
                                     const int&,
                                     const double*,
                                     const int&,
                                     const double*,
                                     const int&,
                                     const double*,
                                     const int&,
                                     const int&,
                                     const int&,
                                     const int&,
                                     const int&,
                                     const int&,
                                     const int&);
    void LOG_UPPER_CONVECTIVE_OP_FC(const double*,
                                    const double*,
                                    const double*,
                                    const double*,
                                    const int&,
                                    const double*,
                                    const int&,
                                    const double*,
                                    const int&,
                                    const double*,
                                    const int&,
                                    const double*,
                                    const int&,
                                    const int&,
                                    const int&,
                                    const int&,
                                    const int&,
                                    const int&,
                                    const int&);
#endif
}

double
node_to_cell(const SAMRAI::pdat::CellIndex<NDIM>& idx, SAMRAI::pdat::NodeData<NDIM, double>& data)
{
#if (NDIM == 2)
    NodeIndex<NDIM> idx_ll(idx, IntVector<NDIM>(0, 0));
    NodeIndex<NDIM> idx_lu(idx, IntVector<NDIM>(0, 1));
    NodeIndex<NDIM> idx_ul(idx, IntVector<NDIM>(1, 0));
    NodeIndex<NDIM> idx_uu(idx, IntVector<NDIM>(1, 1));
    double ll = data(idx_ll), lu = data(idx_lu), ul = data(idx_ul), uu = data(idx_uu);
    return 0.25 * (ll + lu + ul + uu);
#endif
#if (NDIM == 3)
    double val = 0.0;
    for (int x = 0; x < 2; ++x)
    {
        for (int y = 0; y < 2; ++y)
        {
            for (int z = 0; z < 2; ++z)
            {
                NodeIndex<NDIM> n_idx(idx, IntVector<NDIM>(x, y, z));
                val += data(n_idx);
            }
        }
    }
    return val / 8.0;
#endif
}

double
cell_to_side(const SAMRAI::pdat::SideIndex<NDIM>& idx, SAMRAI::pdat::CellData<NDIM, double>& data)
{
    return 0.5 * (data(idx.toCell(0)) + data(idx.toCell(1)));
}

namespace IBAMR
{
CFOneSidedUpperConvectiveOperator::CFOneSidedUpperConvectiveOperator(
    const std::string& object_name,
    Pointer<CellVariable<NDIM, double> > Q_var,
    Pointer<Database> input_db,
    const std::string& difference_form,
    const std::vector<RobinBcCoefStrategy<NDIM>*>& Q_bc_coefs,
    const std::vector<RobinBcCoefStrategy<NDIM>*>& u_bc_coefs)
    : ConvectiveOperator(object_name, UNKNOWN_CONVECTIVE_DIFFERENCING_TYPE),
      d_Q_var(Q_var),
      d_u_adv_var(new SideVariable<NDIM, double>("Complex U var")),
      d_difference_form(difference_form),
      d_Q_bc_coefs(Q_bc_coefs),
      d_u_bc_coefs(u_bc_coefs)
{
    if (input_db)
    {
        d_interp_type = input_db->getStringWithDefault("interp_type", d_interp_type);
        if (input_db->keyExists("evolution_type"))
            d_evolve_type = string_to_enum<TensorEvolutionType>(input_db->getString("evolution_type"));
        if (input_db->keyExists("evolve_type"))
            d_evolve_type = string_to_enum<TensorEvolutionType>(input_db->getString("evolve_type"));
    }
    // Register some scratch variables
    auto var_db = VariableDatabase<NDIM>::getDatabase();
    d_Q_convec_idx = var_db->registerVariableAndContext(
        d_Q_var, var_db->getContext(d_object_name + "::CONVECTIVE"), IntVector<NDIM>(0));
    Pointer<VariableContext> new_cxt = var_db->getContext(d_object_name + "::U_ADV_CXT");
    d_u_scratch_idx = var_db->registerVariableAndContext(d_u_adv_var, new_cxt, IntVector<NDIM>(2));
    Pointer<VariableContext> src_cxt = var_db->getContext(d_object_name + "::SOURCE");
    d_s_idx = var_db->registerVariableAndContext(d_Q_var, src_cxt);
    auto convective_op_manager = AdvDiffConvectiveOperatorManager::getManager();
    d_convec_oper = convective_op_manager->allocateOperator(
        d_difference_form, d_object_name + "::convec_oper", d_Q_var, input_db, ADVECTIVE, d_Q_bc_coefs);
} // Constructor

CFOneSidedUpperConvectiveOperator::~CFOneSidedUpperConvectiveOperator()
{
    deallocateOperatorState();
    return;
} // Destructor

void
CFOneSidedUpperConvectiveOperator::applyConvectiveOperator(int Q_idx, int Y_idx)
{
    if (d_evolve_type != STANDARD)
    {
        TBOX_ERROR("CFOneSidedUpperConvectiveOperator currently only works with STANDARD evolution types.\n");
    }
    if (!d_s_fcn)
    {
        TBOX_ERROR("CFOneSidedUpperConvectiveOperator::applyConvectiveOperator():\n"
                   << "  Source function must be register prior to call to "
                      "applyConvectiveOperator\n");
    }
    if (!d_is_initialized)
    {
        TBOX_ERROR("CFOneSidedUpperConvectiveOperator::applyConvectiveOperator():\n"
                   << "  operator must be initialized prior to call to "
                      "applyConvectiveOperator\n");
    }
    // Set up velocity information:
    d_convec_oper->setAdvectionVelocity(d_u_idx);
    d_convec_oper->setSolutionTime(d_solution_time);
    Pointer<CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();
    // Copy data to side centered velocity field
    //
    // TODO: This is only done because we currently only have operators to
    // fill in physical boundary conditions for cell- and side-centered
    // quantities.
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
            Pointer<FaceData<NDIM, double> > u_f_data = patch->getPatchData(d_u_idx);
            Pointer<SideData<NDIM, double> > u_s_data = patch->getPatchData(d_u_scratch_idx);
            const Box<NDIM> box = patch->getBox();
            for (int axis = 0; axis < NDIM; ++axis)
            {
                for (SideIterator<NDIM> i(box, axis); i; i++)
                {
                    const SideIndex<NDIM>& si = i();
                    FaceIndex<NDIM> fi(si.toCell(0), axis, 1);
                    (*u_s_data)(si) = (*u_f_data)(fi);
                }
            }
        }
    }
    // Fill boundary conditions for side centered velocity
    using InterpolationTransactionComponent = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
    std::vector<InterpolationTransactionComponent> ghost_cell_components(1);
    ghost_cell_components[0] = InterpolationTransactionComponent(
        d_u_scratch_idx, "NONE", true, "CUBIC_COARSEN", "LINEAR", false, d_u_bc_coefs, NULL, d_interp_type);
    HierarchyGhostCellInterpolation ghost_fill_op;
    ghost_fill_op.initializeOperatorState(ghost_cell_components, d_hierarchy);
    StaggeredStokesPhysicalBoundaryHelper::setupBcCoefObjects(d_u_bc_coefs, nullptr, d_u_scratch_idx, -1, false);
    ghost_fill_op.fillData(d_solution_time);
    StaggeredStokesPhysicalBoundaryHelper::resetBcCoefObjects(d_u_bc_coefs, nullptr);

    d_convec_oper->applyConvectiveOperator(Q_idx, d_Q_convec_idx);

    d_s_fcn->setPatchDataIndex(Q_idx);
    d_s_fcn->setDataOnPatchHierarchy(d_s_idx, d_Q_var, d_hierarchy, d_solution_time, false, d_coarsest_ln, d_finest_ln);

    for (int level_num = d_coarsest_ln; level_num <= d_finest_ln; ++level_num)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(level_num);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());

            const Pointer<CartesianPatchGeometry<NDIM> > p_geom = patch->getPatchGeometry();
            const double* const dx = p_geom->getDx();
            const double* const xlow = p_geom->getXLower();
            const Box<NDIM>& patch_box = patch->getBox();
            const hier::Index<NDIM>& idx_low = patch_box.lower();

            Pointer<SideData<NDIM, double> > u_data = patch->getPatchData(d_u_scratch_idx); // velocity data
            Pointer<CellData<NDIM, double> > Q_data = patch->getPatchData(Q_idx);           // Stress data
            Pointer<CellData<NDIM, double> > C_data = patch->getPatchData(d_Q_convec_idx);  // u.grad(C) data
            Pointer<CellData<NDIM, double> > S_data = patch->getPatchData(d_s_idx);         // Relaxation data
            Pointer<CellData<NDIM, double> > ls_data = patch->getPatchData(d_ls_idx);       // Level set data
            Pointer<CellData<NDIM, double> > ret_data = patch->getPatchData(Y_idx);         // Data to fill in

            for (CellIterator<NDIM> ci(patch_box); ci; ci++)
            {
                const CellIndex<NDIM>& idx = ci();

                // Compute gradient of velocity
                const double ls = (*ls_data)(idx);
                VectorNd eval_pt;
                for (int d = 0; d < NDIM; ++d)
                    eval_pt[d] = xlow[d] + dx[d] * (static_cast<double>(idx(d) - idx_low(d)) + 0.5);
                MatrixNd grad_u;
                // Using quadratic polynomials, so use 14 closest points in 2D.
                size_t poly_size = 2 * NDIM + 2;
                size_t stencil_size = 14;
                for (int axis = 0; axis < NDIM; ++axis)
                {
                    // Start with two side indices on cell
                    std::vector<SideIndex<NDIM> > test_idxs = { SideIndex<NDIM>(idx, axis, 0),
                                                                SideIndex<NDIM>(idx, axis, 1) };
                    std::vector<VectorNd> X_vals;
                    std::vector<double> u_vals;
                    unsigned int i = 0;
                    // Use a flooding algorithm to grab points
                    // TODO: This produces a biased stencil. We should reconsider the best thing to do here...
                    while (X_vals.size() < stencil_size)
                    {
                        TBOX_ASSERT(i < test_idxs.size());
                        SideIndex<NDIM> test_idx = test_idxs[i];
                        // Check we are of the same sign as where we want to evaluate the gradient.
                        if (cell_to_side(test_idx, *ls_data) * ls > 0.0)
                        {
                            u_vals.push_back((*u_data)(test_idx));
                            VectorNd x_cent = VectorNd::Zero();
                            for (int d = 0; d < NDIM; ++d)
                                x_cent[d] = xlow[d] + dx[d] * (static_cast<double>(test_idx[d] - idx_low(d)) +
                                                               (d == axis ? 0.0 : 0.5));
                            X_vals.push_back(x_cent);
                        }

                        // Add neighboring points to new_idxs
                        for (int axis = 0; axis < NDIM; ++axis)
                        {
                            for (int upperlower = 0; upperlower < 2; ++upperlower)
                            {
                                IntVector<NDIM> shft(0);
                                shft(axis) = upperlower == 0 ? -1 : 1;
                                SideIndex<NDIM> new_idx(test_idx + shft);
                                // Make sure we haven't added searched this point yet
                                if (std::find(test_idxs.begin(), test_idxs.end(), new_idx) == test_idxs.end())
                                    test_idxs.push_back(new_idx);
                            }
                        }
                        ++i;
                    }

                    // At this point, we have the stencil needed to compute grad(u)
                    // Now fit a spline to the data
                    const int m = u_vals.size();
                    IBTK::MatrixXd A(IBTK::MatrixXd::Zero(m, m));
                    IBTK::MatrixXd B(IBTK::MatrixXd::Zero(m, poly_size));
                    IBTK::MatrixXd U(IBTK::MatrixXd::Zero(m + poly_size, NDIM)); // Note we need three derivatives
                    for (size_t i = 0; i < u_vals.size(); ++i)
                    {
                        for (size_t j = 0; j < u_vals.size(); ++j)
                        {
                            const VectorNd X = X_vals[i] - X_vals[j];
                            A(i, j) = std::pow(X.norm(), 5.0);
                        }
                        // Polynomial terms. Note we shift them to the evaluation point and scale them by dx[0]
                        // Constant
                        B(i, 0) = 1.0;
                        // Linear
                        for (int d = 0; d < NDIM; ++d) B(i, d + 1) = (X_vals[i](d) - eval_pt(d)) / dx[0];
                        // Quadratic
                        B(i, NDIM + 1) = (X_vals[i](0) - eval_pt(0)) * (X_vals[i](0) - eval_pt(0)) / dx[0];
                        B(i, NDIM + 2) = (X_vals[i](1) - eval_pt(1)) * (X_vals[i](1) - eval_pt(1)) / dx[0];
                        B(i, NDIM + 3) = (X_vals[i](0) - eval_pt(0)) * (X_vals[i](1) - eval_pt(1)) / dx[0];
                        // RHS, L(rbf)
                        for (int d = 0; d < NDIM; ++d)
                            U(i, d) = 5.0 * (eval_pt(d) - X_vals[i](d)) * std::pow((eval_pt - X_vals[i]).norm(), 3.0);
                    }

                    // Only non-zero term will be in linear term because of the shift.
                    // Divide by the scale too.
                    for (int d = 0; d < NDIM; ++d) U(stencil_size + d + 1, d) = 1.0 / dx[0];

                    // Now build big matrix
                    IBTK::MatrixXd final_mat(IBTK::MatrixXd::Zero(stencil_size + poly_size, stencil_size + poly_size));
                    final_mat.block(0, 0, stencil_size, stencil_size) = A;
                    final_mat.block(0, stencil_size, stencil_size, poly_size) = B;
                    final_mat.block(stencil_size, 0, poly_size, stencil_size) = B.transpose();
                    final_mat.block(stencil_size, stencil_size, poly_size, poly_size).setZero();
                    IBTK::MatrixXd x = final_mat.colPivHouseholderQr().solve(U);

                    // Now evaluate stencil
                    const IBTK::MatrixXd& weights = x.block(0, 0, stencil_size, NDIM);
                    for (int d = 0; d < NDIM; ++d)
                    {
                        for (size_t i = 0; i < stencil_size; ++i)
                        {
                            grad_u(d, i) = weights(i, d) * u_vals[i];
                        }
                    }
                }
            }
        } // end Patch loop
    }     // end Level loop
} // applyConvectiveOperator

void
CFOneSidedUpperConvectiveOperator::initializeOperatorState(const SAMRAIVectorReal<NDIM, double>& in,
                                                           const SAMRAIVectorReal<NDIM, double>& out)
{
    if (d_is_initialized) deallocateOperatorState();
    // Get Hierarchy Information
    d_hierarchy = in.getPatchHierarchy();
    d_coarsest_ln = in.getCoarsestLevelNumber();
    d_finest_ln = in.getFinestLevelNumber();
#if !defined(NDEBUG)
    TBOX_ASSERT(d_hierarchy == out.getPatchHierarchy());
    TBOX_ASSERT(d_coarsest_ln == out.getCoarsestLevelNumber());
    TBOX_ASSERT(d_finest_ln == out.getFinestLevelNumber());
#else
    NULL_USE(out);
#endif
    d_convec_oper->initializeOperatorState(in, out);
    d_convec_oper->setAdvectionVelocity(d_u_idx);
    // Allocate Patch Data
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        if (!level->checkAllocated(d_u_scratch_idx)) level->allocatePatchData(d_u_scratch_idx);
        if (!level->checkAllocated(d_Q_convec_idx)) level->allocatePatchData(d_Q_convec_idx);
        if (!level->checkAllocated(d_s_idx)) level->allocatePatchData(d_s_idx);
    }
    d_is_initialized = true;
    return;
} // initializeOperatorState

void
CFOneSidedUpperConvectiveOperator::deallocateOperatorState()
{
    if (!d_is_initialized) return;
    d_convec_oper->deallocateOperatorState();
    // Deallocate scratch data
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        if (level->checkAllocated(d_Q_convec_idx)) level->deallocatePatchData(d_Q_convec_idx);
        if (level->checkAllocated(d_u_scratch_idx)) level->deallocatePatchData(d_u_scratch_idx);
        if (level->checkAllocated(d_s_idx)) level->deallocatePatchData(d_s_idx);
    }
    d_is_initialized = false;
    return;
} // deallocateOperatorState

void
CFOneSidedUpperConvectiveOperator::registerSourceFunction(Pointer<CFRelaxationOperator> source_fcn)
{
    d_s_fcn = source_fcn;
}

} // namespace IBAMR
