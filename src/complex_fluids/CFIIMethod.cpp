// ---------------------------------------------------------------------
//
// Copyright (c) 2018 - 2022 by the IBAMR developers
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

#include "ibamr/CFIIMethod.h"
#include "ibamr/INSHierarchyIntegrator.h"
#include "ibamr/ibamr_utilities.h"

#include "ibtk/FEDataInterpolation.h"
#include "ibtk/FEDataManager.h"
#include "ibtk/IBTK_CHKERRQ.h"
#include "ibtk/IndexUtilities.h"
#include "ibtk/LEInteractor.h"
#include "ibtk/SAMRAIDataCache.h"
#include "ibtk/ibtk_utilities.h"
#include "ibtk/libmesh_utilities.h"

#include "BasePatchHierarchy.h"
#include "BasePatchLevel.h"
#include "Box.h"
#include "CartesianGridGeometry.h"
#include "CartesianPatchGeometry.h"
#include "CellData.h"
#include "CellIndex.h"
#include "GriddingAlgorithm.h"
#include "Index.h"
#include "IntVector.h"
#include "LoadBalancer.h"
#include "Patch.h"
#include "PatchData.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "RefineSchedule.h"
#include "SideData.h"
#include "SideGeometry.h"
#include "SideIndex.h"
#include "tbox/Array.h"
#include "tbox/Database.h"
#include "tbox/MathUtilities.h"
#include "tbox/PIO.h"
#include "tbox/Pointer.h"
#include "tbox/RestartManager.h"
#include "tbox/SAMRAI_MPI.h"
#include "tbox/Utilities.h"

#include "libmesh/boundary_info.h"
#include "libmesh/compare_types.h"
#include "libmesh/dense_vector.h"
#include "libmesh/dof_map.h"
#include "libmesh/edge.h"
#include "libmesh/elem.h"
#include "libmesh/enum_fe_family.h"
#include "libmesh/enum_order.h"
#include "libmesh/enum_parallel_type.h"
#include "libmesh/enum_quadrature_type.h"
#include "libmesh/enum_xdr_mode.h"
#include "libmesh/equation_systems.h"
#include "libmesh/face.h"
#include "libmesh/fe_base.h"
#include "libmesh/fe_type.h"
#include "libmesh/fem_context.h"
#include "libmesh/id_types.h"
#include "libmesh/libmesh_common.h"
#include "libmesh/libmesh_config.h"
#include "libmesh/libmesh_version.h"
#include "libmesh/mesh_base.h"
#include "libmesh/node.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/petsc_vector.h"
#include "libmesh/point.h"
#include "libmesh/quadrature.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/system.h"
#include "libmesh/tensor_value.h"
#include "libmesh/type_vector.h"
#include "libmesh/variant_filter_iterator.h"
#include "libmesh/vector_value.h"

#include "petscvec.h"

#include <ibamr/app_namespaces.h>

IBTK_DISABLE_EXTRA_WARNINGS
#include <boost/math/special_functions/round.hpp>
#include <boost/multi_array.hpp>
IBTK_ENABLE_EXTRA_WARNINGS

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <limits>
#include <map>
#include <memory>
#include <ostream>
#include <queue>
#include <string>
#include <utility>
#include <vector>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
static Timer* t_compute_lagrangian_force;
static Timer* t_spread_force;
// Version of CFIIMethod restart file data.
static const int CF_IIM_VERSION = 0;

std::string
libmesh_restart_file_name(const std::string& restart_dump_dirname,
                          unsigned int time_step_number,
                          unsigned int part,
                          const std::string& extension)
{
    std::ostringstream file_name_prefix;
    file_name_prefix << restart_dump_dirname << "/libmesh_data_part_" << part << "." << std::setw(6)
                     << std::setfill('0') << std::right << time_step_number << "." << extension;
    return file_name_prefix.str();
}

static const int TENSOR_SIZE = NDIM * (NDIM + 1) / 2;
} // namespace

const std::string CFIIMethod::EXTRA_STRESS_JUMP_SYSTEM_NAME = "extra stress jump system";
const std::string CFIIMethod::EXTRA_STRESS_IN_SYSTEM_NAME = "extra stress in system";
const std::string CFIIMethod::EXTRA_STRESS_OUT_SYSTEM_NAME = "extra stress out system";

/////////////////////////////// PUBLIC ///////////////////////////////////////

CFIIMethod::CFIIMethod(const std::string& object_name,
                       Pointer<Database> input_db,
                       MeshBase* mesh,
                       int max_level_number,
                       bool register_for_restart,
                       const std::string& restart_read_dirname,
                       unsigned int restart_restore_number)
    : IIMethod(object_name,
               input_db,
               mesh,
               max_level_number,
               register_for_restart,
               restart_read_dirname,
               restart_restore_number)
{
    commonConstructor(object_name,
                      input_db,
                      std::vector<MeshBase*>(1, mesh),
                      max_level_number,
                      register_for_restart,
                      restart_read_dirname,
                      restart_restore_number);
    return;
} // CFIIMethod

CFIIMethod::CFIIMethod(const std::string& object_name,
                       Pointer<Database> input_db,
                       const std::vector<MeshBase*>& meshes,
                       int max_level_number,
                       bool register_for_restart,
                       const std::string& restart_read_dirname,
                       unsigned int restart_restore_number)
    : IIMethod(object_name,
               input_db,
               meshes,
               max_level_number,
               register_for_restart,
               restart_read_dirname,
               restart_restore_number)
{
    commonConstructor(object_name,
                      input_db,
                      meshes,
                      max_level_number,
                      register_for_restart,
                      restart_read_dirname,
                      restart_restore_number);
    return;
} // CFIIMethod

CFIIMethod::~CFIIMethod()
{
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        delete d_equation_systems[part];
    }
    if (d_registered_for_restart)
    {
        RestartManager::getManager()->unregisterRestartItem(d_object_name);
        d_registered_for_restart = false;
    }
    return;
} // ~CFIIMethod

void
CFIIMethod::computeStressJumps(const int sig_idx, const int ls_idx, const double data_time, const unsigned int part)
{
    Pointer<PatchHierarchy<NDIM> > patch_hierarchy = d_fe_data_managers[part]->getPatchHierarchy();

    // Extract the FE systems and DOF maps, and setup the FE object.
    EquationSystems* equation_systems = d_fe_data_managers[part]->getEquationSystems();
    const MeshBase& mesh = equation_systems->get_mesh();
    const unsigned int dim = mesh.mesh_dimension();
    std::unique_ptr<QBase> qrule;

    System& X_system = equation_systems->get_system(COORDS_SYSTEM_NAME);
    DofMap& X_dof_map = X_system.get_dof_map();
    std::vector<std::vector<unsigned int> > X_dof_indices(NDIM);
    FEType X_fe_type = X_dof_map.variable_type(0);
    NumericVector<double>* X_ghost_vec = X_system.current_local_solution.get();
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        TBOX_ASSERT(X_dof_map.variable_type(d) == X_fe_type);
    }
    std::unique_ptr<FEBase> fe_X = FEBase::build(dim, X_fe_type);
    const std::vector<std::vector<double> >& phi_X = fe_X->get_phi();
    std::array<const std::vector<std::vector<double> >*, NDIM - 1> dphi_dxi_X;
    dphi_dxi_X[0] = &fe_X->get_dphidxi();
    if (NDIM > 2) dphi_dxi_X[1] = &fe_X->get_dphideta();

    System& sig_jump_sys = equation_systems->get_system(EXTRA_STRESS_JUMP_SYSTEM_NAME);
    FEDataManager::SystemDofMapCache& sig_jump_dof_map_cache =
        *d_fe_data_managers[part]->getDofMapCache(EXTRA_STRESS_JUMP_SYSTEM_NAME);
    NumericVector<double>* sig_jump_vec = sig_jump_sys.solution.get();

    System& sig_in_sys = equation_systems->get_system(EXTRA_STRESS_IN_SYSTEM_NAME);
    FEDataManager::SystemDofMapCache& sig_in_dof_map_cache =
        *d_fe_data_managers[part]->getDofMapCache(EXTRA_STRESS_IN_SYSTEM_NAME);
    NumericVector<double>* sig_in_vec = sig_in_sys.solution.get();
    std::unique_ptr<NumericVector<double> > sig_in_rhs_vec = std::move(sig_in_vec->zero_clone());

    System& sig_out_sys = equation_systems->get_system(EXTRA_STRESS_OUT_SYSTEM_NAME);
    FEDataManager::SystemDofMapCache& sig_out_dof_map_cache =
        *d_fe_data_managers[part]->getDofMapCache(EXTRA_STRESS_OUT_SYSTEM_NAME);
    std::vector<unsigned int> sig_out_dof_indices;
    NumericVector<double>* sig_out_vec = sig_out_sys.solution.get();
    std::unique_ptr<NumericVector<double> > sig_out_rhs_vec = std::move(sig_out_vec->zero_clone());

    std::unique_ptr<FEBase> fe_sig = FEBase::build(dim, sig_jump_sys.get_dof_map().variable_type(0));
    const std::vector<std::vector<double> >& phi_sig = fe_sig->get_phi();
    const std::vector<double>& JxW_sig = fe_sig->get_JxW();

    const std::vector<std::vector<Elem*> >& active_patch_element_map =
        d_fe_data_managers[part]->getActivePatchElementMap();

    boost::multi_array<double, 2> x_node;
    std::array<VectorValue<double>, 2> dx_dxi;

    // DenseVectors for filling in rhs
    std::array<DenseVector<double>, TENSOR_SIZE> sig_in_rhs_e;
    std::array<DenseVector<double>, TENSOR_SIZE> sig_out_rhs_e;

    // Assume the structure is on the finest level
    Pointer<PatchLevel<NDIM> > level =
        d_hierarchy->getPatchLevel(d_fe_data_managers[part]->getFinestPatchLevelNumber());
    const Pointer<CartesianGridGeometry<NDIM> > grid_geom = level->getGridGeometry();
    VectorValue<double> tau1, tau2, n;
    X_ghost_vec->close();
    int local_patch_num = 0;
    for (PatchLevel<NDIM>::Iterator p(level); p; p++, ++local_patch_num)
    {
        // The relevant collection of elements.
        const std::vector<Elem*>& patch_elems = active_patch_element_map[local_patch_num];
        const size_t num_active_patch_elems = patch_elems.size();
        if (!num_active_patch_elems) continue;
        const Pointer<Patch<NDIM> > patch = level->getPatch(p());
        const hier::Index<NDIM>& idx_low = patch->getBox().lower();
        const Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
        const double* const x_lower = pgeom->getXLower();
        const double* const x_upper = pgeom->getXUpper();
        const double* const dx = pgeom->getDx();
        const double dx_min = *std::min_element(dx, dx + NDIM);

        Pointer<CellData<NDIM, double> > sig_data = patch->getPatchData(sig_idx);

        double diag_dis = 0.0;
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            diag_dis += dx[d] * dx[d];
        }

        if (d_sig_calc_width == 0.0)
        {
            TBOX_ERROR(d_object_name << ": The width for the stress jump hasn't been set up!" << std::endl);
        }
        const double dh = d_sig_calc_width * sqrt(diag_dis);

        // Loop over the elements and compute the positions of the quadrature points.
        unsigned int qp_offset = 0;
        for (const auto& elem : patch_elems)
        {
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                X_dof_map.dof_indices(elem, X_dof_indices[d], d);
            }
            get_values_for_interpolation(x_node, *X_ghost_vec, X_dof_indices);

            const boost::multi_array<dof_id_type, 2>& sig_in_dof_indices = sig_in_dof_map_cache.dof_indices(elem);
            const boost::multi_array<dof_id_type, 2>& sig_out_dof_indices = sig_out_dof_map_cache.dof_indices(elem);

            for (unsigned int d = 0; d < TENSOR_SIZE; ++d)
            {
                sig_in_rhs_e[d].resize(static_cast<int>(sig_in_dof_indices[d].size()));
                sig_out_rhs_e[d].resize(static_cast<int>(sig_out_dof_indices[d].size()));
            }

            const bool qrule_changed =
                FEDataManager::updateInterpQuadratureRule(qrule, d_default_interp_spec, elem, x_node, dx_min);

            if (qrule_changed)
            {
                fe_X->attach_quadrature_rule(qrule.get());
                fe_sig->attach_quadrature_rule(qrule.get());
            }
            fe_X->reinit(elem);
            fe_sig->reinit(elem);

            const unsigned int n_node = elem->n_nodes();
            const unsigned int n_qp = qrule->n_points();

            // Interpolate X at all of the quadrature points
            // via accumulation, i.e., X(qp) = sum_k X_k * phi_k(qp) for
            // each qp.
            for (unsigned int qp = 0; qp < n_qp; ++qp)
            {
                // Zero out the values of X and n prior to accumulation.
                VectorValue<double> x_qp, x_in, x_out, N_qp;
                x_qp.zero();
                std::array<double, TENSOR_SIZE> sig_in, sig_out;
                for (unsigned int k = 0; k < NDIM - 1; ++k)
                {
                    interpolate(dx_dxi[k], qp, x_node, *dphi_dxi_X[k]);
                }
                if (NDIM == 2)
                {
                    dx_dxi[1] = VectorValue<double>(0.0, 0.0, 1.0);
                }
                N_qp = (dx_dxi[0].cross(dx_dxi[1])).unit();

                for (unsigned int i = 0; i < NDIM; ++i)
                {
                    for (unsigned int k = 0; k < n_node; ++k) x_qp(i) += x_node[k][i] * phi_X[k][qp];
                }
                // We don't interpolate the stresses to a small distance away from the quadrature point. This helps
                // regularize the interpolation
                // TODO: Test if this is necessary.
                x_in = x_qp - N_qp * dh;
                x_out = x_qp + N_qp * dh;
                // Interpolate values from the Cartesian grid patch to the quadrature
                // points.
                // Note: Values are interpolated only to those quadrature points that
                // are within the patch interior
                // Start with out
                CellIndex<NDIM> idx_cent = IndexUtilities::getCellIndex(&x_out(0), pgeom, patch->getBox());
                std::vector<CellIndex<NDIM> > stencil;
                std::vector<VectorValue<double> > stencil_pts;
                findInterpPts(idx_cent, ls_idx, patch, d_stencil_size, stencil_pts, stencil);
                std::vector<double> stencil_wgts = findLinearInterpWeights(x_out, stencil_pts, dx_min);
                for (int d = 0; d < TENSOR_SIZE; ++d)
                    sig_out[d] = evaluateInterpolant(*sig_data, stencil, stencil_wgts, d);

                idx_cent = IndexUtilities::getCellIndex(&x_in(0), pgeom, patch->getBox());
                findInterpPts(idx_cent, ls_idx, patch, d_stencil_size, stencil_pts, stencil);
                stencil_wgts = findLinearInterpWeights(x_in, stencil_pts, dx_min);
                for (int d = 0; d < TENSOR_SIZE; ++d)
                    sig_in[d] = evaluateInterpolant(*sig_data, stencil, stencil_wgts, d);

                // We have the value on the quadrature point. Fill in appropriate RHS.
                for (unsigned int k = 0; k < sig_out_dof_indices[0].size(); ++k)
                {
                    const double p_JxW = phi_sig[k][qp] * JxW_sig[qp];
                    for (unsigned int d = 0; d < TENSOR_SIZE; ++d)
                    {
                        sig_in_rhs_e[d](k) += sig_in[d] * p_JxW;
                        sig_out_rhs_e[d](k) += sig_out[d] * p_JxW;
                    }
                }
            }

            for (unsigned int d = 0; d < TENSOR_SIZE; ++d)
            {
                // Fill in the big rhs vector
                std::vector<dof_id_type> dof_id_scr;
                copy_dof_ids_to_vector(d, sig_in_dof_indices, dof_id_scr);
                sig_in_rhs_vec->add_vector(sig_in_rhs_e[d], dof_id_scr);

                copy_dof_ids_to_vector(d, sig_out_dof_indices, dof_id_scr);
                sig_out_rhs_vec->add_vector(sig_out_rhs_e[d], dof_id_scr);
            }
        }
    }

    // We've computed what we need on every element. Now we need to project
    sig_in_rhs_vec->close();
    sig_out_rhs_vec->close();

    d_fe_data_managers[part]->computeL2Projection(
        *sig_in_vec, *sig_in_rhs_vec, EXTRA_STRESS_IN_SYSTEM_NAME, d_default_interp_spec.use_consistent_mass_matrix);
    sig_in_vec->close();

    d_fe_data_managers[part]->computeL2Projection(
        *sig_out_vec, *sig_out_rhs_vec, EXTRA_STRESS_OUT_SYSTEM_NAME, d_default_interp_spec.use_consistent_mass_matrix);
    sig_out_vec->close();

    // Finally, compute the jump vector
    sig_jump_vec->zero();
    sig_jump_vec->add(*sig_in_vec);
    sig_jump_vec->add(-1.0, *sig_out_vec);
    sig_jump_vec->close();

    sig_jump_sys.update();
    sig_in_sys.update();
    sig_out_sys.update();

    X_ghost_vec->close();
    return;

} // computeStressJumps

void
CFIIMethod::computeLagrangianForce(const double data_time)
{
    IBAMR_TIMER_START(t_compute_lagrangian_force);
    TBOX_ASSERT(MathUtilities<double>::equalEps(data_time, d_half_time));
    batch_vec_ghost_update(d_X_half_vecs, INSERT_VALUES, SCATTER_FORWARD);
    for (unsigned part = 0; part < d_num_parts; ++part)
    {
        EquationSystems* equation_systems = d_fe_data_managers[part]->getEquationSystems();
        MeshBase& mesh = equation_systems->get_mesh();
        const unsigned int dim = mesh.mesh_dimension();

        // Setup global and elemental right-hand-side vectors.
        NumericVector<double>* F_vec = d_F_half_vecs[part];
        std::unique_ptr<NumericVector<double> > F_rhs_vec = F_vec->zero_clone();
        std::array<DenseVector<double>, NDIM> F_rhs_e;
        VectorValue<double>& F_integral = d_lag_surface_force_integral[part];
        F_integral.zero();

        NumericVector<double>* X_vec = d_X_half_vecs[part];
        double surface_area = 0.0;

        NumericVector<double>* P_jump_vec = d_use_pressure_jump_conditions ? d_P_jump_half_vecs[part] : nullptr;
        std::unique_ptr<NumericVector<double> > P_jump_rhs_vec;
        DenseVector<double> P_jump_rhs_e;
        if (d_use_pressure_jump_conditions)
        {
            P_jump_rhs_vec = P_jump_vec->zero_clone();
        }
        double P_jump_rhs_integral = 0.0;

        std::array<NumericVector<double>*, NDIM> DU_jump_vec;
        std::array<std::unique_ptr<NumericVector<double> >, NDIM> DU_jump_rhs_vec;
        std::array<std::array<DenseVector<double>, NDIM>, NDIM> DU_jump_rhs_e;
        if (d_use_velocity_jump_conditions)
        {
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                DU_jump_vec[d] = d_DU_jump_half_vecs[part][d];
                DU_jump_rhs_vec[d] = DU_jump_vec[d]->zero_clone();
            }
        }

        // Extract the FE systems and DOF maps, and setup the FE objects.
        System& F_system = equation_systems->get_system(FORCE_SYSTEM_NAME);
        const DofMap& F_dof_map = F_system.get_dof_map();
        FEDataManager::SystemDofMapCache& F_dof_map_cache =
            *d_fe_data_managers[part]->getDofMapCache(FORCE_SYSTEM_NAME);
        FEType F_fe_type = F_dof_map.variable_type(0);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            TBOX_ASSERT(F_dof_map.variable_type(d) == F_fe_type);
        }

        System& X_system = equation_systems->get_system(COORDS_SYSTEM_NAME);
        const DofMap& X_dof_map = X_system.get_dof_map();
        FEDataManager::SystemDofMapCache& X_dof_map_cache =
            *d_fe_data_managers[part]->getDofMapCache(COORDS_SYSTEM_NAME);
        FEType X_fe_type = X_dof_map.variable_type(0);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            TBOX_ASSERT(X_dof_map.variable_type(d) == X_fe_type);
        }
        TBOX_ASSERT(X_fe_type == F_fe_type);
        NumericVector<double>& X0_vec = X_system.get_vector("INITIAL_COORDINATES");
        System* P_jump_system;
        const DofMap* P_jump_dof_map = NULL;
        FEDataManager::SystemDofMapCache* P_jump_dof_map_cache = NULL;
        FEType P_jump_fe_type = INVALID_FE;
        if (d_use_pressure_jump_conditions)
        {
            P_jump_system = &equation_systems->get_system(PRESSURE_JUMP_SYSTEM_NAME);
            P_jump_dof_map = &P_jump_system->get_dof_map();
            P_jump_dof_map_cache = d_fe_data_managers[part]->getDofMapCache(PRESSURE_JUMP_SYSTEM_NAME);
            P_jump_fe_type = P_jump_dof_map->variable_type(0);
        }

        std::array<System*, NDIM> DU_jump_system;
        std::array<DofMap*, NDIM> DU_jump_dof_map;
        std::array<FEDataManager::SystemDofMapCache*, NDIM> DU_jump_dof_map_cache;
        FEType DU_jump_fe_type = INVALID_FE;
        if (d_use_velocity_jump_conditions)
        {
            for (unsigned int i = 0; i < NDIM; ++i)
            {
                DU_jump_system[i] = &equation_systems->get_system(VELOCITY_JUMP_SYSTEM_NAME[i]);
                DU_jump_dof_map[i] = &DU_jump_system[i]->get_dof_map();
                DU_jump_dof_map_cache[i] = d_fe_data_managers[part]->getDofMapCache(VELOCITY_JUMP_SYSTEM_NAME[i]);
            }
            DU_jump_fe_type = DU_jump_dof_map[0]->variable_type(0);
            for (unsigned int i = 0; i < NDIM; ++i)
            {
                for (unsigned int j = 0; j < NDIM; ++j)
                {
                    TBOX_ASSERT(DU_jump_dof_map[i]->variable_type(j) == DU_jump_fe_type);
                }
            }
        }

        // The P_jump_fe_type and DU_jump_fe_type are equal only if we are applying both jumps
        // or none of the jumps.
        if (d_use_pressure_jump_conditions && d_use_velocity_jump_conditions)
        {
            TBOX_ASSERT(P_jump_fe_type == DU_jump_fe_type);
        }

        // Pull out sigma jumps
        System& sig_jump_sys = equation_systems->get_system<System>(EXTRA_STRESS_JUMP_SYSTEM_NAME);
        FEDataManager::SystemDofMapCache& sig_jump_dof_map_cache =
            *d_fe_data_managers[part]->getDofMapCache(EXTRA_STRESS_JUMP_SYSTEM_NAME);
        NumericVector<double>* sig_jump_vec = sig_jump_sys.current_local_solution.get();

        std::unique_ptr<QBase> qrule = QBase::build(d_default_quad_type[part], dim, d_default_quad_order[part]);

        std::unique_ptr<FEBase> fe_X = FEBase::build(dim, X_fe_type);
        fe_X->attach_quadrature_rule(qrule.get());
        const std::vector<double>& JxW = fe_X->get_JxW();
        const std::vector<std::vector<double> >& phi_X = fe_X->get_phi();
        std::array<const std::vector<std::vector<double> >*, NDIM - 1> dphi_dxi_X;
        dphi_dxi_X[0] = &fe_X->get_dphidxi();
        if (NDIM > 2) dphi_dxi_X[1] = &fe_X->get_dphideta();

        FEType fe_jump_type = INVALID_FE;
        if (d_use_pressure_jump_conditions)
        {
            fe_jump_type = P_jump_fe_type;
        }
        else
        {
            fe_jump_type = DU_jump_fe_type;
        }

        std::unique_ptr<FEBase> fe_jump = FEBase::build(dim, fe_jump_type);
        fe_jump->attach_quadrature_rule(qrule.get());
        const std::vector<std::vector<double> >& phi_jump = fe_jump->get_phi();

        FEDataInterpolation fe_interpolator(dim, d_fe_data_managers[part]->getFEData());
        fe_interpolator.attachQuadratureRule(qrule.get());

        std::vector<size_t> surface_force_fcn_system_idxs;
        fe_interpolator.setupInterpolatedSystemDataIndexes(
            surface_force_fcn_system_idxs, d_lag_surface_force_fcn_data[part].system_data, equation_systems);
        std::vector<size_t> surface_pressure_fcn_system_idxs;
        fe_interpolator.setupInterpolatedSystemDataIndexes(
            surface_pressure_fcn_system_idxs, d_lag_surface_pressure_fcn_data[part].system_data, equation_systems);
        fe_interpolator.init();

        std::vector<const std::vector<double>*> surface_force_var_data, surface_pressure_var_data;
        std::vector<const std::vector<VectorValue<double> >*> surface_force_grad_var_data,
            surface_pressure_grad_var_data;

        // Loop over the elements to compute the right-hand side vector.
        boost::multi_array<double, 2> X_node, x_node, sig_jump_node;
        double DU[NDIM][NDIM];
        TensorValue<double> FF, sig_jump;
        VectorValue<double> F, F_b, F_s, F_qp, N, X, n, x;
        std::array<VectorValue<double>, 2> dX_dxi, dx_dxi;
        std::vector<libMesh::dof_id_type> dof_id_scratch;
        const auto el_begin = mesh.active_local_elements_begin();
        const auto el_end = mesh.active_local_elements_end();
        for (auto el_it = el_begin; el_it != el_end; ++el_it)
        {
            auto elem = *el_it;
            const auto& F_dof_indices = F_dof_map_cache.dof_indices(elem);
            const auto& X_dof_indices = X_dof_map_cache.dof_indices(elem);
            const auto& sig_jump_dof_indices = sig_jump_dof_map_cache.dof_indices(elem);

            for (unsigned int d = 0; d < NDIM; ++d)
            {
                F_rhs_e[d].resize(static_cast<int>(F_dof_indices[d].size()));
            }
            if (d_use_pressure_jump_conditions)
            {
                const auto& P_jump_dof_indices = P_jump_dof_map_cache->dof_indices(elem);
                P_jump_rhs_e.resize(static_cast<int>(P_jump_dof_indices[0].size()));
            }

            if (d_use_velocity_jump_conditions)
            {
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    const auto& DU_jump_dof_indices = DU_jump_dof_map_cache[d]->dof_indices(elem);
                    for (unsigned int k = 0; k < NDIM; ++k)
                    {
                        DU_jump_rhs_e[d][k].resize(static_cast<int>(DU_jump_dof_indices[k].size()));
                    }
                }
            }

            fe_X->reinit(elem);

            fe_interpolator.reinit(elem);
            fe_interpolator.collectDataForInterpolation(elem);
            fe_interpolator.interpolate(elem);

            if (d_use_pressure_jump_conditions || d_use_velocity_jump_conditions)
            {
                fe_jump->reinit(elem);
            }

            get_values_for_interpolation(x_node, *X_vec, X_dof_indices);
            get_values_for_interpolation(X_node, X0_vec, X_dof_indices);
            get_values_for_interpolation(sig_jump_node, *sig_jump_vec, sig_jump_dof_indices);
            const unsigned int n_qpoints = qrule->n_points();
            const size_t n_basis = phi_X.size();
            const size_t n_basis2 = phi_jump.size();
            for (unsigned int qp = 0; qp < n_qpoints; ++qp)
            {
                interpolate(X, qp, X_node, phi_X);
                interpolate(x, qp, x_node, phi_X);
                for (unsigned int k = 0; k < NDIM - 1; ++k)
                {
                    interpolate(dX_dxi[k], qp, X_node, *dphi_dxi_X[k]);
                    interpolate(dx_dxi[k], qp, x_node, *dphi_dxi_X[k]);
                }
                if (NDIM == 2)
                {
                    dX_dxi[1] = VectorValue<double>(0.0, 0.0, 1.0);
                    dx_dxi[1] = VectorValue<double>(0.0, 0.0, 1.0);
                }

                std::array<double, TENSOR_SIZE> sig_jump_vals;
                interpolate(&sig_jump_vals[0], qp, sig_jump_node, phi_jump);
                for (unsigned int d = 0; d < TENSOR_SIZE; ++d)
                {
                    std::pair<int, int> t_idx = IBTK::voigt_to_tensor_idx(d);
                    sig_jump(t_idx.first, t_idx.second) = sig_jump_vals[d];
                }

                // Construct unit vectors in the reference and current
                // configurations.
                N = dX_dxi[0].cross(dX_dxi[1]);
                const double dA = N.norm();
                N = N.unit();
                n = dx_dxi[0].cross(dx_dxi[1]);
                const double da = n.norm();
                n = n.unit();

                F.zero();

                if (d_lag_surface_pressure_fcn_data[part].fcn)
                {
                    // Compute the value of the pressure at the quadrature point
                    // and add the corresponding force to the right-hand-side
                    // vector.
                    double P = 0;
                    fe_interpolator.setInterpolatedDataPointers(surface_pressure_var_data,
                                                                surface_pressure_grad_var_data,
                                                                surface_pressure_fcn_system_idxs,
                                                                elem,
                                                                qp);
                    d_lag_surface_pressure_fcn_data[part].fcn(P,
                                                              n,
                                                              N,
                                                              FF,
                                                              x,
                                                              X,
                                                              elem,
                                                              /*side*/ 0,
                                                              surface_pressure_var_data,
                                                              surface_pressure_grad_var_data,
                                                              data_time,
                                                              d_lag_surface_pressure_fcn_data[part].ctx);
                    F -= P * n * da / dA;
                }

                if (d_lag_surface_force_fcn_data[part].fcn)
                {
                    // Compute the value of the surface force at the quadrature
                    // point and add the corresponding force to the
                    // right-hand-side vector.
                    fe_interpolator.setInterpolatedDataPointers(
                        surface_force_var_data, surface_force_grad_var_data, surface_force_fcn_system_idxs, elem, qp);
                    d_lag_surface_force_fcn_data[part].fcn(F_s,
                                                           n,
                                                           N,
                                                           FF,
                                                           x,
                                                           X,
                                                           elem,
                                                           /*side*/ 0,
                                                           surface_force_var_data,
                                                           surface_force_grad_var_data,
                                                           data_time,
                                                           d_lag_surface_force_fcn_data[part].ctx);
                    F += F_s;
                }

                double P_j = 0.0;
                if (d_use_extra_stress_jump_conditions)
                    P_j = F * n * dA / da + (sig_jump * n) * n;
                else
                    P_j = F * n * dA / da;
                VectorValue<double> du_sig_jump = ((sig_jump * n) * n) * n + sig_jump * n;
                if (!d_use_extra_stress_jump_conditions) du_sig_jump = 0.0;
                for (unsigned int i = 0; i < NDIM; ++i)
                    for (unsigned int k = 0; k < NDIM; ++k)
                        DU[i][k] =
                            -(dA / da) * (F(i) - F * n * n(i)) * n(k) - du_sig_jump(i) * n(k); // [Ux] , [Uy], [Uz]

                for (unsigned int d = 0; d < NDIM; ++d) F_integral(d) += F(d) * JxW[qp];

                // Remote the part of the force that has been already included in the jump

                if (d_use_pressure_jump_conditions && !d_use_velocity_jump_conditions)
                {
                    F -= (F * n) * n;
                }
                if (!d_use_pressure_jump_conditions && d_use_velocity_jump_conditions)
                {
                    F = (F * n) * n;
                }
                if (d_use_pressure_jump_conditions && d_use_velocity_jump_conditions)
                {
                    F = 0.0;
                }

                // Add the boundary forces to the right-hand-side vector.
                for (unsigned int k = 0; k < n_basis; ++k)
                {
                    F_qp = F * phi_X[k][qp] * JxW[qp];
                    for (unsigned int i = 0; i < NDIM; ++i)
                    {
                        F_rhs_e[i](k) += F_qp(i);
                    }
                }
                for (unsigned int k = 0; k < n_basis2; ++k)
                {
                    if (d_use_pressure_jump_conditions)
                    {
                        P_jump_rhs_e(k) += P_j * phi_jump[k][qp] * JxW[qp];
                    }
                    if (d_use_velocity_jump_conditions)
                    {
                        for (unsigned int i = 0; i < NDIM; ++i)
                        {
                            for (unsigned int j = 0; j < NDIM; ++j)
                            {
                                DU_jump_rhs_e[i][j](k) += DU[i][j] * phi_jump[k][qp] * JxW[qp];
                            }
                        }
                    }
                }
                if (d_use_pressure_jump_conditions)
                {
                    P_jump_rhs_integral += P_j * JxW[qp];
                    surface_area += JxW[qp];
                }
            }

            // Apply constraints (e.g., enforce periodic boundary conditions)
            // and add the elemental contributions to the global vector.
            for (unsigned int i = 0; i < NDIM; ++i)
            {
                copy_dof_ids_to_vector(i, F_dof_indices, dof_id_scratch);
                F_dof_map.constrain_element_vector(F_rhs_e[i], dof_id_scratch);
                F_rhs_vec->add_vector(F_rhs_e[i], dof_id_scratch);
                if (d_use_velocity_jump_conditions)
                {
                    const auto& DU_jump_dof_indices = DU_jump_dof_map_cache[i]->dof_indices(elem);
                    for (unsigned int k = 0; k < NDIM; ++k)
                    {
                        copy_dof_ids_to_vector(k, DU_jump_dof_indices, dof_id_scratch);
                        DU_jump_dof_map[i]->constrain_element_vector(DU_jump_rhs_e[i][k], dof_id_scratch);
                        DU_jump_rhs_vec[i]->add_vector(DU_jump_rhs_e[i][k], dof_id_scratch);
                    }
                }
            }
            if (d_use_pressure_jump_conditions)
            {
                const auto& P_jump_dof_indices = P_jump_dof_map_cache->dof_indices(elem);
                copy_dof_ids_to_vector(0, P_jump_dof_indices, dof_id_scratch);
                P_jump_dof_map->constrain_element_vector(P_jump_rhs_e, dof_id_scratch);
                P_jump_rhs_vec->add_vector(P_jump_rhs_e, dof_id_scratch);
            }
        }

        SAMRAI_MPI::sumReduction(&F_integral(0), NDIM);

        // Solve for F.
        F_rhs_vec->close();
        d_fe_data_managers[part]->computeL2Projection(
            *F_vec, *F_rhs_vec, FORCE_SYSTEM_NAME, d_default_interp_spec.use_consistent_mass_matrix);
        F_vec->close();
        if (d_use_pressure_jump_conditions)
        {
            P_jump_rhs_vec->close();
            d_fe_data_managers[part]->computeL2Projection(*P_jump_vec,
                                                          *P_jump_rhs_vec,
                                                          PRESSURE_JUMP_SYSTEM_NAME,
                                                          d_default_interp_spec.use_consistent_mass_matrix);
            P_jump_rhs_integral = SAMRAI_MPI::sumReduction(P_jump_rhs_integral);
            surface_area = SAMRAI_MPI::sumReduction(surface_area);
            if (d_normalize_pressure_jump[part]) P_jump_vec->add(-P_jump_rhs_integral / surface_area);
            P_jump_vec->close();
        }
        if (d_use_velocity_jump_conditions)
        {
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                DU_jump_rhs_vec[d]->close();
                d_fe_data_managers[part]->computeL2Projection(*DU_jump_vec[d],
                                                              *DU_jump_rhs_vec[d],
                                                              VELOCITY_JUMP_SYSTEM_NAME[d],
                                                              d_default_interp_spec.use_consistent_mass_matrix);
                DU_jump_vec[d]->close();
            }
        }
    }
    IBAMR_TIMER_STOP(t_compute_lagrangian_force);
    return;
} // computeLagrangianForce

void
CFIIMethod::spreadForce(const int f_data_idx,
                        RobinPhysBdryPatchStrategy* f_phys_bdry_op,
                        const std::vector<Pointer<RefineSchedule<NDIM> > >& /*f_prolongation_scheds*/,
                        const double data_time)
{
    TBOX_ASSERT(MathUtilities<double>::equalEps(data_time, d_half_time));
    IBAMR_TIMER_START(t_spread_force);

    std::vector<std::vector<libMesh::PetscVector<double>*> > vec_collection_update = {
        d_X_IB_ghost_vecs, d_X_half_vecs, d_F_IB_ghost_vecs, d_F_half_vecs
    };

    if (d_use_pressure_jump_conditions)
    {
        vec_collection_update.push_back(d_P_jump_IB_ghost_vecs);
        vec_collection_update.push_back(d_P_jump_half_vecs);
    }

    if (d_use_velocity_jump_conditions)
    {
        for (unsigned part = 0; part < d_num_parts; ++part)
        {
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                vec_collection_update.push_back({ d_DU_jump_half_vecs[part][d], d_DU_jump_IB_ghost_vecs[part][d] });
            }
        }
    }

    batch_vec_ghost_update(vec_collection_update, INSERT_VALUES, SCATTER_FORWARD);

    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        PetscVector<double>* X_vec = d_X_half_vecs[part];
        PetscVector<double>* X_ghost_vec = d_X_IB_ghost_vecs[part];
        PetscVector<double>* F_vec = d_F_half_vecs[part];
        PetscVector<double>* F_ghost_vec = d_F_IB_ghost_vecs[part];
        X_vec->localize(*X_ghost_vec);
        F_vec->localize(*F_ghost_vec);
        d_fe_data_managers[part]->spread(
            f_data_idx, *F_ghost_vec, *X_ghost_vec, FORCE_SYSTEM_NAME, f_phys_bdry_op, data_time);
        PetscVector<double>* P_jump_vec;
        PetscVector<double>* P_jump_ghost_vec = NULL;
        std::array<PetscVector<double>*, NDIM> DU_jump_ghost_vec;
        std::array<PetscVector<double>*, NDIM> DU_jump_vec;
        if (d_use_pressure_jump_conditions)
        {
            P_jump_vec = d_P_jump_half_vecs[part];
            P_jump_ghost_vec = d_P_jump_IB_ghost_vecs[part];
            P_jump_vec->localize(*P_jump_ghost_vec);
        }
        if (d_use_velocity_jump_conditions)
        {
            for (auto d = 0; d < NDIM; ++d)
            {
                DU_jump_ghost_vec[d] = d_DU_jump_IB_ghost_vecs[part][d];
                DU_jump_vec[d] = d_DU_jump_half_vecs[part][d];
                DU_jump_vec[d]->localize(*DU_jump_ghost_vec[d]);
            }
        }

        if (d_use_pressure_jump_conditions || d_use_velocity_jump_conditions)
        {
            imposeJumpConditions(f_data_idx, *P_jump_ghost_vec, DU_jump_ghost_vec, *X_ghost_vec, data_time, part);
        }
        if (d_use_velocity_jump_conditions)
        {
            for (unsigned int d = 0; d < NDIM; ++d) d_DU_jump_IB_ghost_vecs[part][d]->close();
        }
        if (d_use_pressure_jump_conditions)
        {
            d_P_jump_IB_ghost_vecs[part]->close();
        }

        d_F_IB_ghost_vecs[part]->close();
        d_X_IB_ghost_vecs[part]->close();
    }
    IBAMR_TIMER_STOP(t_spread_force);
    return;
} // spreadForce

void
CFIIMethod::initializeFEEquationSystems()
{
    if (d_fe_equation_systems_initialized) return;

    const bool from_restart = RestartManager::getManager()->isFromRestart();

    // Create the FE data managers that manage mappings between the FE mesh
    // parts and the Cartesian grid.
    d_equation_systems.resize(d_num_parts, nullptr);
    d_fe_data_managers.resize(d_num_parts, nullptr);
    IntVector<NDIM> min_ghost_width(0);
    if (!d_eulerian_data_cache) d_eulerian_data_cache.reset(new SAMRAIDataCache());
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        // Create FE equation systems objects and corresponding variables.
        d_equation_systems[part] = new EquationSystems(*d_meshes[part]);
        EquationSystems* equation_systems = d_equation_systems[part];
        auto fe_data = std::make_shared<FEData>(
            d_object_name + "::FEdata::" + std::to_string(part), *equation_systems, d_registered_for_restart);

        // Create FE data managers.
        const std::string manager_name = "CFIIMethod FEDataManager::" + std::to_string(part);
        Pointer<InputDatabase> fe_data_manager_db(new InputDatabase(manager_name + "::input_db"));

        d_fe_data_managers[part] = FEDataManager::getManager(fe_data,
                                                             manager_name,
                                                             fe_data_manager_db,
                                                             d_max_level_number + 1,
                                                             d_interp_spec[part],
                                                             d_spread_spec[part],
                                                             d_default_workload_spec,
                                                             min_ghost_width,
                                                             d_eulerian_data_cache);
        d_ghosts = IntVector<NDIM>::max(d_ghosts, d_fe_data_managers[part]->getGhostCellWidth());
        d_fe_data_managers[part]->COORDINATES_SYSTEM_NAME = COORDS_SYSTEM_NAME;
        if (from_restart)
        {
            const std::string& file_name = libmesh_restart_file_name(
                d_libmesh_restart_read_dir, d_libmesh_restart_restore_number, part, d_libmesh_restart_file_extension);
            const XdrMODE xdr_mode = (d_libmesh_restart_file_extension == "xdr" ? DECODE : READ);
            const int read_mode =
                EquationSystems::READ_HEADER | EquationSystems::READ_DATA | EquationSystems::READ_ADDITIONAL_DATA;
            equation_systems->read(file_name, xdr_mode, read_mode, /*partition_agnostic*/ true);
        }
        else
        {
            // vector FE systems:
            std::vector<std::string> vector_system_names{
                COORDS_SYSTEM_NAME,          COORD_MAPPING_SYSTEM_NAME,       VELOCITY_SYSTEM_NAME,
                NORMAL_VELOCITY_SYSTEM_NAME, TANGENTIAL_VELOCITY_SYSTEM_NAME, FORCE_SYSTEM_NAME
            };
            std::vector<std::string> vector_variable_prefixes{ "X", "dX", "U", "U_n", "U_t", "F" };
            std::vector<libMesh::FEFamily> vector_fe_family(vector_system_names.size(), d_fe_family[part]);

            if (d_use_velocity_jump_conditions)
            {
                const auto jump_family =
                    d_use_discon_elem_for_jumps[part] ? d_velocity_jump_fe_family : d_fe_family[part];
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    vector_system_names.push_back(VELOCITY_JUMP_SYSTEM_NAME[d]);
                    vector_variable_prefixes.push_back("DU_jump_" + std::to_string(d));
                    vector_fe_family.push_back(jump_family);
                }

                vector_system_names.push_back(WSS_IN_SYSTEM_NAME);
                vector_variable_prefixes.push_back("WSS_in");
                vector_fe_family.push_back(jump_family);

                vector_system_names.push_back(WSS_OUT_SYSTEM_NAME);
                vector_variable_prefixes.push_back("WSS_out");
                vector_fe_family.push_back(jump_family);

                if (d_use_pressure_jump_conditions)
                {
                    vector_system_names.push_back(TAU_IN_SYSTEM_NAME);
                    vector_variable_prefixes.push_back("TAU_IN");
                    vector_fe_family.push_back(jump_family);

                    vector_system_names.push_back(TAU_OUT_SYSTEM_NAME);
                    vector_variable_prefixes.push_back("TAU_OUT");
                    vector_fe_family.push_back(jump_family);
                }
            }

            for (std::size_t i = 0; i < vector_system_names.size(); ++i)
            {
                auto& system = equation_systems->add_system<System>(vector_system_names[i]);
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    system.add_variable(
                        vector_variable_prefixes[i] + "_" + std::to_string(d), d_fe_order[part], vector_fe_family[i]);
                }
            }

            equation_systems->get_system(COORDS_SYSTEM_NAME)
                .add_vector("INITIAL_COORDINATES", /*projections*/ true, GHOSTED);

            // scalar FE systems:
            if (d_use_pressure_jump_conditions)
            {
                System& P_jump_system = equation_systems->add_system<System>(PRESSURE_JUMP_SYSTEM_NAME);
                System& P_in_system = equation_systems->add_system<System>(PRESSURE_IN_SYSTEM_NAME);
                System& P_out_system = equation_systems->add_system<System>(PRESSURE_OUT_SYSTEM_NAME);
                if (d_use_discon_elem_for_jumps[part])
                {
                    P_jump_system.add_variable("P_jump_", d_fe_order[part], d_pressure_jump_fe_family);
                    P_in_system.add_variable("P_in_", d_fe_order[part], d_pressure_jump_fe_family);
                    P_out_system.add_variable("P_out_", d_fe_order[part], d_pressure_jump_fe_family);
                }
                else
                {
                    P_jump_system.add_variable("P_jump_", d_fe_order[part], d_fe_family[part]);
                    P_in_system.add_variable("P_in_", d_fe_order[part], d_fe_family[part]);
                    P_out_system.add_variable("P_out_", d_fe_order[part], d_fe_family[part]);
                }
            }

            // Extra stress FE system:
            System& stress_jump_system = equation_systems->add_system<System>(EXTRA_STRESS_JUMP_SYSTEM_NAME);
            for (int d = 0; d < TENSOR_SIZE; ++d)
                stress_jump_system.add_variable(
                    "Sig_jump_" + std::to_string(d), d_fe_order[part], d_sig_jump_fe_family);

            System& stress_in_system = equation_systems->add_system<System>(EXTRA_STRESS_IN_SYSTEM_NAME);
            for (int d = 0; d < TENSOR_SIZE; ++d)
                stress_in_system.add_variable("Sig_in_" + std::to_string(d), d_fe_order[part], d_sig_jump_fe_family);

            System& stress_out_system = equation_systems->add_system<System>(EXTRA_STRESS_OUT_SYSTEM_NAME);
            for (int d = 0; d < TENSOR_SIZE; ++d)
                stress_out_system.add_variable("Sig_out_" + std::to_string(d), d_fe_order[part], d_sig_jump_fe_family);
        }

        const std::vector<std::string> system_names{ COORDS_SYSTEM_NAME };
        const std::vector<std::string> vector_names{ "INITIAL_COORDINATES" };
        IBTK::setup_system_vectors(equation_systems, system_names, vector_names, from_restart);
    }
    d_fe_equation_systems_initialized = true;
    return;
} // initializeFEEquationSystems

void
CFIIMethod::initializeFEData()
{
    if (d_fe_data_initialized) return;
    initializeFEEquationSystems();
    const bool from_restart = RestartManager::getManager()->isFromRestart();
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        // Initialize FE equation systems.
        EquationSystems* equation_systems = d_equation_systems[part];
        if (from_restart)
        {
            equation_systems->reinit();
        }
        else
        {
            equation_systems->init();
            initializeCoordinates(part);
            initializeVelocity(part);
        }
        updateCoordinateMapping(part);

        // Assemble systems.
        auto& X_system = equation_systems->get_system<System>(COORDS_SYSTEM_NAME);
        auto& dX_system = equation_systems->get_system<System>(COORD_MAPPING_SYSTEM_NAME);
        auto& U_system = equation_systems->get_system<System>(VELOCITY_SYSTEM_NAME);
        auto& U_n_system = equation_systems->get_system<System>(NORMAL_VELOCITY_SYSTEM_NAME);
        auto& U_t_system = equation_systems->get_system<System>(TANGENTIAL_VELOCITY_SYSTEM_NAME);
        auto& F_system = equation_systems->get_system<System>(FORCE_SYSTEM_NAME);

        X_system.assemble_before_solve = false;
        X_system.assemble();

        dX_system.assemble_before_solve = false;
        dX_system.assemble();

        U_system.assemble_before_solve = false;
        U_system.assemble();

        U_n_system.assemble_before_solve = false;
        U_n_system.assemble();

        U_t_system.assemble_before_solve = false;
        U_t_system.assemble();

        F_system.assemble_before_solve = false;
        F_system.assemble();

        if (d_use_pressure_jump_conditions)
        {
            System& P_jump_system = equation_systems->get_system<System>(PRESSURE_JUMP_SYSTEM_NAME);
            P_jump_system.assemble_before_solve = false;
            P_jump_system.assemble();

            System& P_in_system = equation_systems->get_system<System>(PRESSURE_IN_SYSTEM_NAME);
            P_in_system.assemble_before_solve = false;
            P_in_system.assemble();

            System& P_out_system = equation_systems->get_system<System>(PRESSURE_OUT_SYSTEM_NAME);
            P_out_system.assemble_before_solve = false;
            P_out_system.assemble();
        }

        if (d_use_velocity_jump_conditions)
        {
            std::array<System*, NDIM> DU_jump_system;
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                DU_jump_system[d] = &equation_systems->get_system<System>(VELOCITY_JUMP_SYSTEM_NAME[d]);
                DU_jump_system[d]->assemble_before_solve = false;
                DU_jump_system[d]->assemble();
            }
            System& WSS_in_system = equation_systems->get_system<System>(WSS_IN_SYSTEM_NAME);
            WSS_in_system.assemble_before_solve = false;
            WSS_in_system.assemble();

            System& WSS_out_system = equation_systems->get_system<System>(WSS_OUT_SYSTEM_NAME);
            WSS_out_system.assemble_before_solve = false;
            WSS_out_system.assemble();
        }

        System& sig_jump_sys = equation_systems->get_system<System>(EXTRA_STRESS_JUMP_SYSTEM_NAME);
        sig_jump_sys.assemble_before_solve = false;
        sig_jump_sys.assemble();

        System& sig_in_sys = equation_systems->get_system<System>(EXTRA_STRESS_IN_SYSTEM_NAME);
        sig_in_sys.assemble_before_solve = false;
        sig_in_sys.assemble();

        System& sig_out_sys = equation_systems->get_system<System>(EXTRA_STRESS_OUT_SYSTEM_NAME);
        sig_out_sys.assemble_before_solve = false;
        sig_out_sys.assemble();

        if (d_use_pressure_jump_conditions && d_use_velocity_jump_conditions)
        {
            System& TAU_in_system = equation_systems->get_system<System>(TAU_IN_SYSTEM_NAME);
            TAU_in_system.assemble_before_solve = false;
            TAU_in_system.assemble();

            System& TAU_out_system = equation_systems->get_system<System>(TAU_OUT_SYSTEM_NAME);
            TAU_out_system.assemble_before_solve = false;
            TAU_out_system.assemble();
        }
    }
    d_fe_data_initialized = true;
    return;
} // initializeFEData

void
CFIIMethod::registerEulerianVariables()
{
    IIMethod::registerEulerianVariables();
    d_sig_var = new CellVariable<NDIM, double>(d_object_name + "::SIG", TENSOR_SIZE);
    registerVariable(d_sig_scratch_idx, d_sig_var, d_ghosts);
    return;
} // registerEulerianVariables

/////////////////////////////// PROTECTED ////////////////////////////////////

namespace
{
struct IndexOrder
{
    inline bool operator()(const SAMRAI::hier::Index<NDIM>& lhs, const SAMRAI::hier::Index<NDIM>& rhs) const
    {
        return (lhs(0) < rhs(0)
#if (NDIM > 1)
                || (lhs(0) == rhs(0) && lhs(1) < rhs(1))
#if (NDIM > 2)
                || (lhs(0) == rhs(0) && lhs(1) == rhs(1) && lhs(2) < rhs(2))
#endif
#endif
        );
    }
};
} // namespace

void
CFIIMethod::imposeJumpConditions(const int f_data_idx,
                                 PetscVector<double>& P_jump_ghost_vec,
                                 std::array<PetscVector<double>*, NDIM>& DU_jump_ghost_vec,
                                 PetscVector<double>& X_ghost_vec,
                                 const double /*data_time*/,
                                 const unsigned int part)
{
    // Extract the mesh.
    EquationSystems* equation_systems = d_fe_data_managers[part]->getEquationSystems();
    const MeshBase& mesh = equation_systems->get_mesh();
    const unsigned int dim = mesh.mesh_dimension();

    // Extract the FE systems and DOF maps, and setup the FE object
    System& X_system = equation_systems->get_system(COORDS_SYSTEM_NAME);
    const DofMap& X_dof_map = X_system.get_dof_map();
    FEDataManager::SystemDofMapCache& X_dof_map_cache = *d_fe_data_managers[part]->getDofMapCache(COORDS_SYSTEM_NAME);
    FEType X_fe_type = X_dof_map.variable_type(0);
    for (unsigned int d = 0; d < NDIM; ++d) TBOX_ASSERT(X_dof_map.variable_type(d) == X_fe_type);

    System* P_jump_system;
    const DofMap* P_jump_dof_map;
    FEDataManager::SystemDofMapCache* P_jump_dof_map_cache = NULL;
    FEType P_jump_fe_type;
    if (d_use_pressure_jump_conditions)
    {
        P_jump_system = &equation_systems->get_system(PRESSURE_JUMP_SYSTEM_NAME);
        P_jump_dof_map = &P_jump_system->get_dof_map();
        P_jump_dof_map_cache = d_fe_data_managers[part]->getDofMapCache(PRESSURE_JUMP_SYSTEM_NAME);
        P_jump_fe_type = P_jump_dof_map->variable_type(0);
    }

    std::array<DofMap*, NDIM> DU_jump_dof_map;
    std::array<FEDataManager::SystemDofMapCache*, NDIM> DU_jump_dof_map_cache;
    std::array<System*, NDIM> DU_jump_system;
    FEType DU_jump_fe_type;
    if (d_use_velocity_jump_conditions)
    {
        for (unsigned int i = 0; i < NDIM; ++i)
        {
            DU_jump_system[i] = &equation_systems->get_system(VELOCITY_JUMP_SYSTEM_NAME[i]);
            DU_jump_dof_map_cache[i] = d_fe_data_managers[part]->getDofMapCache(VELOCITY_JUMP_SYSTEM_NAME[i]);
            DU_jump_dof_map[i] = &DU_jump_system[i]->get_dof_map();
            DU_jump_fe_type = DU_jump_dof_map[i]->variable_type(0);
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                TBOX_ASSERT(DU_jump_dof_map[i]->variable_type(d) == DU_jump_fe_type);
            }
        }
    }

    System& sig_jump_system = d_equation_systems[part]->get_system(EXTRA_STRESS_JUMP_SYSTEM_NAME);
    NumericVector<double>* sig_jump_vec = sig_jump_system.current_local_solution.get();
    FEDataManager::SystemDofMapCache* sig_jump_dof_map_cache =
        d_fe_data_managers[part]->getDofMapCache(EXTRA_STRESS_JUMP_SYSTEM_NAME);

    FEType fe_type = X_fe_type;
    std::unique_ptr<FEBase> fe_X = FEBase::build(dim, fe_type);
    std::array<const std::vector<std::vector<double> >*, NDIM - 1> dphi_dxi;
    dphi_dxi[0] = &fe_X->get_dphidxi();
    if (NDIM > 2) dphi_dxi[1] = &fe_X->get_dphideta();

    std::unique_ptr<FEBase> fe_P_jump = FEBase::build(dim, P_jump_fe_type);
    const std::vector<std::vector<double> >& phi_P_jump = fe_P_jump->get_phi();

    // Loop over the patches to impose jump conditions on the Eulerian grid.
    const std::vector<std::vector<Elem*> >& active_patch_element_map =
        d_fe_data_managers[part]->getActivePatchElementMap();
    const int level_num = d_fe_data_managers[part]->getFinestPatchLevelNumber();
    boost::multi_array<double, 1> P_jump_node;
    boost::multi_array<double, 2> x_node;
    std::array<boost::multi_array<double, 2>, NDIM> DU_jump_node;
    std::array<boost::multi_array<double, 1>, TENSOR_SIZE> sig_jump_node;
    std::array<VectorValue<double>, 2> dx_dxi;
    VectorValue<double> n, jn;
    std::vector<libMesh::Point> X_node_cache, x_node_cache;
    IBTK::Point x_min, x_max;
    Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(level_num);
    const IntVector<NDIM>& ratio = level->getRatio();
    const Pointer<CartesianGridGeometry<NDIM> > grid_geom = level->getGridGeometry();
    int local_patch_num = 0;
    for (PatchLevel<NDIM>::Iterator p(level); p; p++, ++local_patch_num)
    {
        // The relevant collection of elements.
        const std::vector<Elem*>& patch_elems = active_patch_element_map[local_patch_num];
        const size_t num_active_patch_elems = patch_elems.size();
        if (num_active_patch_elems == 0) continue;

        const Pointer<Patch<NDIM> > patch = level->getPatch(p());
        Pointer<SideData<NDIM, double> > f_data = patch->getPatchData(f_data_idx);
        const Box<NDIM>& patch_box = patch->getBox();
        const CellIndex<NDIM>& patch_lower = patch_box.lower();
        std::array<Box<NDIM>, NDIM> side_ghost_boxes;
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            side_ghost_boxes[d] = SideGeometry<NDIM>::toSideBox(f_data->getGhostBox(), d);
        }

        Box<NDIM> side_boxes[NDIM];
        for (int d = 0; d < NDIM; ++d)
        {
            side_boxes[d] = SideGeometry<NDIM>::toSideBox(patch_box, d);
        }

        const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
        const double* const x_lower = patch_geom->getXLower();
        const double* const dx = patch_geom->getDx();

        std::array<std::map<hier::Index<NDIM>, std::vector<libMesh::Point>, IndexOrder>, NDIM> intersection_points,
            intersection_ref_coords;
        std::array<std::map<hier::Index<NDIM>, std::vector<VectorValue<double> >, IndexOrder>, NDIM>
            intersection_normals;

        std::array<std::map<hier::Index<NDIM>, std::vector<libMesh::Point>, IndexOrder>, NDIM> intersection_u_points,
            intersection_u_ref_coords;
        std::array<std::map<hier::Index<NDIM>, std::vector<VectorValue<double> >, IndexOrder>, NDIM>
            intersection_u_normals;

        std::array<std::array<std::map<hier::Index<NDIM>, std::vector<libMesh::Point>, IndexOrder>, NDIM>, NDIM>
            intersectionSide_u_points, intersectionSide_u_ref_coords;
        std::array<std::array<std::map<hier::Index<NDIM>, std::vector<VectorValue<double> >, IndexOrder>, NDIM>, NDIM>
            intersectionSide_u_normals;

        // Loop over the elements.
        std::vector<libMesh::dof_id_type> dof_id_scratch;
        for (size_t e_idx = 0; e_idx < num_active_patch_elems; ++e_idx)
        {
            Elem* const elem = patch_elems[e_idx];
            const auto& X_dof_indices = X_dof_map_cache.dof_indices(elem);
            get_values_for_interpolation(x_node, X_ghost_vec, X_dof_indices);
            if (d_use_pressure_jump_conditions)
            {
                const auto& P_jump_dof_indices = P_jump_dof_map_cache->dof_indices(elem);
                copy_dof_ids_to_vector(0, P_jump_dof_indices, dof_id_scratch);
                get_values_for_interpolation(P_jump_node, P_jump_ghost_vec, dof_id_scratch);
            }
            if (d_use_velocity_jump_conditions)
            {
                for (unsigned int axis = 0; axis < NDIM; ++axis)
                {
                    const auto& DU_jump_dof_indices = DU_jump_dof_map_cache[axis]->dof_indices(elem);
                    get_values_for_interpolation(DU_jump_node[axis], *DU_jump_ghost_vec[axis], DU_jump_dof_indices);
                }
            }
            if (d_use_extra_stress_jump_conditions)
            {
                for (unsigned int d = 0; d < TENSOR_SIZE; ++d)
                {
                    std::vector<dof_id_type> sig_dof_indices;
                    sig_jump_dof_map_cache->dof_indices(elem, sig_dof_indices, d);
                    get_values_for_interpolation(sig_jump_node[d], *sig_jump_vec, sig_dof_indices);
                }
            }

            // Cache the nodal and physical coordinates of the side element,
            // determine the bounding box of the current configuration of the
            // element, and set the nodal coordinates to correspond to the
            // physical coordinates.
            const unsigned int n_nodes = elem->n_nodes();
            X_node_cache.resize(n_nodes);
            x_node_cache.resize(n_nodes);
            x_min = IBTK::Point::Constant(std::numeric_limits<double>::max());
            x_max = IBTK::Point::Constant(-std::numeric_limits<double>::max());
            for (unsigned int k = 0; k < n_nodes; ++k)
            {
                X_node_cache[k] = elem->point(k);
                libMesh::Point& x = x_node_cache[k];
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    x(d) = x_node[k][d];
                }
                if (d_perturb_fe_mesh_nodes)
                {
                    // Perturb the mesh configuration to keep the FE mesh nodes
                    // away from cell edges, nodes, and centers.
                    //
                    // This implies that we only have to deal with multiple
                    // intersections along element edges, and not at element
                    // nodes.
                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        const int i_s = boost::math::iround(((x(d) - x_lower[d]) / dx[d]) - 0.5) + patch_lower[d];
                        for (int shift = 0; shift <= 2; ++shift)
                        {
                            const double x_s =
                                x_lower[d] + dx[d] * (static_cast<double>(i_s - patch_lower[d]) + 0.5 * shift);
                            const double tol = 1.0e-4 * dx[d];
                            if (x(d) <= x_s) x(d) = std::min(x_s - tol, x(d));
                            if (x(d) >= x_s) x(d) = std::max(x_s + tol, x(d));
                        }
                    }
                }
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    x_min[d] = std::min(x_min[d], x(d));
                    x_max[d] = std::max(x_max[d], x(d));
                }
                elem->point(k) = x;
            }
            Box<NDIM> box(IndexUtilities::getCellIndex(&x_min[0], grid_geom, ratio),
                          IndexUtilities::getCellIndex(&x_max[0], grid_geom, ratio));
            box.grow(IntVector<NDIM>(1));
            box = box * patch_box;

            // Loop over coordinate directions and look for intersections with
            // the background fluid grid.
            // The axis variable determines the direction of the derivative for which we apply corrections.
            // E.g. axis = 0 --> need to apply corrections for dp/dx, d^2u/dx^2, d^2v/dx^2, d^2w/dx^2.
            for (unsigned int axis = 0; axis < NDIM; ++axis)
            {
                Box<NDIM> extended_box = patch_box;
                extended_box.grow(IntVector<NDIM>(1));
                Box<NDIM> extended_side_box = patch_box;
                extended_side_box.grow(IntVector<NDIM>(2));
                if (patch_geom->getTouchesRegularBoundary(axis, 1)) extended_box.upper(axis) += 1;

                // These are never used???
                Box<NDIM> side_u_boxes[NDIM];
                for (int d = 0; d < NDIM; ++d)
                {
                    side_u_boxes[d] = SideGeometry<NDIM>::toSideBox(extended_side_box, d);
                }

                // Setup a unit vector pointing in the coordinate direction of
                // interest.
                VectorValue<double> q;
                q(axis) = 1.0;

                // Loop over the relevant range of indices.
                // We search in the axis direction for the corrections that we need to apply.
                // If we are searching for corrections when axis = 0 (so in the x-direction), we can check all the
                // x-indices at once because the line in one cell is the same as the next. So we only loop over y (and
                // z) indices.
                Box<NDIM> axis_box = box;
                axis_box.lower(axis) = 0;
                axis_box.upper(axis) = 0;

                // Note sure why we need a double array. We only index SideDim with the first argument being [axis]
                // The second argument contains the other two axes.
                unsigned int SideDim[NDIM][NDIM - 1];
                for (unsigned int d = 0; d < NDIM; ++d)
                    for (unsigned int l = 0; l < NDIM - 1; ++l) SideDim[d][l] = (d + l + 1) % NDIM;

                for (BoxIterator<NDIM> b(axis_box); b; b++)
                {
                    const hier::Index<NDIM>& i_c = b();
                    libMesh::Point r;
                    std::array<libMesh::Point, NDIM - 1> rs;

                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        // Grab the point on the cell side.
                        // Need this for pressure derivative and normal u derivative.
                        // If axis = 0, this is needed for dp/dx and d^2u/dx^2.
                        // Note that the component in the axis direction is not needed.
                        r(d) = (d == axis ? 0.0 :
                                            x_lower[d] + dx[d] * (static_cast<double>(i_c(d) - patch_lower[d]) + 0.5));

                        for (unsigned int l = 0; l < NDIM - 1; ++l)
                        {
                            // Grab the point close to which tangential components would live
                            // If axis = 0, this is needed for d^2v/dx^2 and d^2w/dx^2.
                            // Note that the component in the axis direction is not needed.
                            rs[l](d) =
                                (d == axis ? 0.0 :
                                             d == SideDim[axis][l] ?
                                             x_lower[d] + dx[d] * (static_cast<double>(i_c(d) - patch_lower[d])) :
                                             x_lower[d] + dx[d] * (static_cast<double>(i_c(d) - patch_lower[d]) + 0.5));
                        }
                    }

                    std::vector<std::pair<double, libMesh::Point> > intersections;
                    std::array<std::vector<std::pair<double, libMesh::Point> >, NDIM - 1> intersectionsSide;

                    static const double tolerance = sqrt(std::numeric_limits<double>::epsilon());

#if (NDIM == 2)
                    // Search for intersections with the point r and direction q (axis) with the current element.
                    // If axis = 0, this is used to compute dp/dx and d^2u/dx^2
                    intersect_line_with_edge(intersections, static_cast<Edge*>(elem), r, q, tolerance);
#endif
#if (NDIM == 3)
                    intersect_line_with_face(intersections, static_cast<Face*>(elem), r, q, tolerance);
#endif
                    for (unsigned int l = 0; l < NDIM - 1; ++l)
                    {
#if (NDIM == 2)
                        // Search for intersections with the point rs[l] and direction q (axis) with the current
                        // element. If axis = 0, this is used to compute d^2v/dx^2 and d^2w/dx^2.
                        intersect_line_with_edge(intersectionsSide[l], static_cast<Edge*>(elem), rs[l], q, tolerance);
#endif
#if (NDIM == 3)
                        intersect_line_with_face(intersectionsSide[l], static_cast<Face*>(elem), rs[l], q, tolerance);
#endif
                    }

                    if (d_use_pressure_jump_conditions)
                    {
                        for (unsigned int k = 0; k < intersections.size(); ++k)
                        {
                            // Note that r[axis] has the correct component. The other component is found from the
                            // location of the intersection.
                            const libMesh::Point x = r + intersections[k].first * q;
                            const libMesh::Point& xi = intersections[k].second;
                            SideIndex<NDIM> i_s(i_c, axis, 0);
                            i_s(axis) = boost::math::iround((x(axis) - x_lower[axis]) / dx[axis]) + patch_lower[axis];
                            if (extended_box.contains(i_s))
                            {
                                std::vector<libMesh::Point> ref_coords(1, xi);
                                fe_X->reinit(elem, &ref_coords);
                                fe_P_jump->reinit(elem, &ref_coords);
                                for (unsigned int l = 0; l < NDIM - 1; ++l)
                                {
                                    interpolate(dx_dxi[l], 0, x_node, *dphi_dxi[l]);
                                }
                                if (NDIM == 2)
                                {
                                    dx_dxi[1] = VectorValue<double>(0.0, 0.0, 1.0);
                                }
                                n = (dx_dxi[0].cross(dx_dxi[1])).unit();

                                // Make sure we haven't already found this
                                // intersection.
                                //
                                // (Because we are doing this in floating point
                                // arithmetic, we can't even count on the
                                // intersection being assigned to the same index!)
                                bool found_same_intersection_point = false;
                                for (int shift = -1; shift <= 1; ++shift)
                                {
                                    SideIndex<NDIM> i_s_prime = i_s;
                                    i_s_prime(axis) += shift;
                                    const std::vector<libMesh::Point>& candidate_coords =
                                        intersection_points[axis][i_s_prime];
                                    const std::vector<libMesh::Point>& candidate_ref_coords =
                                        intersection_ref_coords[axis][i_s_prime];
                                    const std::vector<VectorValue<double> >& candidate_normals =
                                        intersection_normals[axis][i_s_prime];

                                    found_same_intersection_point =
                                        checkDoubleCountingIntersection(axis,
                                                                        dx,
                                                                        n,
                                                                        x,
                                                                        xi,
                                                                        i_s,
                                                                        i_s_prime,
                                                                        candidate_coords,
                                                                        candidate_ref_coords,
                                                                        candidate_normals);
                                    if (found_same_intersection_point) break;
                                }

                                if (!found_same_intersection_point)
                                {
                                    // Evaluate the jump conditions and apply them
                                    // to the Eulerian grid.
                                    if (side_ghost_boxes[axis].contains(i_s))
                                    {
                                        // There are contributions from stress and pressure.
                                        const double C_p = interpolate(0, P_jump_node, phi_P_jump);
                                        const double C_sig = interpolate(0, sig_jump_node[axis], phi_P_jump);
                                        const double sgn = n(axis) > 0.0 ? 1.0 : n(axis) < 0.0 ? -1.0 : 0.0;
                                        (*f_data)(i_s) += sgn * ((C_p + C_sig) / dx[axis]);
                                    }

                                    // Keep track of the positions where we have
                                    // imposed jump conditions.
                                    intersection_points[axis][i_s].push_back(x);
                                    intersection_ref_coords[axis][i_s].push_back(xi);
                                    intersection_normals[axis][i_s].push_back(n);
                                }
                            }
                        }
                    }

                    if (d_use_velocity_jump_conditions)
                    {
                        for (unsigned int k = 0; k < intersections.size(); ++k)
                        {
                            // Note that r[axis] has the correct component. The other component is found from the
                            // location of the intersection.
                            libMesh::Point xu = r + intersections[k].first * q;
                            const libMesh::Point& xui = intersections[k].second;
                            // Note there are two stencils that contain this intersection.
                            SideIndex<NDIM> i_s_um(i_c, axis, 0);
                            hier::Index<NDIM> i_c_neighbor = i_c;
                            i_c_neighbor(axis) += 1;
                            // This index is created from i_c_neighbor, which is a shift of i_c in the axis direction,
                            // but we overwrite that value anyway?
                            SideIndex<NDIM> i_s_up(i_c_neighbor, axis, 0);
                            i_s_up(axis) =
                                boost::math::iround((xu(axis) - x_lower[axis]) / dx[axis] + 0.5) + patch_lower[axis];
                            i_s_um(axis) =
                                boost::math::iround((xu(axis) - x_lower[axis]) / dx[axis] - 0.5) + patch_lower[axis];

                            if (extended_box.contains(i_s_up) && extended_box.contains(i_s_um))
                            {
                                std::vector<libMesh::Point> ref_coords(1, xui);
                                fe_X->reinit(elem, &ref_coords);
                                fe_P_jump->reinit(elem, &ref_coords);
                                for (unsigned int l = 0; l < NDIM - 1; ++l)
                                {
                                    interpolate(dx_dxi[l], 0, x_node, *dphi_dxi[l]);
                                }
                                if (NDIM == 2)
                                {
                                    dx_dxi[1] = VectorValue<double>(0.0, 0.0, 1.0);
                                }
                                n = (dx_dxi[0].cross(dx_dxi[1])).unit();

                                bool found_same_intersection_point = false;

                                for (int shift = -1; shift <= 1; ++shift)
                                {
                                    SideIndex<NDIM> i_s_prime = i_s_um;
                                    i_s_prime(axis) += shift;
                                    const std::vector<libMesh::Point>& candidate_coords =
                                        intersection_u_points[axis][i_s_prime];
                                    const std::vector<libMesh::Point>& candidate_ref_coords =
                                        intersection_u_ref_coords[axis][i_s_prime];
                                    const std::vector<VectorValue<double> >& candidate_normals =
                                        intersection_u_normals[axis][i_s_prime];

                                    found_same_intersection_point =
                                        checkDoubleCountingIntersection(axis,
                                                                        dx,
                                                                        n,
                                                                        xu,
                                                                        xui,
                                                                        i_s_um,
                                                                        i_s_prime,
                                                                        candidate_coords,
                                                                        candidate_ref_coords,
                                                                        candidate_normals);
                                    if (found_same_intersection_point) break;
                                }

                                if (!found_same_intersection_point)
                                {
                                    // imposed jump conditions.
                                    TBOX_ASSERT(i_s_um.getAxis() == i_s_up.getAxis());
                                    TBOX_ASSERT(i_s_up(axis) - i_s_um(axis) == 1);
                                    const double x_cell_bdry_um =
                                        x_lower[axis] +
                                        static_cast<double>(i_s_um(axis) - patch_lower[axis]) * dx[axis];

                                    const double x_cell_bdry_up =
                                        x_lower[axis] +
                                        static_cast<double>(i_s_up(axis) - patch_lower[axis]) * dx[axis];
                                    const double sdh_um = ((xu(axis) - x_cell_bdry_um)); // Signed Distance h

                                    const double sdh_up = ((xu(axis) - x_cell_bdry_up)); // Signed Distance h
                                    TBOX_ASSERT((sdh_um) < dx[axis] && sdh_um > 0);
                                    TBOX_ASSERT(fabs(sdh_up) < dx[axis] && sdh_up < 0);
                                    // Make sure these indices are within the ghost box for the patch data.
                                    if (side_ghost_boxes[axis].contains(i_s_up) &&
                                        side_ghost_boxes[axis].contains(i_s_um))
                                    {
                                        double C_u_um = 0;
                                        double C_u_up = 0;

                                        interpolate(&jn(0), 0, DU_jump_node[axis], phi_P_jump);
                                        C_u_up = sdh_up * jn(axis);
                                        C_u_um = sdh_um * jn(axis);

                                        const double sgn = n(axis) > 0.0 ? 1.0 : n(axis) < 0.0 ? -1.0 : 0.0;
                                        // Note that the corrections are applied to opposite sides
                                        // No extra stress jumps because they were handled with the pressure term.
                                        (*f_data)(i_s_up) -= sgn * (C_u_um / (dx[axis] * dx[axis]));
                                        (*f_data)(i_s_um) += sgn * (C_u_up / (dx[axis] * dx[axis]));
                                    }

                                    // Keep track of the positions where we have
                                    // imposed jump conditions.
                                    intersection_u_points[axis][i_s_um].push_back(xu);
                                    intersection_u_ref_coords[axis][i_s_um].push_back(xui);
                                    intersection_u_normals[axis][i_s_um].push_back(n);
                                }
                            }
                        }

                        // Loop over the other axes
                        for (unsigned int j = 0; j < NDIM - 1; ++j)
                        {
                            for (unsigned int k = 0; k < intersectionsSide[j].size(); ++k)
                            {
                                // Note that rs[axis] has the correct component. The other component is found from the
                                // location of the intersection.
                                libMesh::Point xu = rs[j] + intersectionsSide[j][k].first * q;
                                const libMesh::Point& xui = intersectionsSide[j][k].second;
                                // Two different stencils are pierced by this intersection.
                                SideIndex<NDIM> i_s_up;
                                SideIndex<NDIM> i_s_um;

                                // Check that the intersection location is actually on this patch
                                if (xu(axis) - x_lower[axis] > 0.0)
                                {
                                    // Check to see which half of the cell the intersection is located.
                                    // Note that this is important for both the sign of the correction and to find the
                                    // correct stencils to be corrected. The second block in this if statement shifts
                                    // the cell back half a grid cell.
                                    if (fmod(xu(axis) - x_lower[axis], dx[axis]) >= 0.5 * dx[axis])
                                    {
                                        SideIndex<NDIM> i_side_um(i_c, SideDim[axis][j], 0);
                                        hier::Index<NDIM> i_c_neighbor = i_c;
                                        i_c_neighbor(axis) += 1;

                                        SideIndex<NDIM> i_side_up(i_c_neighbor, SideDim[axis][j], 0);

                                        i_side_up(axis) = boost::math::iround((xu(axis) - x_lower[axis]) / dx[axis]) +
                                                          patch_lower[axis];
                                        i_side_um(axis) =
                                            boost::math::iround((xu(axis) - x_lower[axis]) / dx[axis] - 0.5) +
                                            patch_lower[axis];
                                        i_s_up = i_side_up;
                                        i_s_um = i_side_um;
                                    }
                                    else if (fmod((xu(axis) - x_lower[axis]), dx[axis]) < 0.5 * dx[axis])
                                    {
                                        SideIndex<NDIM> i_side_up(i_c, SideDim[axis][j], 0);
                                        hier::Index<NDIM> i_c_neighbor = i_c;
                                        i_c_neighbor(axis) -= 1;
                                        SideIndex<NDIM> i_side_um(i_c_neighbor, SideDim[axis][j], 0);
                                        i_side_up(axis) =
                                            boost::math::iround((xu(axis) - x_lower[axis]) / dx[axis] - 0.5) +
                                            patch_lower[axis];
                                        i_side_um(axis) =
                                            boost::math::iround((xu(axis) - x_lower[axis]) / dx[axis] - 1.0) +
                                            patch_lower[axis];
                                        i_s_up = i_side_up;
                                        i_s_um = i_side_um;
                                    }
                                    else
                                    {
                                        continue;
                                    }
                                }
                                // What is the point of this if statement? The code blocks are exactly the same
                                else if (xu(axis) - x_lower[axis] < 0.0)
                                {
                                    if (fmod(fabs(xu(axis) - x_lower[axis]), dx[axis]) < 0.5 * dx[axis])
                                    {
                                        SideIndex<NDIM> i_side_um(i_c, SideDim[axis][j], 0);
                                        hier::Index<NDIM> i_c_neighbor = i_c;
                                        i_c_neighbor(axis) += 1;

                                        SideIndex<NDIM> i_side_up(i_c_neighbor, SideDim[axis][j], 0);

                                        i_side_up(axis) = boost::math::iround((xu(axis) - x_lower[axis]) / dx[axis]) +
                                                          patch_lower[axis];
                                        i_side_um(axis) =
                                            boost::math::iround((xu(axis) - x_lower[axis]) / dx[axis] - 0.5) +
                                            patch_lower[axis];
                                        i_s_up = i_side_up;
                                        i_s_um = i_side_um;
                                    }
                                    else
                                    {
                                        SideIndex<NDIM> i_side_up(i_c, SideDim[axis][j], 0);
                                        hier::Index<NDIM> i_c_neighbor = i_c;
                                        i_c_neighbor(axis) -= 1;
                                        SideIndex<NDIM> i_side_um(i_c_neighbor, SideDim[axis][j], 0);
                                        i_side_up(axis) =
                                            boost::math::iround((xu(axis) - x_lower[axis]) / dx[axis] - 0.5) +
                                            patch_lower[axis];
                                        i_side_um(axis) =
                                            boost::math::iround((xu(axis) - x_lower[axis]) / dx[axis] - 1.0) +
                                            patch_lower[axis];
                                        i_s_up = i_side_up;
                                        i_s_um = i_side_um;
                                    }
                                }
                                else
                                {
                                    TBOX_ERROR(d_object_name << ":  Restart file version different than class version."
                                                             << std::endl);
                                }

                                if (extended_side_box.contains(i_s_up) && extended_side_box.contains(i_s_um))
                                {
                                    TBOX_ASSERT(i_s_up(axis) - i_s_um(axis) == 1);
                                    std::vector<libMesh::Point> ref_coords(1, xui);
                                    fe_X->reinit(elem, &ref_coords);
                                    fe_P_jump->reinit(elem, &ref_coords);
                                    for (unsigned int l = 0; l < NDIM - 1; ++l)
                                    {
                                        interpolate(dx_dxi[l], 0, x_node, *dphi_dxi[l]);
                                    }
                                    if (NDIM == 2)
                                    {
                                        dx_dxi[1] = VectorValue<double>(0.0, 0.0, 1.0);
                                    }
                                    n = (dx_dxi[0].cross(dx_dxi[1])).unit();

                                    bool found_same_intersection_point = false;

                                    for (int shift = -1; shift <= 1; ++shift)
                                    {
                                        SideIndex<NDIM> i_s_prime = i_s_um;
                                        i_s_prime(SideDim[axis][j]) += shift;
                                        const std::vector<libMesh::Point>& candidate_coords =
                                            intersectionSide_u_points[j][axis][i_s_prime];
                                        const std::vector<libMesh::Point>& candidate_ref_coords =
                                            intersectionSide_u_ref_coords[j][axis][i_s_prime];
                                        const std::vector<VectorValue<double> >& candidate_normals =
                                            intersectionSide_u_normals[j][axis][i_s_prime];

                                        found_same_intersection_point =
                                            checkDoubleCountingIntersection(axis,
                                                                            dx,
                                                                            n,
                                                                            xu,
                                                                            xui,
                                                                            i_s_um,
                                                                            i_s_prime,
                                                                            candidate_coords,
                                                                            candidate_ref_coords,
                                                                            candidate_normals);
                                        if (found_same_intersection_point) break;
                                    }

                                    if (!found_same_intersection_point)
                                    {
                                        // Evaluate the jump conditions and apply them
                                        // to the Eulerian grid.

                                        const double x_mid_side_up =
                                            x_lower[axis] +
                                            static_cast<double>(i_s_up(axis) - patch_lower[axis] + 0.5) * dx[axis];

                                        const double x_mid_side_um =
                                            x_lower[axis] +
                                            static_cast<double>(i_s_um(axis) - patch_lower[axis] + 0.5) * dx[axis];

                                        TBOX_ASSERT(xu(axis) <= x_mid_side_up);
                                        TBOX_ASSERT(xu(axis) > x_mid_side_um);

                                        const double sdh_up = xu(axis) - x_mid_side_up; // Signed Distance h
                                        const double sdh_um = xu(axis) - x_mid_side_um;
                                        // Make sure these indices are within the ghost box for the patch data.
                                        if (side_ghost_boxes[SideDim[axis][j]].contains(i_s_up) &&
                                            side_ghost_boxes[SideDim[axis][j]].contains(i_s_um))
                                        {
                                            double C_u_um = 0;
                                            double C_u_up = 0;

                                            interpolate(&jn(0), 0, DU_jump_node[SideDim[axis][j]], phi_P_jump);
                                            C_u_um = sdh_um * jn(axis);
                                            C_u_up = sdh_up * jn(axis);

                                            int sig_comp = 0;
#if (NDIM == 2)
                                            // There's only one off diagonal term.
                                            sig_comp = 2;
#endif
#if (NDIM == 3)
                                            if (axis == 0)
                                            {
                                                if (SideDim[axis][j] == 1)
                                                    sig_comp = 3;
                                                else if (SideDim[axis][j] == 2)
                                                    sig_comp = 4;
                                                else
                                                    TBOX_ERROR("SHOULD NOT GET HERE!");
                                            }
                                            else if (axis == 1)
                                            {
                                                if (SideDim[axis][j] == 0)
                                                    sig_comp = 5;
                                                else if (SideDim[axis][j] == 2)
                                                    sig_comp = 3;
                                                else
                                                    TBOX_ERROR("SHOULD NOT GET HERE!");
                                            }
                                            else if (axis == 2)
                                            {
                                                if (SideDim[axis][j] == 0)
                                                    sig_comp = 4;
                                                else if (SideDim[axis][j] == 1)
                                                    sig_comp = 3;
                                                else
                                                    TBOX_ERROR("SHOULD NOT GET HERE!");
                                            }
                                            else
                                            {
                                                TBOX_ERROR("SHOULD NOT GET HERE!");
                                            }
#endif
                                            double sig_jump = 0.0;
                                            interpolate(sig_jump, 0, sig_jump_node[sig_comp], phi_P_jump);

                                            const double sgn = n(axis) > 0.0 ? 1.0 : n(axis) < 0.0 ? -1.0 : 0.0;

                                            (*f_data)(i_s_um) +=
                                                sgn * (C_u_up / (dx[axis] * dx[axis]) + sig_jump / dx[axis]);
                                            (*f_data)(i_s_up) -=
                                                sgn * (C_u_um / (dx[axis] * dx[axis]) + sig_jump / dx[axis]);
                                        }
                                        intersectionSide_u_points[j][axis][i_s_um].push_back(xu);
                                        intersectionSide_u_ref_coords[j][axis][i_s_um].push_back(xui);
                                        intersectionSide_u_normals[j][axis][i_s_um].push_back(n);
                                    }
                                }
                            }
                        }
                    }
                }
            }

            // Restore the element coordinates.
            for (unsigned int k = 0; k < n_nodes; ++k)
            {
                elem->point(k) = X_node_cache[k];
            }
        }
    }

    return;
} // imposeJumpConditions

/////////////////////////////// PRIVATE //////////////////////////////////////

void
CFIIMethod::commonConstructor(const std::string& object_name,
                              Pointer<Database> input_db,
                              const std::vector<libMesh::MeshBase*>& meshes,
                              int max_level_number,
                              bool register_for_restart,
                              const std::string& restart_read_dirname,
                              unsigned int restart_restore_number)
{
    // Initialize object with data read from the input and restart databases.
    bool from_restart = RestartManager::getManager()->isFromRestart();
    if (from_restart) getFromRestart();
    if (input_db) getFromInput(input_db, from_restart);

    // Setup timers.
    auto set_timer = [&](const char* name) { return TimerManager::getManager()->getTimer(name); };
    IBAMR_DO_ONCE(t_compute_lagrangian_force = set_timer("IBAMR::CFIIMMethod::computeLagrangianForce()");
                  t_spread_force = set_timer("IBAMR::CFIIMMethod::spreadForce()"););

    return;
} // commonConstructor

void
CFIIMethod::getFromInput(Pointer<Database> db, bool /*is_from_restart*/)
{
    d_use_extra_stress_jump_conditions = db->getBool("use_extra_stress_jump_conditions");
    return;
} // getFromInput

void
CFIIMethod::getFromRestart()
{
    Pointer<Database> restart_db = RestartManager::getManager()->getRootDatabase();
    Pointer<Database> db;
    if (restart_db->isDatabase(d_object_name))
    {
        db = restart_db->getDatabase(d_object_name);
    }
    else
    {
        TBOX_ERROR(d_object_name << ":  Restart database corresponding to " << d_object_name
                                 << " not found in restart file." << std::endl);
    }
    int ver = db->getInteger("CF_IIM_VERSION");
    if (ver != CF_IIM_VERSION)
    {
        TBOX_ERROR(d_object_name << ":  Restart file version different than class version." << std::endl);
    }
    return;
} // getFromRestart

template <typename Vector>
void
CFIIMethod::findInterpPts(const SAMRAI::pdat::CellIndex<NDIM>& eval_idx,
                          const int ls_idx,
                          SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                          const int stencil_size,
                          std::vector<Vector>& xpts,
                          std::vector<SAMRAI::pdat::CellIndex<NDIM> >& stencil)
{
    Pointer<CellData<NDIM, double> > ls_data = patch->getPatchData(ls_idx);
    const double ls_eval = (*ls_data)(eval_idx);
    Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
    const double* const dx = pgeom->getDx();
    const double* const xlow = pgeom->getXLower();
    const hier::Index<NDIM>& idx_low = patch->getBox().lower();
    std::queue<CellIndex<NDIM> > test_stencils;
    test_stencils.push(eval_idx);
    stencil.clear();
    stencil.reserve(stencil_size);
    xpts.clear();
    xpts.reserve(stencil_size);
    size_t i = 0;
    // TODO: This isn't really a flooding algorithm... All the points we append to test_stencils are qualified, so
    // there's no need to do this.
    while (i < stencil_size)
    {
        // Test the top of test_stencils
        CellIndex<NDIM> test_idx = test_stencils.front();
        // If this index doesn't exist in the stencil and is the same sign as the base index, add it
        if (std::find(stencil.begin(), stencil.end(), test_idx) == stencil.end() &&
            (((*ls_data)(test_idx)*ls_eval) > 0.0))
        {
            Vector X;
            for (int d = 0; d < NDIM; ++d)
                X(d) = xlow[d] + dx[d] * (static_cast<double>(test_idx(d) - idx_low(d)) + 0.5);
            stencil.push_back(test_idx);
            xpts.push_back(X);
            i++;
        }
        // Now append more indices to test_stencils
        for (int d = 0; d < NDIM; ++d)
        {
            IntVector<NDIM> dir(0);
            dir(d) = 1;
            test_stencils.push(test_idx + dir);
            test_stencils.push(test_idx - dir);
        }
        test_stencils.pop();
    }
}

template <typename Vector>
std::vector<double>
CFIIMethod::findLinearInterpWeights(const Vector& eval_pt, const std::vector<Vector>& xpts, const double scale)
{
    // Form a polyharmonic spline interpolant with linear polynomials
    const int m = xpts.size();
    IBTK::MatrixXd A(MatrixXd::Zero(m, m));
    IBTK::MatrixXd B(MatrixXd::Zero(m, NDIM + 1));
    IBTK::VectorXd U(VectorXd::Zero(m + NDIM + 1));
    for (size_t i = 0; i < m; ++i)
    {
        for (size_t j = 0; j < m; ++j)
        {
            // Use a third order rbf
            A(i, j) = std::pow((xpts[i] - xpts[j]).norm(), 3.0);
        }
        // Now evaluate polynomial. We use shifted and scaled monomials.
        B(i, 0) = 1.0;
        for (int d = 0; d < NDIM; ++d) B(i, d + 1) = (xpts[i](d) - eval_pt(d)) / scale;
        // RHS
        U(i) = std::pow((xpts[i] - eval_pt).norm(), 3.0);
    }
    // Finish RHS, should be polynomials evaluated at xpt. Note, since they are shifted and scaled, only the constant
    // term is non-zero.
    U(m) = 1.0;
    // Now build the final matrix
    MatrixXd final_mat(MatrixXd::Zero(m + NDIM + 1, m + NDIM + 1));
    final_mat.block(0, 0, m, m) = A;
    final_mat.block(0, m, m, NDIM + 1) = B;
    final_mat.block(m, 0, NDIM + 1, m) = B.transpose();
    VectorXd x = final_mat.colPivHouseholderQr().solve(U);

    // Now evaluate the stencil at x_in
    std::vector<double> wgts(m);
    for (int i = 0; i < m; ++i) wgts[i] = x[i];
    return wgts;
}

double
CFIIMethod::evaluateInterpolant(const CellData<NDIM, double>& data,
                                const std::vector<CellIndex<NDIM> >& idxs,
                                const std::vector<double>& wgts,
                                const int depth)
{
    double ret = 0.0;
#ifndef NDEBUG
    TBOX_ASSERT(idxs.size() == wgts.size());
#endif
    for (size_t i = 0; i < idxs.size(); ++i)
    {
#ifndef NDEBUG
        data.getGhostBox().contains(idxs[i]);
#endif
        ret += data(idxs[i], depth) * wgts[i];
    }
    return ret;
}

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
