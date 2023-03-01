#include "ibtk/IndexUtilities.h"
#include "ibtk/InterfaceLocator.h"
#include "ibtk/ls_utilities.h"

#include <CartesianGridGeometry.h>
#include <CellData.h>
#include <CellVariable.h>
#include <CoarsenAlgorithm.h>
#include <CoarsenOperator.h>
#include <CoarsenSchedule.h>
#include <SideData.h>
#include <SideVariable.h>

#include "ibtk/app_namespaces.h"

namespace IBTK
{

InterfaceLocator::InterfaceLocator(std::string object_name,
                                   Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
                                   std::vector<FEDataManager*> fe_data_managers)
    : d_object_name(std::move(object_name)),
      d_hierarchy(patch_hierarchy),
      d_fe_data_managers(std::move(fe_data_managers))
{
    d_flip_normal.resize(d_fe_data_managers.size(), 0);
}

void
InterfaceLocator::resetGlobalLSSign(const int ls_idx,
                                    Pointer<hier::Variable<NDIM> > ls_var,
                                    const double base_val,
                                    const double time)
{
    // First reset all the points
    int coarsest_ln = 0;
    int finest_ln = d_hierarchy->getFinestLevelNumber();
    TBOX_ASSERT(ls_idx != IBTK::invalid_index);
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            Pointer<CellData<NDIM, double> > ls_cc_data = patch->getPatchData(ls_idx);
            Pointer<SideData<NDIM, double> > ls_sc_data = patch->getPatchData(ls_idx);
            if (ls_cc_data) ls_cc_data->fillAll(base_val + ln + 1);
            if (ls_sc_data) ls_sc_data->fillAll(base_val + ln + 1);
        }
    }

    // Update the local sign
    updateLocalLSSign(ls_idx, base_val + finest_ln, time);

    // Now do a flood filling algorithm
    Pointer<CellVariable<NDIM, double> > ls_cc_var = ls_var;
    Pointer<SideVariable<NDIM, double> > ls_sc_var = ls_var;
    // We assume the finest level has been filled so that we have an interface.
    if (ls_cc_var)
        flood_fill_on_level_cell(
            ls_idx, d_hierarchy->getPatchLevel(finest_ln), static_cast<double>(base_val + finest_ln + 1), time);
    if (ls_sc_var)
        flood_fill_on_level_side(
            ls_idx, d_hierarchy->getPatchLevel(finest_ln), static_cast<double>(base_val + finest_ln + 1), time);
    // Now that the finest level has the correct sign, we fill coarser levels.
    Pointer<CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();
    Pointer<CoarsenOperator<NDIM> > coarsen_op = grid_geom->lookupCoarsenOperator(ls_var, "CONSERVATIVE_COARSEN");
    Pointer<CoarsenAlgorithm<NDIM> > coarsen_alg = new CoarsenAlgorithm<NDIM>();
    coarsen_alg->registerCoarsen(ls_idx, ls_idx, coarsen_op);
    std::vector<Pointer<CoarsenSchedule<NDIM> > > coarsen_scheds(finest_ln + 1);
    for (int ln = finest_ln; ln > 0; --ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        Pointer<PatchLevel<NDIM> > coarser_level = d_hierarchy->getPatchLevel(ln - 1);
        coarsen_scheds[ln] = coarsen_alg->createSchedule(coarser_level, level);
        coarsen_scheds[ln]->coarsenData();
        if (ls_cc_var) flood_fill_on_level_cell(ls_idx, coarser_level, base_val + ln + 1, time);
        if (ls_sc_var) flood_fill_on_level_side(ls_idx, coarser_level, base_val + ln + 1, time);
    }
}

void
InterfaceLocator::updateLocalLSSign(const int ls_idx, const double base_val, const double time)
{
    unsigned int num_parts = d_fe_data_managers.size();
    int coarsest_ln = 0;
    int finest_ln = d_hierarchy->getFinestLevelNumber();
#ifndef NDEBUG
    // First make sure ls_idx has data allocated
    TBOX_ASSERT(ls_idx != IBTK::invalid_index);
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        TBOX_ASSERT(level->checkAllocated(ls_idx));
    }
#endif

    // Loop through elements, find bounding box. For each bounding box, determine which side of the boundary we are on.
    for (unsigned int part = 0; part < num_parts; ++part)
    {
        EquationSystems* eq_sys = d_fe_data_managers[part]->getEquationSystems();
        const MeshBase& mesh = eq_sys->get_mesh();

        System& X_sys = eq_sys->get_system(d_fe_data_managers[part]->COORDINATES_SYSTEM_NAME);
        NumericVector<double>* X_vec = X_sys.current_local_solution.get();
        NumericVector<double>* X_ghost_vec = d_fe_data_managers[part]->buildGhostedCoordsVector(false);
        copy_and_synch(*X_vec, *X_ghost_vec);
        const DofMap& X_dof_map = X_sys.get_dof_map();
        FEDataManager::SystemDofMapCache& X_dof_map_cache =
            *d_fe_data_managers[part]->getDofMapCache(d_fe_data_managers[part]->COORDINATES_SYSTEM_NAME);
        FEType X_fe_type = X_dof_map.variable_type(0);

        // We need any point on the structure, we'll just use a first order nodal rule (which should be a single point,
        // the midpoint).
        // TODO: This won't work with quadratic or higher structures. Unsure what to do there...
        std::unique_ptr<QBase> qrule = QBase::build("QNodal", mesh.mesh_dimension(), FIRST);
        std::unique_ptr<FEBase> fe = FEBase::build(mesh.mesh_dimension(), X_fe_type);
        const std::vector<std::vector<double> >& phi = fe->get_phi();
        std::array<const std::vector<std::vector<double> >*, NDIM - 1> dphi_dxi;
        dphi_dxi[0] = &fe->get_dphidxi();
        if (NDIM > 2) dphi_dxi[1] = &fe->get_dphideta();
        fe->attach_quadrature_rule(qrule.get());

        boost::multi_array<double, 2> x_node, X_node, n_qp_node;
        std::array<VectorValue<double>, 2> dx_dxi;
        VectorValue<double> n, N, x, X;
        IBTK::Point x_min, x_max;

        // TODO: We need to be careful if the structure is not on the finest level.
        TBOX_ASSERT(d_fe_data_managers[part]->getFinestPatchLevelNumber() == finest_ln);
        Pointer<PatchLevel<NDIM> > level =
            d_hierarchy->getPatchLevel(d_fe_data_managers[part]->getFinestPatchLevelNumber());
        const Pointer<CartesianGridGeometry<NDIM> > grid_geom = level->getGridGeometry();
        const IntVector<NDIM>& ratio = level->getRatio();
        int local_patch_num = 0;
        for (PatchLevel<NDIM>::Iterator p(level); p; p++, ++local_patch_num)
        {
            const std::vector<Elem*>& patch_elems =
                d_fe_data_managers[part]->getActivePatchElementMap()[local_patch_num];
            const size_t num_active_patch_elems = patch_elems.size();
            if (num_active_patch_elems == 0) continue;

            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            Pointer<CellData<NDIM, double> > ls_cc_data = patch->getPatchData(ls_idx);
            Pointer<SideData<NDIM, double> > ls_sc_data = patch->getPatchData(ls_idx);
            Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
            const double* const dx = pgeom->getDx();
            const double* const xlow = pgeom->getXLower();
            const double dx_min = *std::min_element(dx, dx + NDIM);

            const Box<NDIM>& patch_box = patch->getBox();
            const hier::Index<NDIM>& idx_low = patch_box.lower();

            for (const auto& elem : patch_elems)
            {
                const auto& X_dof_indices = X_dof_map_cache.dof_indices(elem);
                get_values_for_interpolation(x_node, *X_ghost_vec, X_dof_indices);
                const unsigned int n_nodes = elem->n_nodes();
                x_min.setConstant(std::numeric_limits<double>::max());
                x_max.setConstant(std::numeric_limits<double>::lowest());
                // Determine element node locations
                for (unsigned int k = 0; k < n_nodes; ++k)
                {
                    for (unsigned int d = 0; d < NDIM; ++d) x(d) = x_node[k][d];
                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        x_min[d] = std::min(x_min[d], x(d));
                        x_max[d] = std::max(x_max[d], x(d));
                    }
                }
                // We use a first order quadrature rule to find the midpoint of the element.
                fe->reinit(elem);
                interpolate(x, 0, x_node, phi);
                for (unsigned int l = 0; l < NDIM - 1; ++l) interpolate(dx_dxi[l], 0, x_node, *dphi_dxi[l]);
                if (NDIM == 2) dx_dxi[1] = VectorValue<double>(0.0, 0.0, 1.0);
                // We need the normal to detect sides.
                n = (dx_dxi[0].cross(dx_dxi[1])).unit();
                // NOTE: While we don't really care which side is which, we need to make sure sides are consistent with
                // other parts. Flipping the normal is necessary in cases where there are multiple parts.
                if (d_flip_normal[part] == 1) n *= -1.0;

                // Compute bounding box that is interior to the patch box.
                Box<NDIM> box(IndexUtilities::getCellIndex(x_min, grid_geom, ratio),
                              IndexUtilities::getCellIndex(x_max, grid_geom, ratio));
                box.grow(IntVector<NDIM>(1));
                box = box * patch_box;

                // Now loop over cell indices and compute signs.
                // This is a very easy check if the element is linear. We compute the normal of the element, create the
                // vector from ANY point on the element to the index. If the dot product is positive, these vectors face
                // the same direction. Otherwise they face a different direction. Then set the value of the level set to
                // the sign of the dot product.
                // TODO: In finite precision, we should really check something close to the normal projection.
                // TODO: What happens if the element is quadratic or higher?
                // TODO: I think we need to compute distances for when points are on different sides of the element.
                if (ls_cc_data)
                {
                    for (CellIterator<NDIM> ci(box); ci; ci++)
                    {
                        const CellIndex<NDIM>& idx = ci();
                        VectorValue<double> xc(0.0);
                        for (int d = 0; d < NDIM; ++d)
                            xc(d) = xlow[d] + dx[d] * (static_cast<double>(idx(d) - idx_low(d)) + 0.5);
                        VectorValue<double> x_to_xc = xc - x;
                        double signed_val = ((x_to_xc * n) > 0.0 ? 1.0 : -1.0) * base_val;
                        (*ls_cc_data)(idx) = signed_val;
                    }
                }
                else if (ls_sc_data)
                {
                    for (int axis = 0; axis < NDIM; ++axis)
                    {
                        for (SideIterator<NDIM> si(box, axis); si; si++)
                        {
                            const SideIndex<NDIM>& idx = si();
                            VectorValue<double> xc(0.0);
                            for (int d = 0; d < NDIM; ++d)
                                xc(d) = xlow[d] +
                                        dx[d] * (static_cast<double>(idx(d) - idx_low(d)) + (d == axis ? 0.0 : 0.5));
                            VectorValue<double> x_to_xc = xc - x;
                            (*ls_sc_data)(idx) = ((x_to_xc * n) > 0.0 ? 1.0 : -1.0) * base_val;
                        }
                    }
                }
            }

            // Fill in boundaries
            if (d_bdry_fcn)
            {
                std::vector<Box<NDIM> > fill_boxes;
                IntVector<NDIM> gcw_to_fill;
                if (ls_cc_data) gcw_to_fill = ls_cc_data->getGhostCellWidth();
                if (ls_sc_data) gcw_to_fill = ls_sc_data->getGhostCellWidth();
                for (int d = 1; d <= NDIM; ++d)
                {
                    const tbox::Array<BoundaryBox<NDIM> >& bdry_boxes = pgeom->getCodimensionBoundaries(d);
                    for (int i = 0; i < bdry_boxes.size(); ++i)
                    {
                        const BoundaryBox<NDIM>& bdry_box = bdry_boxes[i];
                        const int location_index = bdry_box.getLocationIndex();
                        const int axis = location_index % 2;
                        const int upper_lower = location_index / 2;
                        if (pgeom->getTouchesRegularBoundary(axis, upper_lower))
                            fill_boxes.push_back(pgeom->getBoundaryFillBox(bdry_box, patch->getBox(), gcw_to_fill));
                    }
                }
                for (const auto& box : fill_boxes)
                {
                    if (ls_cc_data)
                    {
                        for (CellIterator<NDIM> ci(box); ci; ci++)
                        {
                            const CellIndex<NDIM>& idx = ci();
                            VectorNd x;
                            for (int d = 0; d < NDIM; ++d)
                                x[d] = xlow[d] + dx[d] * (static_cast<double>(idx(d) - idx_low(d)) + 0.5);
                            (*ls_cc_data)(idx) = d_bdry_fcn(x, time, d_bdry_ctx);
                        }
                    }
                    else if (ls_sc_data)
                    {
                        for (int axis = 0; axis < NDIM; ++axis)
                        {
                            for (SideIterator<NDIM> si(box, axis); si; si++)
                            {
                                const SideIndex<NDIM>& idx = si();
                                VectorNd x;
                                for (int d = 0; d < NDIM; ++d)
                                    x[d] = xlow[d] +
                                           dx[d] * (static_cast<double>(idx(d) - idx_low(d)) + (d == axis ? 0.0 : 0.5));
                                (*ls_sc_data)(idx) = d_bdry_fcn(x, time, d_bdry_ctx);
                            }
                        }
                    }
                }
            }
        }
    }
    return;
}
} // namespace IBTK
