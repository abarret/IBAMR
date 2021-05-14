#include "ibamr/cut_cells/VolumeBoundaryMeshMapping.h"

#include "ibtk/IBTK_MPI.h"
#include "ibtk/IndexUtilities.h"

#include "libmesh/explicit_system.h"

#include "ibamr/app_namespaces.h"

namespace LS
{
VolumeBoundaryMeshMapping::VolumeBoundaryMeshMapping(std::string object_name,
                                                     Pointer<Database> input_db,
                                                     Pointer<PatchHierarchy<NDIM> > hierarchy,
                                                     MeshBase* mesh,
                                                     FEDataManager* fe_data_manager,
                                                     const std::vector<std::set<boundary_id_type> >& bdry_id)
    : d_object_name(std::move(object_name)),
      d_hierarchy(hierarchy),
      d_vol_meshes({ mesh }),
      d_vol_fe_data_managers({ fe_data_manager })
{
    std::vector<unsigned int> parts(bdry_id.size(), 0);

    commonConstructor(bdry_id, parts, input_db);
}

VolumeBoundaryMeshMapping::VolumeBoundaryMeshMapping(std::string object_name,
                                                     Pointer<Database> input_db,
                                                     Pointer<PatchHierarchy<NDIM> > hierarchy,
                                                     const std::vector<MeshBase*>& meshes,
                                                     const std::vector<FEDataManager*>& fe_data_managers,
                                                     const std::vector<std::set<boundary_id_type> >& bdry_ids,
                                                     const std::vector<unsigned int>& part)
    : d_object_name(std::move(object_name)),
      d_hierarchy(hierarchy),
      d_vol_meshes(meshes),
      d_vol_fe_data_managers(fe_data_managers)
{
    commonConstructor(bdry_ids, part, input_db);
}

void
VolumeBoundaryMeshMapping::commonConstructor(const std::vector<std::set<boundary_id_type> >& bdry_ids,
                                             const std::vector<unsigned int>& parts,
                                             Pointer<Database> input_db)
{
    for (const auto& part : parts)
    {
        std::unique_ptr<BoundaryMesh> bdry_mesh =
            libmesh_make_unique<BoundaryMesh>(d_vol_meshes[part]->comm(), d_vol_meshes[part]->spatial_dimension() - 1);
        d_vol_meshes[part]->boundary_info->sync(bdry_ids[part], *bdry_mesh);
        d_bdry_meshes.push_back(std::move(bdry_mesh));
        d_bdry_eq_sys_vec.push_back(std::move(libmesh_make_unique<EquationSystems>(*d_bdry_meshes[part])));
        d_fe_data.push_back(
            std::make_shared<FEData>(d_object_name + "::" + std::to_string(part), *d_bdry_eq_sys_vec[part], true));

        // TODO: Need to read this from restart files
        auto& X_sys = d_bdry_eq_sys_vec[part]->add_system<ExplicitSystem>(d_coords_sys_name);
        for (unsigned int d = 0; d < NDIM; ++d) X_sys.add_variable("X_" + std::to_string(d), FEType());
        auto& dX_sys = d_bdry_eq_sys_vec[part]->add_system<ExplicitSystem>(d_disp_sys_name);
        for (unsigned int d = 0; d < NDIM; ++d) dX_sys.add_variable("dX_" + std::to_string(d), FEType());
        X_sys.assemble_before_solve = false;
        X_sys.assemble();
        dX_sys.assemble_before_solve = false;
        dX_sys.assemble();
        d_bdry_mesh_partitioners.push_back(
            std::make_shared<FEMeshPartitioner>(d_object_name + "_" + std::to_string(part),
                                                input_db,
                                                input_db->getInteger("max_level"),
                                                IntVector<NDIM>(0),
                                                d_fe_data[part],
                                                d_coords_sys_name,
                                                false));
    }
    return;
}

void
VolumeBoundaryMeshMapping::matchBoundaryToVolume()
{
    for (unsigned int part = 0; part < d_bdry_meshes.size(); ++part) matchBoundaryToVolume(part);
    return;
}

void
VolumeBoundaryMeshMapping::matchBoundaryToVolume(unsigned int part)
{
    FEDataManager* fe_data_manager = d_vol_fe_data_managers[part];
    EquationSystems* eq_sys = fe_data_manager->getEquationSystems();

    System& X_system = eq_sys->get_system(fe_data_manager->COORDINATES_SYSTEM_NAME);
    const DofMap& X_dof_map = X_system.get_dof_map();
    NumericVector<double>* X_vec = X_system.solution.get();
    auto X_petsc_vec = dynamic_cast<PetscVector<double>*>(X_vec);
    TBOX_ASSERT(X_petsc_vec != nullptr);
    const double* const X_local_soln = X_petsc_vec->get_array_read();
    FEDataManager::SystemDofMapCache& X_dof_map_cache =
        *fe_data_manager->getDofMapCache(fe_data_manager->COORDINATES_SYSTEM_NAME);

    System& X_bdry_sys = d_bdry_eq_sys_vec[part]->get_system(d_coords_sys_name);
    const DofMap& X_bdry_dof_map = X_bdry_sys.get_dof_map();
    NumericVector<double>* X_bdry_vec = X_bdry_sys.solution.get();

    System& dX_bdry_sys = d_bdry_eq_sys_vec[part]->get_system(d_disp_sys_name);
    const DofMap& dX_bdry_dof_map = X_bdry_sys.get_dof_map();
    NumericVector<double>* dX_bdry_vec = dX_bdry_sys.solution.get();

    std::map<dof_id_type, dof_id_type> node_id_map;
    std::map<dof_id_type, unsigned char> side_id_map;
    d_vol_meshes[part]->boundary_info->get_side_and_node_maps(*d_bdry_meshes[part], node_id_map, side_id_map);
    auto node_it = d_bdry_meshes[part]->local_nodes_begin();
    auto node_end = d_bdry_meshes[part]->local_nodes_end();
    for (; node_it != node_end; ++node_it)
    {
        Node* node = *node_it;
        dof_id_type bdry_node_id = node->id();
        auto vol_iter = std::find_if(node_id_map.begin(), node_id_map.end(), [bdry_node_id](const auto& obj) {
            return obj.second == bdry_node_id;
        });
        dof_id_type vol_node_id = vol_iter->first;
        // Grab current position of volumetric mesh.
        std::vector<dof_id_type> X_dof_indices, X_bdry_dof_indices;
        for (int d = 0; d < NDIM; ++d)
        {
            X_dof_map.dof_indices(d_vol_meshes[part]->node_ptr(vol_node_id), X_dof_indices, d);
            X_bdry_dof_map.dof_indices(node, X_bdry_dof_indices, d);
            X_bdry_vec->set(X_bdry_dof_indices[0], (*X_vec)(X_dof_indices[0]));
            dX_bdry_vec->set(X_bdry_dof_indices[0], (*X_vec)(X_dof_indices[0]) - (*node)(d));
        }
    }
    X_petsc_vec->restore_array();
    X_bdry_vec->close();
    dX_bdry_vec->close();
    return;
}

void
VolumeBoundaryMeshMapping::initializeEquationSystems()
{
    for (unsigned int part = 0; part < d_bdry_meshes.size(); ++part)
    {
        d_bdry_eq_sys_vec[part]->init();

        matchBoundaryToVolume(part);
    }
    return;
}
} // namespace LS
