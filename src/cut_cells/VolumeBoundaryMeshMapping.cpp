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
                                             const std::vector<unsigned int>& vol_parts,
                                             Pointer<Database> input_db)
{
    unsigned int bdry_mesh_id = 0;
    for (const auto& vol_part : vol_parts)
    {
        d_vol_id_vec.push_back(vol_part);
        std::unique_ptr<BoundaryMesh> bdry_mesh = libmesh_make_unique<BoundaryMesh>(
            d_vol_meshes[vol_part]->comm(), d_vol_meshes[vol_part]->spatial_dimension() - 1);
        d_vol_meshes[vol_part]->boundary_info->sync(bdry_ids[bdry_mesh_id], *bdry_mesh);
        d_bdry_meshes.push_back(std::move(bdry_mesh));
        d_bdry_eq_sys_vec.push_back(std::move(libmesh_make_unique<EquationSystems>(*d_bdry_meshes[bdry_mesh_id])));
        d_fe_data.push_back(std::make_shared<FEData>(
            d_object_name + "::" + std::to_string(bdry_mesh_id), *d_bdry_eq_sys_vec[bdry_mesh_id], true));

        // TODO: Need to read this from restart files
        auto& X_sys = d_bdry_eq_sys_vec[bdry_mesh_id]->add_system<ExplicitSystem>(d_coords_sys_name);
        for (unsigned int d = 0; d < NDIM; ++d) X_sys.add_variable("X_" + std::to_string(d), FEType());
        auto& dX_sys = d_bdry_eq_sys_vec[bdry_mesh_id]->add_system<ExplicitSystem>(d_disp_sys_name);
        for (unsigned int d = 0; d < NDIM; ++d) dX_sys.add_variable("dX_" + std::to_string(d), FEType());
        X_sys.assemble_before_solve = false;
        X_sys.assemble();
        dX_sys.assemble_before_solve = false;
        dX_sys.assemble();
        d_bdry_mesh_partitioners.push_back(
            std::make_shared<FEMeshPartitioner>(d_object_name + "_" + std::to_string(bdry_mesh_id),
                                                input_db,
                                                input_db->getInteger("max_level"),
                                                IntVector<NDIM>(0),
                                                d_fe_data[bdry_mesh_id],
                                                d_coords_sys_name,
                                                false));
        ++bdry_mesh_id;
    }
    return;
}

void
VolumeBoundaryMeshMapping::matchBoundaryToVolume(std::string sys_name)
{
    for (unsigned int part = 0; part < d_bdry_meshes.size(); ++part) matchBoundaryToVolume(part, sys_name);
    return;
}

void
VolumeBoundaryMeshMapping::matchBoundaryToVolume(unsigned int part, std::string sys_name)
{
    FEDataManager* fe_data_manager = d_vol_fe_data_managers[d_vol_id_vec[part]];
    EquationSystems* eq_sys = fe_data_manager->getEquationSystems();

    System& X_system = eq_sys->get_system(fe_data_manager->COORDINATES_SYSTEM_NAME);
    const DofMap& X_dof_map = X_system.get_dof_map();
    NumericVector<double>* X_vec;
    if (sys_name == "")
        X_vec = X_system.solution.get();
    else
        X_vec = &X_system.get_vector(sys_name);

    System& X_bdry_sys = d_bdry_eq_sys_vec[part]->get_system(d_coords_sys_name);
    const DofMap& X_bdry_dof_map = X_bdry_sys.get_dof_map();
    NumericVector<double>* X_bdry_vec = X_bdry_sys.solution.get();

    System& dX_bdry_sys = d_bdry_eq_sys_vec[part]->get_system(d_disp_sys_name);
    NumericVector<double>* dX_bdry_vec = dX_bdry_sys.solution.get();

    std::map<dof_id_type, dof_id_type> node_id_map;
    std::map<dof_id_type, unsigned char> side_id_map;
    d_vol_meshes[d_vol_id_vec[part]]->boundary_info->get_side_and_node_maps(
        *d_bdry_meshes[part], node_id_map, side_id_map);
    auto node_it = d_bdry_meshes[part]->local_nodes_begin();
    auto node_end = d_bdry_meshes[part]->local_nodes_end();
    for (; node_it != node_end; ++node_it)
    {
        Node* node = *node_it;
        dof_id_type bdry_node_id = node->id();
        // This is potentially expensive. We should cache our own map between bdry nodes and volumetric nodes.
        auto vol_iter = std::find_if(node_id_map.begin(), node_id_map.end(), [bdry_node_id](const auto& obj) {
            return obj.second == bdry_node_id;
        });
        dof_id_type vol_node_id = vol_iter->first;
        // Grab current position of volumetric mesh.
        std::vector<dof_id_type> X_dof_indices, X_bdry_dof_indices;
        for (int d = 0; d < NDIM; ++d)
        {
            X_dof_map.dof_indices(d_vol_meshes[d_vol_id_vec[part]]->node_ptr(vol_node_id), X_dof_indices, d);
            X_bdry_dof_map.dof_indices(node, X_bdry_dof_indices, d);
            X_bdry_vec->set(X_bdry_dof_indices[0], (*X_vec)(X_dof_indices[0]));
            dX_bdry_vec->set(X_bdry_dof_indices[0], (*X_vec)(X_dof_indices[0]) - (*node)(d));
        }
    }
    X_bdry_vec->close();
    dX_bdry_vec->close();
    X_bdry_sys.update();
    dX_bdry_sys.update();
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
