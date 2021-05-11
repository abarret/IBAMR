#ifndef included_CutCellVolumeMeshMapping
#define included_CutCellVolumeMeshMapping
#include "ibamr/cut_cells/CutCellMeshMapping.h"
#include "ibamr/cut_cells/FEMeshPartitioner.h"
#include "ibamr/cut_cells/ls_functions.h"
#include "ibamr/cut_cells/ls_utilities.h"

#include "ibtk/FEDataManager.h"

#include "libmesh/boundary_mesh.h"
#include "libmesh/mesh.h"

namespace LS
{
/*!
 * CutCellVolumeMeshMapping maintains a description of the Lagrangian mesh from the point of view of the background
 * mesh. We maintain a mapping from each cut cell index to a vector of element and element parent pairs.
 */
class CutCellVolumeMeshMapping : public CutCellMeshMapping
{
public:
    /*!
     * \brief Constructor.
     */
    CutCellVolumeMeshMapping(std::string object_name,
                             SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                             libMesh::MeshBase* vol_mesh,
                             IBTK::FEDataManager* vol_fe_data_manager,
                             const std::set<libMesh::boundary_id_type>& bdry_id);

    /*!
     * \brief Constructor.
     */
    CutCellVolumeMeshMapping(std::string object_name,
                             SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                             const std::vector<libMesh::MeshBase*>& vol_meshes,
                             const std::vector<IBTK::FEDataManager*>& vol_fe_data_managers,
                             const std::vector<std::set<libMesh::boundary_id_type> >& bdry_ids);

    /*!
     * \brief Default deconstructor.
     */
    virtual ~CutCellVolumeMeshMapping();

    /*!
     * \brief Deleted default constructor.
     */
    CutCellVolumeMeshMapping() = delete;

    /*!
     * \brief Deleted copy constructor.
     */
    CutCellVolumeMeshMapping(const CutCellVolumeMeshMapping& from) = delete;

    /*!
     * \brief Deleted assignment operator.
     */
    CutCellVolumeMeshMapping& operator=(const CutCellVolumeMeshMapping& that) = delete;

    void generateCutCellMappings() override;

    inline void registerMappingFcn(MappingFcn fcn, unsigned int part = 0)
    {
        d_mapping_fcns[part] = fcn;
    }

    inline void registerBdryIdToSkip(unsigned int bdry_id, unsigned int part = 0)
    {
        d_bdry_id_to_skip_vec[part].insert(bdry_id);
    }

    inline std::shared_ptr<FEMeshPartitioner>& getMeshPartitioner(unsigned int part = 0)
    {
        return d_bdry_mesh_partitioners[part];
    }

private:
    void commonConstructor(const std::vector<std::set<libMesh::boundary_id_type> >& bdry_ids,
                           SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db);

    bool findIntersection(libMesh::Point& p, libMesh::Elem* elem, libMesh::Point r, libMesh::VectorValue<double> q);

    void matchBoundaryToVolume(unsigned int part = 0);

    std::vector<libMesh::MeshBase*> d_vol_meshes;
    std::vector<IBTK::FEDataManager*> d_vol_fe_data_managers;

    std::vector<std::shared_ptr<IBTK::FEData> > d_fe_data;
    std::vector<std::shared_ptr<FEMeshPartitioner> > d_bdry_mesh_partitioners;
    std::vector<std::unique_ptr<libMesh::BoundaryMesh> > d_bdry_meshes;
    std::vector<std::unique_ptr<libMesh::EquationSystems> > d_bdry_eq_sys_vec;
    std::string d_coords_sys_name = "COORDINATES_SYSTEM";
    std::string d_disp_sys_name = "DISPLACEMENT_SYSTEM";

    std::vector<MappingFcn> d_mapping_fcns;

    std::vector<std::set<unsigned int> > d_bdry_id_to_skip_vec;

    std::vector<std::vector<libMesh::Elem*> > d_active_patch_elem_map;
};

} // namespace LS
#endif
