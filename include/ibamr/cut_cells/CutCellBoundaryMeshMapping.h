#ifndef included_CutCellBoundaryMeshMapping
#define included_CutCellBoundaryMeshMapping
#include "ibamr/cut_cells/CutCellMeshMapping.h"
#include "ibamr/cut_cells/ls_functions.h"
#include "ibamr/cut_cells/ls_utilities.h"

#include "ibtk/FEDataManager.h"

#include "libmesh/mesh.h"

namespace LS
{
/*!
 * CutCellMeshMapping maintains a description of the Lagrangian mesh from the point of view of the background mesh. We
 * maintain a mapping from each cut cell index to a vector of element and element parent pairs.
 */
class CutCellBoundaryMeshMapping : public CutCellMeshMapping
{
public:
    /*!
     * \brief Constructor.
     */
    CutCellBoundaryMeshMapping(std::string object_name,
                               SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                               libMesh::MeshBase* mesh,
                               IBTK::FEDataManager* fe_data_manager);

    /*!
     * \brief Constructor.
     */
    CutCellBoundaryMeshMapping(std::string object_name,
                               SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                               const std::vector<libMesh::MeshBase*>& meshes,
                               const std::vector<IBTK::FEDataManager*>& fe_data_managers);

    /*!
     * \brief Default deconstructor.
     */
    virtual ~CutCellBoundaryMeshMapping();

    /*!
     * \brief Deleted default constructor.
     */
    CutCellBoundaryMeshMapping() = delete;

    /*!
     * \brief Deleted copy constructor.
     */
    CutCellBoundaryMeshMapping(const CutCellBoundaryMeshMapping& from) = delete;

    /*!
     * \brief Deleted assignment operator.
     */
    CutCellBoundaryMeshMapping& operator=(const CutCellBoundaryMeshMapping& that) = delete;

    void generateCutCellMappings() override;

    inline void registerMappingFcn(MappingFcn fcn, unsigned int part = 0)
    {
        d_mapping_fcns[part] = fcn;
    }

    inline void registerBdryIdToSkip(unsigned int bdry_id, unsigned int part = 0)
    {
        d_bdry_id_to_skip_vec[part].insert(bdry_id);
    }

private:
    bool findIntersection(libMesh::Point& p, libMesh::Elem* elem, libMesh::Point r, libMesh::VectorValue<double> q);

    std::vector<libMesh::MeshBase*> d_meshes;
    std::vector<IBTK::FEDataManager*> d_fe_data_managers;

    std::vector<MappingFcn> d_mapping_fcns;
    std::vector<std::set<unsigned int> > d_bdry_id_to_skip_vec;
};

} // namespace LS
#endif
