#ifndef included_IBTK_InterfaceLocator
#define included_IBTK_InterfaceLocator

#include "ibtk/FEDataManager.h"

#include <tbox/Pointer.h>

#include <PatchHierarchy.h>
#include <Variable.h>

#include <vector>

namespace IBTK
{
/*!
 * InterfaceLocator is a class that can create an indicator function from a vector of meshes. The
 */
class InterfaceLocator
{
public:
    InterfaceLocator(std::string object_name,
                     SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > patch_hierarchy,
                     std::vector<FEDataManager*> fe_data_managers);

    virtual ~InterfaceLocator() = default;

    /*!
     * \brief Update the sign of the level set at points that are close to the boundary.
     */
    virtual void updateLocalLSSign(int ls_idx, double base_val, double time);

    /*!
     * \brief Update the sign of the level set at all Eulerian points.
     *
     * \note This uses a flooding algorithm, which can potentially require several sweeps through the computational
     * domain. This can be expensive and should only be done when necesary (e.g. regriddings, initializations).
     */
    virtual void resetGlobalLSSign(int ls_idx,
                                   SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > ls_var,
                                   double base_val,
                                   double time);

    inline void flipNormal(int part = 0)
    {
        TBOX_ASSERT(d_flip_normal.size() > part);
        d_flip_normal[part] = 1;
    }

    using BdryFcn = std::function<double(const IBTK::VectorNd&, double, void*)>;

    inline void registerBdryFcn(BdryFcn fcn, void* ctx = nullptr)
    {
        d_bdry_fcn = fcn;
        d_bdry_ctx = ctx;
    }

private:
    std::string d_object_name;
    std::vector<FEDataManager*> d_fe_data_managers;
    std::vector<int> d_flip_normal;
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > d_hierarchy;

    BdryFcn d_bdry_fcn;
    void* d_bdry_ctx = nullptr;
};
} // namespace IBTK

#endif
