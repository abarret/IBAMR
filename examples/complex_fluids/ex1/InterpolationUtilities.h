#ifndef included_InterpolationUtilities
#define included_InterpolationUtilities

#include "IBAMR_config.h"
#include "IBTK_config.h"

#include "ibamr/namespaces.h"
#include <ibamr/app_namespaces.h>

#include "ibtk/HierarchyGhostCellInterpolation.h"
#include "ibtk/IndexUtilities.h"

#include "Box.h"
#include "CartesianGridGeometry.h"
#include "CartesianPatchGeometry.h"
#include "CellIndex.h"
#include "Index.h"
#include "IntVector.h"
#include "RobinBcCoefStrategy.h"
#include "SAMRAI_config.h"
#include "tbox/Pointer.h"
#include <tbox/SAMRAI_MPI.h>

#include <Patch.h>
#include <PatchLevel.h>

#include <vector>

namespace SAMRAI
{
namespace hier
{
template <int DIM>
class Patch;
template <int DIM>
class PatchHierarchy;
template <int DIM>
class Variable;
} // namespace hier
} // namespace SAMRAI

namespace IBTK
{
class InterpolationUtilities
{
public:
    // Use polynomial interpolation
    static double interpolate(const std::vector<double>& X,
                              const int data_idx,
                              SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> Q_var,
                              SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> patch_hierarchy,
                              const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& bc_coefs,
                              const double data_time,
                              const int depth = 0);
    // Use least squares interpolation
    static double interpolateL2(const std::vector<double>& X,
                                const int data_idx,
                                SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> Q_var,
                                SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> patch_hierarchy,
                                const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& bc_coefs,
                                const double data_time,
                                const int depth = 0);
    static double weightFcn(const std::vector<double>&, const std::vector<double>&);

private:
    static double interpolate(const double& x, const std::vector<int>& xi, const std::vector<double>& yi);

    static double interpolateInBoxes(const CellIndex<NDIM>& idx,
                                     const std::vector<double>& X,
                                     SAMRAI::pdat::CellData<NDIM, int>& r_data,
                                     SAMRAI::pdat::CellData<NDIM, double>& q_data,
                                     Pointer<CartesianPatchGeometry<NDIM>> pgeom,
                                     const Box<NDIM>& pbox,
                                     int dim,
                                     int cycle,
                                     std::vector<int>& completed_dims);
};

} // namespace IBTK
#endif
