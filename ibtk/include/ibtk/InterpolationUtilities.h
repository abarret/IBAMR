#ifndef included_InterpolationUtilities
#define included_InterpolationUtilities

#include "IBTK_config.h"
#include "SAMRAI_config.h"
#include <vector>

#include "Box.h"
#include "CartesianGridGeometry.h"
#include "CartesianPatchGeometry.h"
#include "CellIndex.h"
#include "Index.h"
#include "IntVector.h"
#include "RobinBcCoefStrategy.h"
#include "ibtk/HierarchyGhostCellInterpolation.h"
#include "tbox/Pointer.h"
#include <Patch.h>
#include <PatchLevel.h>
#include <tbox/SAMRAI_MPI.h>

#include "ibtk/IndexUtilities.h"

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
    /*
     * Does a weighted least squares quadratic approximation to the cell centered data in data_idx at the point
     * specified in X.
     * Uses the weights as calculated in weightFcn.
     */
    static double interpolate(const std::vector<double>& X,
                              const int data_idx,
                              SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > Q_var,
                              SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > patch_hierarchy,
                              const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& bc_coefs,
                              const double data_time,
                              const int depth = 0);
    /*
     * Weight function for weighted least squares.
     */
    static double weightFcn(const std::vector<double>&, const std::vector<double>&);

    /*
     * Does a weighted least squares quadratic approximation to the gradient of the side centered data in data_idx at
     * the point specified in X.
     * Uses the weights as calculated in weightFcn.
     */
    static void gradient(std::vector<double>& ret,
                         const std::vector<double>& Q,
                         const int data_idx,
                         SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > Q_var,
                         SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > patch_hierarchy,
                         const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& bc_coefs,
                         const double data_time,
                         const int depth = 0);

    /*
     * Does a weighted least squares quadratic approximation to the gradient of the side centered data in data_idx
     * Uses the weights as calculated in weightFcn.
     */
    static void gradient(const int grad_idx,
                         const int data_idx,
                         SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > Q_var,
                         SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > patch_hierarchy,
                         const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& bc_coefs,
                         const double data_time,
                         const int depth = 0);
};

} // namespace IBTK
#endif
