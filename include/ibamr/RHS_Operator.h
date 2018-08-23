
#ifndef included_RHS_Operator
#define included_RHS_Operator
#include <limits>
#include <math.h>
#include <ostream>
#include <stddef.h>
#include <string>
#include <vector>
//#include "AdvDiffComplexFluidConvectiveOperator.h"
#include <ibamr/app_namespaces.h>
#include <ibtk/muParserRobinBcCoefs.h>

#include "ArrayData.h"
#include "BoundaryBox.h"
#include "Box.h"
#include "CartesianPatchGeometry.h"
#include "CellData.h"
#include "CellIndex.h"
#include "CellVariable.h"
#include "EdgeData.h"     // IWYU pragma: keep
#include "EdgeGeometry.h" // IWYU pragma: keep
#include "EdgeIndex.h"    // IWYU pragma: keep
#include "HierarchyDataOpsManager.h"
#include "HierarchyDataOpsReal.h"
#include "IBAMR_config.h"
#include "Index.h"
#include "IntVector.h"
#include "LocationIndexRobinBcCoefs.h"
#include "MultiblockDataTranslator.h"
#include "NodeData.h"     // IWYU pragma: keep
#include "NodeGeometry.h" // IWYU pragma: keep
#include "NodeIndex.h"    // IWYU pragma: keep
#include "NodeVariable.h"
#include "Patch.h"
#include "PatchGeometry.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "RobinBcCoefStrategy.h"
#include "SideData.h"
#include "Variable.h"
#include "VariableContext.h"
#include "VariableDatabase.h"
#include "ibamr/AdvDiffSemiImplicitHierarchyIntegrator.h"
#include "ibamr/INSStaggeredHierarchyIntegrator.h"
#include "ibamr/RNG.h"
#include "ibamr/StokesSpecifications.h"
#include "ibamr/ibamr_enums.h"
#include "ibamr/namespaces.h" // IWYU pragma: keep
#include "ibtk/HierarchyGhostCellInterpolation.h"
#include "ibtk/PhysicalBoundaryUtilities.h"
#include "tbox/Array.h"
#include "tbox/Database.h"
#include "tbox/Pointer.h"
#include "tbox/Utilities.h"
#if (NDIM == 2)
#include "ibtk/NodeDataSynchronization.h"
#endif
#if (NDIM == 3)
#include "ibtk/EdgeDataSynchronization.h"
#endif
#include "CellVariable.h"
#include "EdgeVariable.h" // IWYU pragma: keep
#include "IntVector.h"
#include "NodeVariable.h" // IWYU pragma: keep
#include "PatchLevel.h"
#include "VariableContext.h"
#include "VisItDataWriter.h"
#include "ibamr/INSStaggeredHierarchyIntegrator.h"
#include "ibamr/ibamr_enums.h"
#include "ibtk/CartGridFunction.h"
#include "ibtk/muParserCartGridFunction.h"
#include "tbox/Array.h"
#include "tbox/Pointer.h"
#include <stddef.h>
#include <string>
#include <unsupported/Eigen/MatrixFunctions>
#include <vector>
namespace IBAMR
{
class ComplexFluidForcing;
} // namespace IBAMR
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
namespace tbox
{
class Database;
} // namespace tbox
} // namespace SAMRAI
/////////////////////////////// INCLUDES /////////////////////////////////////

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class INSStaggeredStochasticForcing provides an interface for
 * specifying a stochastic forcing term for a staggered-grid incompressible
 * Navier-Stokes solver.
 */
class RHS_Operator : public IBTK::CartGridFunction
{
public:
    /*!
     * \brief This constructor creates Variable and VariableContext objects for
     * storing the stochastic stresses at the centers and nodes of the Cartesian
     * grid.
     */
    RHS_Operator(const std::string& object_name);

    /*!
     * \brief Empty destructor.
     */
    virtual ~RHS_Operator();

    /*!
     * \name Methods to set patch data.
     */
    //\{

    /*!
     * \brief Indicates whether the concrete INSStaggeredStochasticForcing object is
     * time-dependent.
     */
    virtual bool isTimeDependent() const = 0;

    /*!
     * \brief Evaluate the function on the patch interiors on the specified
     * levels of the patch hierarchy.
     */
    virtual void setDataOnPatchHierarchy(const int data_idx,
                                         SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > var,
                                         SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                                         const double data_time,
                                         const bool initial_time = false,
                                         const int coarsest_ln = -1,
                                         const int finest_ln = -1);

    /*!
     * \brief Evaluate the function on the patch interior.
     */
    virtual void setDataOnPatch(const int data_idx,
                                SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > var,
                                SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                                const double data_time,
                                const bool initial_time = false,
                                SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > patch_level =
                                    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> >(NULL)) = 0;

    void setPatchDataIndex(const int data_idx);

    //\}

protected:
    /*!
     * The object name is used for error/warning reporting.
     */
    std::string d_object_name;
    double d_W_cc_idx;

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    RHS_Operator();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    RHS_Operator(const RHS_Operator& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    RHS_Operator& operator=(const RHS_Operator& that);
}; // Private

} // Namespace IBAMR
#endif
