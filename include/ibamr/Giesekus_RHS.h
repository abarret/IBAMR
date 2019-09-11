
#ifndef included_Giesekus_RHS
#define included_Giesekus_RHS
#include "ibamr/AdvDiffComplexFluidConvectiveOperator.h"
#include <ibamr/app_namespaces.h>
#include <ibtk/muParserRobinBcCoefs.h>
#include <limits>
#include <math.h>
#include <ostream>
#include <stddef.h>
#include <string>
#include <vector>

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
#include "ibamr/RHS_Operator.h"
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
class Giesekus_RHS : public IBAMR::RHS_Operator
{
public:
    /*!
     * \brief This constructor creates Variable and VariableContext objects for
     * storing the stochastic stresses at the centers and nodes of the Cartesian
     * grid.
     */
    Giesekus_RHS(const std::string& object_name,
                 Pointer<Database> input_db,
                 Pointer<CellVariable<NDIM, double> > Q_var,
                 Pointer<AdvDiffSemiImplicitHierarchyIntegrator> adv_diff_integrator);

    /*!
     * \brief Empty destructor.
     */
    ~Giesekus_RHS();

    /*!
     * \name Methods to set patch data.
     */
    //\{

    /*!
     * \brief Indicates whether the concrete INSStaggeredStochasticForcing object is
     * time-dependent.
     */
    bool isTimeDependent() const;

    /*!
     * \brief Evaluate the function on the patch interior.
     */
    void setDataOnPatch(const int data_idx,
                        SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > var,
                        SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                        const double data_time,
                        const bool initial_time = false,
                        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > patch_level =
                            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> >(NULL));

    //\}

protected:
    /*!
     * The object name is used for error/warning reporting.
     */
    std::string d_object_name;

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    Giesekus_RHS();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    Giesekus_RHS(const Giesekus_RHS& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    Giesekus_RHS& operator=(const Giesekus_RHS& that);
    Pointer<CellVariable<NDIM, double> > d_W_cc_var;
    double d_alpha;
    const AdvDiffSemiImplicitHierarchyIntegrator* const d_adv_diff_integrator;
    bool d_sqr_root, d_log_conform;
    double d_lambda;
}; // Private

} // Namespace IBAMR
#endif
