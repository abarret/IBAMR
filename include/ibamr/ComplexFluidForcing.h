
#ifndef included_ComplexFluidForcing
#define included_ComplexFluidForcing
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
#include "tbox/PIO.h"
#include "tbox/Pointer.h"
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
#include "ibamr/Giesekus_RHS.h"
#include "ibamr/INSStaggeredHierarchyIntegrator.h"
#include "ibamr/OldroydB_RHS.h"
#include "ibamr/RoliePoly_RHS.h"
#include "ibamr/ibamr_enums.h"
#include "ibtk/CartGridFunction.h"
#include "ibtk/muParserCartGridFunction.h"
#include "tbox/Array.h"
#include "tbox/Pointer.h"
#include <Eigen/Cholesky>
#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <boost/algorithm/string.hpp>
#include <mpi.h>
#include <petscsys.h>
#include <stddef.h>
#include <string>
#include <unsupported/Eigen/MatrixFunctions>
#include <vector>

// extern std::pair< iterator, iterator > r;
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
class ComplexFluidForcing : public IBTK::CartGridFunction
{
public:
    /*!
     * \brief This constructor creates Variable and VariableContext objects for storing the viscoleastic stresses at the
     * cneters of the Cartesian grid. Sets up the advection diffusion solver to use the velocity function prescribed.
     */
    ComplexFluidForcing(const std::string& object_name,
                        SAMRAI::tbox::Pointer<Database> input_db,
                        SAMRAI::tbox::Pointer<IBTK::CartGridFunction> u_fcn,
                        SAMRAI::tbox::Pointer<CartesianGridGeometry<NDIM> > grid_geometry,
                        SAMRAI::tbox::Pointer<IBAMR::AdvDiffSemiImplicitHierarchyIntegrator> adv_diff_integrator,
                        Pointer<VisItDataWriter<NDIM> > visit_data_writer);

    /*!
     * \brief This constructor creates Variable and VariableContext objects for
     * storing the stochastic stresses at the centers and nodes of the Cartesian
     * grid.
     */
    ComplexFluidForcing(const std::string& object_name,
                        SAMRAI::tbox::Pointer<Database> app_initializer,
                        const SAMRAI::tbox::Pointer<IBAMR::INSHierarchyIntegrator> fluid_solver,
                        SAMRAI::tbox::Pointer<CartesianGridGeometry<NDIM> > grid_geometry,
                        SAMRAI::tbox::Pointer<IBAMR::AdvDiffSemiImplicitHierarchyIntegrator> adv_diff_integrator,
                        Pointer<VisItDataWriter<NDIM> > visit_data_writer);

    SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > getVariable()
    {
        return d_W_cc_var;
    }
    int getStressTensorDrawVar()
    {
        return d_conform_idx_draw;
    }
    int getVariableIdx()
    {
        return d_W_cc_idx;
    }
    void registerRHSTerm(SAMRAI::tbox::Pointer<RHS_Operator> rhs);

    /*!
     * \brief Empty destructor.
     */
    ~ComplexFluidForcing();

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
     * \brief Evaluate the function on the patch interiors on the specified
     * levels of the patch hierarchy.
     */
    void setDataOnPatchHierarchy(const int data_idx,
                                 SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > var,
                                 SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                                 const double data_time,
                                 const bool initial_time = false,
                                 const int coarsest_ln = -1,
                                 const int finest_ln = -1);

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

    void setDataOnPatchLevel(const int data_idx,
                             SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > var,
                             SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level,
                             const double data_time,
                             const bool initial_time);
    //\}

    void checkPositiveDefiniteOnLevel(const int data_idx,
                                      const SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > var,
                                      const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level,
                                      const double data_time,
                                      const bool initial_time);
    void checkPositiveDefiniteOnPatch(const int data_idx,
                                      const SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > var,
                                      const SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                                      const double data_time,
                                      const bool initial_time);

    void applyGradientDetector(SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
                               int level_number,
                               double error_data_time,
                               int tag_index,
                               bool initial_time,
                               bool /*richardson_extrapolation_too*/);

    static void applyGradientDetectorCallback(SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
                                              int level_number,
                                              double error_data_time,
                                              int tag_index,
                                              bool initial_time,
                                              bool richardson_extrapolation_too,
                                              void* ctx);

    inline double getViscosity()
    {
        return d_eta;
    }

    inline double getRelaxationTime()
    {
        return d_lambda;
    }

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
    ComplexFluidForcing();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    ComplexFluidForcing(const ComplexFluidForcing& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    ComplexFluidForcing& operator=(const ComplexFluidForcing& that);

    /*!
     * \brief Get values from input file
     *
     * \param
     *
     */
    void getFromInput(const SAMRAI::tbox::Pointer<Database> input_db,
                      SAMRAI::tbox::Pointer<VisItDataWriter<NDIM> > visit_data_writer);

    void findDeterminant(const int data_idx,
                         const SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > var,
                         const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level,
                         const double data_time,
                         const bool initial_time);

    void squareMatrix(const int data_idx,
                      const SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > var,
                      const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                      const double data_time,
                      const bool initial_time,
                      const int coarsest_ln,
                      const int finest_ln);

    void exponentiateMatrix(const int data_idx,
                            const SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > var,
                            const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                            const double data_time,
                            const bool initial_time,
                            const int coarsest_ln,
                            const int finest_ln);

    void projectMatrix(const int data_idx,
                       const SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > var,
                       const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                       const double data_time,
                       const bool initial_time,
                       const int coarsest_ln,
                       const int finest_ln);
    bool d_project_conform;

    /*!
     * VariableContext and Variable objects for storing the components of the
     * stochastic stresses.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext> d_context;
    SAMRAI::tbox::Pointer<IBTK::muParserCartGridFunction> init_conds;
    AdvDiffSemiImplicitHierarchyIntegrator* const d_adv_diff_integrator;
    SAMRAI::tbox::Pointer<AdvDiffComplexFluidConvectiveOperator> d_convec_oper;
    IBAMR::CFConvectiveOperatorType d_convec_oper_type;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_W_cc_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_conform_var_draw, d_stress_var_draw,
        d_divW_var_draw;
    int d_conform_idx_draw, d_stress_idx_draw, d_divW_idx_draw;
    bool d_conform_draw, d_stress_draw, d_divW_draw;
    double d_lambda, d_eta;
    int d_W_cc_idx;
    int d_W_cc_scratch_idx;
    std::vector<RobinBcCoefStrategy<NDIM>*> d_conc_bc_coefs;
    std::string d_fluid_model;
    bool d_sqr_root, d_log_conform;
    SAMRAI::tbox::Pointer<PatchHierarchy<NDIM> > d_hierarchy;
    // Temp values for logging
    double d_max_det, d_min_det;
    bool d_log_det;
    bool d_positive_def, d_error_on_spd;
    double d_min_norm, d_max_norm;
    bool d_log_divW;

    // AMR stuffz
    SAMRAI::tbox::Array<double> d_divW_rel_thresh, d_divW_abs_thresh;
    bool d_divW_rel_tag, d_divW_abs_tag;

    // Velocity stuff
    SAMRAI::tbox::Pointer<IBTK::CartGridFunction> d_u_fcn;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM, double> > d_u_var;
}; // Private

} // Namespace IBAMR
#endif
