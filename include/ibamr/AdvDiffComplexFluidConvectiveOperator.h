#ifndef included_AdvDiffOldroydBConvectiveOperator
#define included_AdvDiffOldroydBConvectiveOperator

#include <string>
#include <vector>

#include "CellVariable.h"
#include "CoarsenAlgorithm.h"
#include "FaceVariable.h"
#include "IntVector.h"
#include "PatchHierarchy.h"
#include "RefineAlgorithm.h"
#include "RefinePatchStrategy.h"
#include "ibamr/ConvectiveOperator.h"
#include "ibamr/ibamr_enums.h"
#include "tbox/Database.h"
#include "tbox/Pointer.h"
#include <ostream>
#include <stddef.h>
#include <string>
#include <vector>

#include "BoundaryBox.h"
#include "Box.h"
#include "CartesianGridGeometry.h"
#include "CartesianPatchGeometry.h"
#include "CellData.h"
#include "CellDataFactory.h"
#include "CellVariable.h"
#include "CoarsenAlgorithm.h"
#include "CoarsenOperator.h"
#include "CoarsenSchedule.h"
#include "FaceData.h"
#include "FaceVariable.h"
#include "IBAMR_config.h"
#include "Index.h"
#include "IntVector.h"
#include "MultiblockDataTranslator.h"
#include "Patch.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "RefineAlgorithm.h"
#include "RefineOperator.h"
#include "RefinePatchStrategy.h"
#include "RefineSchedule.h"
#include "SAMRAIVectorReal.h"
#include "Variable.h"
#include "VariableContext.h"
#include "VariableDatabase.h"
#include "ibamr/AdvDiffCenteredConvectiveOperator.h"
#include "ibamr/AdvDiffPPMConvectiveOperator.h"
#include "ibamr/AdvDiffPhysicalBoundaryUtilities.h"
#include "ibamr/AdvDiffWavePropConvectiveOperator.h"
#include "ibamr/ConvectiveOperator.h"
#include "ibamr/INSCollocatedHierarchyIntegrator.h"
#include "ibamr/INSHierarchyIntegrator.h"
#include "ibamr/INSStaggeredHierarchyIntegrator.h"
#include "ibamr/RHS_Operator.h"
#include "ibamr/ibamr_enums.h"
#include "ibamr/ibamr_utilities.h"
#include "ibamr/namespaces.h" // IWYU pragma: keep
#include "ibtk/CartExtrapPhysBdryOp.h"
#include "ibtk/CartSideRobinPhysBdryOp.h"
#include "ibtk/HierarchyGhostCellInterpolation.h"
#include "ibtk/PhysicalBoundaryUtilities.h"
#include "tbox/Array.h"
#include "tbox/Database.h"
#include "tbox/Pointer.h"
#include "tbox/Timer.h"
#include "tbox/TimerManager.h"
#include "tbox/Utilities.h"
#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <boost/concept_check.hpp>
#include <unsupported/Eigen/MatrixFunctions>

namespace SAMRAI
{
namespace solv
{
template <int DIM, class TYPE>
class SAMRAIVectorReal;
template <int DIM>
class RobinBcCoefStrategy;
} // namespace solv
namespace xfer
{
template <int DIM>
class CoarsenSchedule;
template <int DIM>
class RefineSchedule;
} // namespace xfer
} // namespace SAMRAI

namespace IBAMR
{
enum CFConvectiveOperatorType
{
    CENTERED,
    WAVE_PROP,
    CF_PPM,
    UNKNOWN = -1
};
template <>
inline std::string
enum_to_string<CFConvectiveOperatorType>(CFConvectiveOperatorType val)
{
    if (val == CENTERED) return "CENTERED";
    if (val == WAVE_PROP) return "WAVE_PROP";
    if (val == CF_PPM) return "PPM";
    return "UNKNOWN";
}
template <>
inline CFConvectiveOperatorType
string_to_enum<CFConvectiveOperatorType>(const std::string& str)
{
    if (strcasecmp(str.c_str(), "CENTERED") == 0) return CENTERED;
    if (strcasecmp(str.c_str(), "WAVE_PROP") == 0) return WAVE_PROP;
    if (strcasecmp(str.c_str(), "PPM") == 0) return CF_PPM;
    return UNKNOWN;
}
class AdvDiffComplexFluidConvectiveOperator : public ConvectiveOperator
{
public:
    // Constructor
    AdvDiffComplexFluidConvectiveOperator(const std::string& object_name,
                                          SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > Q_var,
                                          SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                                          CFConvectiveOperatorType difference_form,
                                          const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& conc_bc_coefs,
                                          SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM, double> > u_var);
    // Destructor
    ~AdvDiffComplexFluidConvectiveOperator();

    static SAMRAI::tbox::Pointer<ConvectiveOperator>
    allocate_operator(const std::string& object_name,
                      SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > Q_var,
                      SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                      CFConvectiveOperatorType difference_form,
                      const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& conc_bc_coefs,
                      SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM, double> > u_var)
    {
        return new AdvDiffComplexFluidConvectiveOperator(
            object_name, Q_var, input_db, difference_form, conc_bc_coefs, u_var);
    }
    void registerRHSTerm(Pointer<RHS_Operator> rhs);

    void registerRelaxationFunction(Pointer<CartGridFunction> fcn);

    void registerRateFunction(Pointer<CartGridFunction> fcn);

    void applyConvectiveOperator(int Q_idx, int Y_idx);

    void initializeOperatorState(const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& in,
                                 const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& out);

    void deallocateOperatorState();

private:
    AdvDiffComplexFluidConvectiveOperator();

    AdvDiffComplexFluidConvectiveOperator(const AdvDiffComplexFluidConvectiveOperator& from);

    AdvDiffComplexFluidConvectiveOperator& operator=(const AdvDiffComplexFluidConvectiveOperator& that);

    // Data communication algorithms, operators, and schedules.
    SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenAlgorithm<NDIM> > d_coarsen_alg_Q;
    SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenAlgorithm<NDIM> > d_coarsen_alg_u;
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenSchedule<NDIM> > > d_coarsen_scheds_Q;
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenSchedule<NDIM> > > d_coarsen_scheds_u;
    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineAlgorithm<NDIM> > d_ghostfill_alg_Q;
    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineAlgorithm<NDIM> > d_ghostfill_alg_u;
    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefinePatchStrategy<NDIM> > d_ghostfill_strategy_Q;
    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefinePatchStrategy<NDIM> > d_ghostfill_strategy_u;
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > > d_ghostfill_scheds_Q;
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > > d_ghostfill_scheds_u;
    std::string d_outflow_bdry_extrap_type;

    // Hierarchy configuration.
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > d_hierarchy;
    int d_coarsest_ln, d_finest_ln;

    // Scratch data.
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_Q_var;
    unsigned int d_Q_data_depth;
    int d_Q_scratch_idx;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM, double> > d_u_adv_var;
    int d_u_scratch_idx;

    // Convective Operator
    CFConvectiveOperatorType d_difference_form;
    SAMRAI::tbox::Pointer<IBAMR::ConvectiveOperator> d_convec_oper;
    int d_Q_convec_idx;
    const std::vector<RobinBcCoefStrategy<NDIM>*> d_conc_bc_coefs;
    Pointer<RHS_Operator> d_rhs;
    int d_R_idx;
    bool d_conform_tens, d_sqr_root, d_log_conform;
    double d_lambda, d_eta;
};
} // namespace IBAMR

#endif
