#ifndef included_AdvDiffOldroydBConvectiveOperator
#define included_AdvDiffOldroydBConvectiveOperator

#include "IBAMR_config.h"

#include "ibamr/AdvDiffCenteredConvectiveOperator.h"
#include "ibamr/AdvDiffPPMConvectiveOperator.h"
#include "ibamr/AdvDiffPhysicalBoundaryUtilities.h"
#include "ibamr/AdvDiffWavePropConvectiveOperator.h"
#include "ibamr/CFRelaxationOperator.h"
#include "ibamr/ConvectiveOperator.h"
#include "ibamr/INSHierarchyIntegrator.h"
#include "ibamr/StaggeredStokesPhysicalBoundaryHelper.h"
#include "ibamr/ibamr_enums.h"
#include "ibamr/ibamr_utilities.h"
#include "ibamr/namespaces.h" // IWYU pragma: keep

#include "ibtk/CartExtrapPhysBdryOp.h"
#include "ibtk/HierarchyGhostCellInterpolation.h"
#include "ibtk/ibtk_utilities.h"

#include "Box.h"
#include "CartesianGridGeometry.h"
#include "CartesianPatchGeometry.h"
#include "CellData.h"
#include "CellVariable.h"
#include "FaceData.h"
#include "FaceVariable.h"
#include "Index.h"
#include "IntVector.h"
#include "Patch.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "SAMRAIVectorReal.h"
#include "Variable.h"
#include "VariableContext.h"
#include "VariableDatabase.h"
#include "tbox/Database.h"
#include "tbox/Pointer.h"
#include "tbox/Utilities.h"

#include <string>
#include <vector>

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
/*!
 * \brief Class AdvDiffComplexFluidConvectiveOperator is a concrete
 * ConvectiveOperator that implements a convective operator for Oldroyd-B type models.
 * A relaxation function for the stress must be registered before use of the operator.
 */

class AdvDiffComplexFluidConvectiveOperator : public ConvectiveOperator
{
public:
    // Constructor
    AdvDiffComplexFluidConvectiveOperator(const std::string& object_name,
                                          SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> Q_var,
                                          SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                                          CFConvectiveOperatorType difference_form,
                                          const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& conc_bc_coefs,
                                          const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& u_bc_coefs);
    // Destructor
    ~AdvDiffComplexFluidConvectiveOperator();

    static SAMRAI::tbox::Pointer<ConvectiveOperator>
    allocate_operator(const std::string& object_name,
                      SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> Q_var,
                      SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                      CFConvectiveOperatorType difference_form,
                      const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& conc_bc_coefs,
                      const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& u_bc_coefs)
    {
        return new AdvDiffComplexFluidConvectiveOperator(
            object_name, Q_var, input_db, difference_form, conc_bc_coefs, u_bc_coefs);
    }

    void registerRelaxationOperator(Pointer<CFRelaxationOperator> rhs);

    void applyConvectiveOperator(int Q_idx, int Y_idx);

    void initializeOperatorState(const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& in,
                                 const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& out);

    void deallocateOperatorState();

private:
    AdvDiffComplexFluidConvectiveOperator();

    AdvDiffComplexFluidConvectiveOperator(const AdvDiffComplexFluidConvectiveOperator& from);

    AdvDiffComplexFluidConvectiveOperator& operator=(const AdvDiffComplexFluidConvectiveOperator& that);

    // Hierarchy configuration.
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> d_hierarchy;
    int d_coarsest_ln = IBTK::invalid_level_number, d_finest_ln = IBTK::invalid_level_number;

    // Scratch data.
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> d_Q_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double>> d_u_adv_var;
    int d_u_scratch_idx = IBTK::invalid_index;

    // Convective Operator
    CFConvectiveOperatorType d_difference_form;
    SAMRAI::tbox::Pointer<IBAMR::ConvectiveOperator> d_convec_oper;
    int d_Q_convec_idx = IBTK::invalid_index;
    const std::vector<RobinBcCoefStrategy<NDIM>*> d_conc_bc_coefs;
    const std::vector<RobinBcCoefStrategy<NDIM>*> d_u_bc_coefs;
    Pointer<CFRelaxationOperator> d_rhs;
    int d_R_idx = IBTK::invalid_index;
    bool d_conform_tens = true, d_sqr_root = false, d_log_conform = false;
    std::string d_interp_type = "LINEAR";
};
} // namespace IBAMR

#endif
