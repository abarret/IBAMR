#include "ibamr/RHS_Operator.h"

// Namespace
namespace IBAMR
{
// Constructor
RHS_Operator::RHS_Operator(const std::string& object_name) : d_object_name(object_name), d_W_cc_idx(-1)
{
    // intentionally blank
    return;
} // End Constructor

RHS_Operator::~RHS_Operator()
{
    // intentionally blank
    return;
}

void
RHS_Operator::setDataOnPatchHierarchy(const int data_idx,
                                      Pointer<Variable<NDIM> > var,
                                      Pointer<PatchHierarchy<NDIM> > hierarchy,
                                      const double data_time,
                                      const bool initial_time,
                                      const int coarsest_ln_in,
                                      const int finest_ln_in)
{
    const int coarsest_ln = (coarsest_ln_in == -1 ? 0 : coarsest_ln_in);
    const int finest_ln = (finest_ln_in == -1 ? hierarchy->getFinestLevelNumber() : finest_ln_in);
    // Compute RHS on each patch level
    for (int level_num = coarsest_ln; level_num <= finest_ln; ++level_num)
    {
        setDataOnPatchLevel(data_idx, var, hierarchy->getPatchLevel(level_num), data_time, initial_time);
    }
    return;
} // End setDataOnPatchHierarchy

void
RHS_Operator::setPatchDataIndex(const int data_idx)
{
    d_W_cc_idx = data_idx;
}

} // namespace IBAMR
