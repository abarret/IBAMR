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
    return;
}

void
RHS_Operator::setPatchDataIndex(const int data_idx)
{
    d_W_cc_idx = data_idx;
}

} // namespace IBAMR
