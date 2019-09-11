#include "ibamr/CFRelaxationOperator.h"

// Namespace
namespace IBAMR
{
// Constructor
CFRelaxationOperator::CFRelaxationOperator(const std::string& object_name) : d_object_name(object_name)
{
    // intentionally blank
    return;
} // End Constructor

bool
CFRelaxationOperator::isTimeDependent() const
{
    return true;
}

void
CFRelaxationOperator::setPatchDataIndex(const int data_idx)
{
    d_W_cc_idx = data_idx;
}

} // namespace IBAMR
