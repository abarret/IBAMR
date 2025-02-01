#ifndef included_IBAMR_IBRedundantInitializer_inl
#define included_IBAMR_IBRedundantInitializer_inl

#include <ibamr/IBRedundantInitializer.h>

namespace IBAMR
{
template <typename StreamableType>
inline void
IBRedundantInitializer::registerStreamable(StreamableInitFcn fcn, void* ctx)
{
    StreamableType::registerWithStreamableManager();
    d_streamables.push_back(StreamableData(fcn, ctx));
}
} // namespace IBAMR

#endif
