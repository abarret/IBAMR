#ifndef included_IBTK_ls_utilities
#define included_IBTK_ls_utilities

#include <tbox/Pointer.h>

#include <PatchHierarchy.h>
#include <PatchLevel.h>
#include <Variable.h>

namespace IBTK
{
void flood_fill_on_level_cell(const int ls_idx,
                              SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level,
                              double base_val,
                              double time);
void flood_fill_on_level_side(int ls_idx,
                              SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level,
                              double base_val,
                              double time);

} // namespace IBTK
#endif
