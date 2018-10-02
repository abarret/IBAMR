#include "ibamr/DetectCenterFcn.h"

// Namespace
namespace IBAMR
{
// Constructor
DetectCenterFcn::DetectCenterFcn(const std::string& object_name, Pointer<Database> /*input_db*/)
    : d_object_name(object_name)
{
    return;
} // End Constructor
// Destructor
DetectCenterFcn::~DetectCenterFcn()
{
    // intentionally blank
    return;
} // End Destructor

// Time Dependent?
bool
DetectCenterFcn::isTimeDependent() const
{
    return true;
} // End Time Dependent?

void
DetectCenterFcn::setDataOnPatch(const int data_idx,
                                Pointer<Variable<NDIM> > var,
                                Pointer<Patch<NDIM> > patch,
                                const double /*data_time*/,
                                const bool initial_time,
                                Pointer<PatchLevel<NDIM> > /*patch_level*/)
{
    const Box<NDIM>& patch_box = patch->getBox();
    const Pointer<CartesianPatchGeometry<NDIM> > p_geom = patch->getPatchGeometry();
    Pointer<CellData<NDIM, int> > ret_data = patch->getPatchData(data_idx);
#if !defined(NDEBUG)
    TBOX_ASSERT(ret_data);
    TBOX_ASSERT(ret_data->getDepth() == NDIM);
#endif
    ret_data->fillAll(0);
    // We need to detect if we are near the cylinder boundary and determine what kind of stencil we should use.
    // ret_data should be cell centered and integer valued with depth 2.
    // -1 -> standard stencil is fine.
    // 1st bit 1 -> dudx should be in positive x-direction.
    // 1st bit 0 -> dudx should be in negative x-direction.
    // 2nd bit 1 -> dudy should be in positive y-direction.
    // 2nd bit 0 -> dudy should be in negative y-direction.
    // 3rd bit 0 -> dudy can include center point
    if (initial_time) return;
    const double* const dx = p_geom->getDx();
    TBOX_ASSERT(MathUtilities<double>::equalEps(dx[0], dx[1]));
    double dX = dx[0];
    const double* const x_low = p_geom->getXLower();
    CellIndex<NDIM> l_idx = patch_box.lower();
    double r;
    std::vector<double> X(NDIM);
    for (CellIterator<NDIM> i(patch_box); i; i++)
    {
        CellIndex<NDIM> idx = i();
        r = 0.0;
        for (int d = 0; d < NDIM; ++d)
        {
            X[d] = x_low[d] + dx[d] * (idx(d) - l_idx(d) + 0.5);
            r += X[d] * X[d];
        }
        r = sqrt(r);
        if ((r < 1.0) || (r > (1.0 + dX)))
        {
            (*ret_data)(idx, 0) = (*ret_data)(idx, 1) = -1;
            continue;
        }
        for (int a = 0; a < NDIM; ++a)
        {
            for (int d = 0; d < NDIM; ++d)
            {
                // check in positive direction
                std::vector<double> Xn(X);
                Xn[d] += dX;
                r = sqrt(Xn[0] * Xn[0] + Xn[1] * Xn[1]);
                if (r < 1.0)
                {
                    (*ret_data)(idx, a) += 1 << d;
                }
            }
            // Now we check if we can use center point for stencils.
            int d = a == 0 ? 1 : 0;
            // Check off side point. Which side we need to choose is given by value in bit 2
            int side = (*ret_data)(idx, a) / 2 == 0 ? 1 : 0;
            std::vector<double> Xn(X);
            if (side == 1)
            {
                Xn[d] += 0.5 * dX;
            }
            else
            {
                Xn[d] -= 0.5 * dX;
            }
            r = sqrt(Xn[0] * Xn[0] + Xn[1] * Xn[1]);
            if (r < 1.0)
            {
                (*ret_data)(idx, a) += 8;
            }
        }
    }
} // end setDataOnPatch

} // namespace IBAMR
