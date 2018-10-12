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
    return false;
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
    // 1st trit for dudx:
    //    0 : standard stencil
    //    1 : positive biased
    //    2 : negative biased
    // 2nd trit for dudy:
    //    0 : standard stencil
    //    1 : positive biased
    //    2 : negative biased
    // 3rd trit for interpolation
    //    0 : standard interpolation
    //    1 : positive biased interp
    //    2 : negative biased interp
    if (initial_time) return;
    const double* const dx = p_geom->getDx();
    TBOX_ASSERT(MathUtilities<double>::equalEps(dx[0], dx[1]));
    double dX = dx[0];
    const double* const x_low = p_geom->getXLower();
    CellIndex<NDIM> l_idx = patch_box.lower();
    double r, rx;
    std::vector<double> X(NDIM);
    for (CellIterator<NDIM> i(patch_box); i; i++)
    {
        CellIndex<NDIM> idx = i();
        rx = 0.0;
        for (int d = 0; d < NDIM; ++d)
        {
            X[d] = x_low[d] + dx[d] * (idx(d) - l_idx(d) + 0.5);
            rx += X[d] * X[d];
        }
        rx = sqrt(rx);
        if ((rx < 1.0 - dX) || (rx > (1.0 + dX)))
        {
            (*ret_data)(idx, 0) = (*ret_data)(idx, 1) = 0;
            continue;
        }
        for (int a = 0; a < NDIM; ++a)
        {
            for (int d = 0; d < NDIM; ++d)
            {
                // We have a compact stencil, things are easy.
                // check in positive direction
                std::vector<double> Xn(X);
                Xn[d] += dX;
                r = sqrt(Xn[0] * Xn[0] + Xn[1] * Xn[1]);
                if ((r < 1.0 && rx > 1.0) || (r > 1.0 && rx < 1.0))
                {
                    // Stencil crosses a boundary. We need to go in negative direction.
                    (*ret_data)(idx, a) += 2 * std::pow(3, d);
//                     if (d != a)
//                     {
//                         // Stencil is more spread out, we need to be more careful.
//                         // Check negative side.
//                         Xn[(d == 0 ? 1 : 0)] -= dX;
//                         double rr;
//                         rr = sqrt(Xn[0]* Xn[0] + Xn[1]*Xn[1]);
//                         if ((rr > 1.0 && r < 1.0) || (rr < 1.0 && r > 1.0))
//                         {
//                             // Negative point crosses a boundary. Need to go in positive side.
//                             (*ret_data)(idx,a) += std::pow(3,3);
//                         }
//                         Xn[(d == 0 ? 1 : 0)] += 1.5*dX;
//                         rr = sqrt(Xn[0] * Xn[0] + Xn[1] * Xn[1]);
//                         if ((rr > 1.0 && r < 1.0) || (rr < 1.0 && r > 1.0))
//                         {
//                             // Positive point crosses a boundary. Need to go in negative side.
//                             (*ret_data)(idx,a) += 2*std::pow(3,3);
//                         }
//                     }
                }
                else
                {
                    Xn = X;
                    Xn[d] -= dX;
                    r = sqrt(Xn[0] * Xn[0] + Xn[1] * Xn[1]);
                    if ((r > 1.0 && rx < 1.0) || (r < 1.0 && rx > 1.0))
                    {
                        // Stencil crosses a boundary. We need to go in positive direction.
                        (*ret_data)(idx, a) += std::pow(3, d);
//                         if (d != a)
//                         {
//                             // Stencil is more spread out, we need to be more careful.
//                             // Check negative side.
//                             Xn[(d == 0 ? 1 : 0)] -= dX;
//                             double rr;
//                             rr = sqrt(Xn[0]* Xn[0] + Xn[1]*Xn[1]);
//                             if ((rr > 1.0 && r < 1.0) || (rr < 1.0 && r > 1.0))
//                             {
//                                 // Negative point crosses a boundary. Need to go in positive side.
//                                 (*ret_data)(idx,a) += std::pow(3,3);
//                             }
//                             Xn[(d == 0 ? 1 : 0)] += 1.5*dX;
//                             rr = sqrt(Xn[0] * Xn[0] + Xn[1] * Xn[1]);
//                             if ((rr > 1.0 && r < 1.0) || (rr < 1.0 && r > 1.0))
//                             {
//                                 // Positive point crosses a boundary. Need to go in negative side.
//                                 (*ret_data)(idx,a) += 2*std::pow(3,3);
//                             }
//                         }
                    }
                }
            }
        }
    }
} // end setDataOnPatch

} // namespace IBAMR
