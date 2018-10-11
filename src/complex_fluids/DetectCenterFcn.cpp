#include "ibamr/DetectCenterFcn.h"

// Namespace
namespace IBAMR
{
// Constructor
DetectCenterFcn::DetectCenterFcn(const std::string& object_name, Pointer<Database> /*input_db*/)
    : d_object_name(object_name),
      d_cell_var(new CellVariable<NDIM, bool>(object_name + "::RCell")),
      d_side_var(new SideVariable<NDIM, bool>(object_name + "::RSide")),
      d_cell_idx(-1),
      d_side_idx(-1),
      d_set(false)
{
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<VariableContext> ctx = var_db->getContext(object_name + "::Context");
    d_cell_idx = var_db->registerVariableAndContext(d_cell_var, ctx);
    d_side_idx = var_db->registerVariableAndContext(d_side_var, ctx);

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
DetectCenterFcn::setDataOnPatchHierarchy(const int data_idx,
                                         Pointer<Variable<NDIM> > var,
                                         Pointer<Patch<NDIM> > patch,
                                         const double data_time,
                                         const bool initial_time,
                                         Pointer<PatchHierarchy<NDIM> > patch_hierarchy)
{
    int coarsest_ln = 0;
    int finest_ln = patch_hierarchy->getFinestLevelNumber();

    for (int ln = coarsest_ln; ln < finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
        if (!level->checkAllocated(d_cell_idx)) level->allocatePatchData(d_cell_idx);
        if (!level->checkAllocated(d_side_idx)) level->allocatePatchData(d_side_idx);
    }
    if (!d_set)
    {
        for (int ln = coarsest_ln; ln < finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM> > patch = level->getPatch(p());
                const Box<NDIM>& box = patch->getBox();
                Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
                const double* const dx = patch_geom->getDx();
                const double* const xlow = patch_geom->getXLower();
                const CellIndex<NDIM>& low_idx_c = box.lower();
                Pointer<CellData<NDIM, bool> > cell_data = patch->getPatchData(d_cell_idx);
                for (CellIterator<NDIM> ci(box); ci; ci++)
                {
                    CellIndex idx = ci();
                    std::vector<double> X(NDIM);
                    for (int d = 0; d < NDIM; ++d) X[d] = xlow[d] + (idx(d) - low_idx_c(d) + 0.5) * dx[d];
                    (*cell_data)(idx) = sqrt(X[0] * X[0] + X[1] * X[1]) > 1.0 ? true : false;
                }
                Pointer<SideData<NDIM, bool> > side_data = patch->getPatchData(d_side_idx);
                const SideIndex<NDIM>& low_idx_s = box.lower();
                for (int axis = 0; axis < NDIM; ++axis)
                {
                    for (SideIterator<NDIM> si(box, axis); si; si++)
                    {
                        SideIndex<NDIM> idx = si();
                        std::vector<double> X(NDIM);
                        for (int d = 0; d < NDIM; ++d)
                            X[d] = xlow[d] + (idx(d) - low_idx_s(d) + (axis == d ? 0.0 : 0.5)) * dx[d];
                        (*side_data)(idx) = sqrt(X[0] * X[0] + X[1] * X[1]) > 1.0 ? true : false;
                    }
                }
            }
        }
        d_set = true;
    }
    for (int ln = coarsest_ln; ln < finest_ln; ++ln)
    {
        setDataOnPatchLevel(data_idx, var, patch, data_time, initial_time, patch_hierarchy->getPatchLevel(ln));
    }
    return;
}

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
                        continue;
                    }
                }
                if (d != a)
                {
                    // The stencil is spread out, we need to be more careful.
                    // check in positive direction
                    std::vector<double> Xn(X);
                    double rl, rr;
                    Xn[d] += dX;
                    r = sqrt(Xn[0] * Xn[0] + Xn[1] * Xn[1]);
                    if ((r < 1.0 && rx > 1.0) || (r > 1.0 && rx < 1.0))
                    {
                        // Stencil crosses a boundary. We need to go in negative direction.
                        // Check left and right point to see one sided interpolation.
                        // Left
                        Xn[d % a] -= 0.5 * dX;
                        rl = sqrt(Xn[0] * Xn[0] + Xn[1] * Xn[1]);
                        if ((rl > 1.0 && r < 1.0) || (rl < 1.0 && r > 1.0))
                        {
                            // Left point crosses a boundary.
                            (*ret_data)(idx, a) += std::pow(3, 3);
                        }
                        // Right
                        Xn[d % a] += 1.5 * dX;
                        rr = sqrt(Xn[0] * Xn[0] + Xn[1] * Xn[1]);
                        if ((rr > 1.0 && r > 1.0) || (rr < 1.0 && r > 1.0))
                        {
                            // Right point crosses a boundary.
                            (*ret_data)(idx, a) += 2 * std::pow(3, 3);
                        }
                    }
                    Xn = X;
                    Xn[d] -= dX;
                    r = sqrt(Xn[0] * Xn[0] + Xn[1] * Xn[1]);
                    if ((r < 1.0 && rx > 1.0) || (r > 1.0 && rx < 1.0))
                    {
                        // Stencil crosses a boundary. We need to go in positive direction
                    }
            }
        }
    }
} // end setDataOnPatch

} // namespace IBAMR
