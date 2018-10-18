#include "InterpolationUtilities.h"
#include "tbox/Pointer.h"
#include <HierarchyCellDataOpsReal.h>
#include <Patch.h>
#include <PatchLevel.h>
#include <Eigen/Core>
#include <Eigen/QR>
#include <boost/concept_check.hpp>

using Eigen::Matrix3d;
using Eigen::Vector3d;

namespace IBTK
{
double
InterpolationUtilities::interpolate(const std::vector<double>& X,
                                    const int data_idx,
                                    Pointer<CellVariable<NDIM, double> > Q_var,
                                    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > patch_hierarchy,
                                    const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& bc_coefs,
                                    const double data_time,
                                    const int depth)
{
    double q_val = 0.0;
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    const int data_idx_temp =
        var_db->registerVariableAndContext(Q_var, var_db->getContext("Interpolation"), IntVector<NDIM>(3));
    for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
        level->allocatePatchData(data_idx_temp);
    }
    HierarchyCellDataOpsReal<NDIM, double> hier_cc_data_ops(patch_hierarchy);
    hier_cc_data_ops.copyData(data_idx_temp, data_idx);
    typedef HierarchyGhostCellInterpolation::InterpolationTransactionComponent InterpolationTransactionComponent;
    std::vector<InterpolationTransactionComponent> ghost_cell_components(1);
    ghost_cell_components[0] = InterpolationTransactionComponent(
        data_idx_temp, "CONSERVATIVE_LINEAR_REFINE", false, "CONSERVATIVE_COARSEN", "LINEAR", false, bc_coefs, NULL);
    HierarchyGhostCellInterpolation ghost_fill_op;
    ghost_fill_op.initializeOperatorState(ghost_cell_components, patch_hierarchy);
    ghost_fill_op.fillData(data_time);
    bool done = false;
    for (int ln = patch_hierarchy->getFinestLevelNumber(); ln >= 0 && !done; --ln)
    {
        // Start at the finest level...
        Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
        CellIndex<NDIM> idx = IndexUtilities::getCellIndex(X, level->getGridGeometry(), level->getRatio());
        for (PatchLevel<NDIM>::Iterator p(level); p && !done; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Pointer<CartesianPatchGeometry<NDIM> > p_geom = patch->getPatchGeometry();
            const double* const dx = p_geom->getDx();
            const double* const x_lower = p_geom->getXLower();
            const Box<NDIM>& patch_box = patch->getBox();
            Pointer<CellData<NDIM, double> > S_data = patch->getPatchData(data_idx_temp);
            if (patch_box.contains(idx))
            {
                // Great. The patch is currently on this level
                // Let's create a box that contains this data
                Box<NDIM> box(idx, idx);
                // Grow it by some number of grid cells
                box.grow(2);
                // Loop through the box, make sure the point is located
                // OUTSIDE the disk
                double r;
                std::vector<double> x(NDIM);
                CellData<NDIM, int> i_data(box, 1, IntVector<NDIM>(0));
                const CellIndex<NDIM> ci_l = patch_box.lower();
                int num = 0;
                for (CellIterator<NDIM> i(box); i; i++)
                {
                    CellIndex<NDIM> ci = i();
                    for (int d = 0; d < NDIM; ++d) x[d] = x_lower[d] + dx[d] * (ci(d) - ci_l(d) + 0.5);
                    double r = sqrt(x[0] * x[0] + x[1] * x[1]);
                   if (r > 1.0)
                   {
                        i_data(ci) = 1;
                        num++;
                   }
                }
                // Moving least squares
                VectorXd rhs = VectorXd::Zero(6);
                VectorXd soln = VectorXd::Zero(6);
                MatrixXd mat = MatrixXd::Zero(6,6);
                for (CellIterator<NDIM> i(box); i; i++)
                {
                    CellIndex<NDIM> ci = i();
                    if (i_data(ci) == 1)
                    {
                        for (int d = 0; d < NDIM; ++d) x[d] = x_lower[d] + dx[d] * (ci(d) - ci_l(d)+0.5);
                        // Fill in RHS and MATRIX values
                        double w = weightFcn(X, x);
                        double f = (*S_data)(ci, depth);
                        rhs(0) += w*f;
                        rhs(1) += w*f*x[0];
                        rhs(2) += w*f*x[1];
                        rhs(3) += w*f*x[0]*x[0];
                        rhs(4) += w*f*x[1]*x[1];
                        rhs(5) += w*f*x[0]*x[1];
                        mat(0,0) += w; /*mat(0,1) += x[0]*w; mat(0,2) += x[1]*w; mat(0,3) += x[0]*x[0]*w; mat(0,4) += x[1]*x[1]*w; mat(0,5) += x[0]*x[1]*w;*/
                        mat(1,0) += x[0]*w; mat(1,1) += x[0]*x[0]*w;/* mat(1,2) += x[0]*x[1]*w; mat(2,3) += x[0]*x[0]*x[0]*w; mat(2,4) += x[0]*x[1]*x[1]*w; mat(3,5) += x[0]*x[0]*x[1]*w;*/
                        mat(2,0) += x[1]*w; mat(2,1) += x[1]*x[0]*w; mat(2,2) += x[1]*x[1]*w;
                        mat(3,0) += x[0]*x[0]*w; mat(3,1) += x[0]*x[0]*x[0]*w; mat(3,2) += x[0]*x[0]*x[1]*w; mat(3,3) += x[0]*x[0]*x[0]*x[0]*w;
                        mat(4,0) += x[1]*x[1]*w; mat(4,1) += x[1]*x[1]*x[0]*w; mat(4,2) += x[1]*x[1]*x[1]*w; mat(4,3) += x[1]*x[1]*x[0]*x[0]*w; mat(4,4) += x[1]*x[1]*x[1]*x[1]*w;
                        mat(5,0) += x[0]*x[1]*w; mat(5,1) += x[0]*x[1]*x[0]*w; mat(5,2) += x[0]*x[1]*x[1]*w; mat(5,3) += x[0]*x[1]*x[0]*x[0]*w; mat(5,4) += x[0]*x[1]*x[1]*x[1]*w; mat(5,5) += x[0]*x[1]*x[0]*x[1]*w;
                    }
                }
//                mat(0,1) = mat(1,0); mat(0,2) = mat(2,0); mat(0,3) = mat(3,0); mat(0,4) = mat(4,0); mat(0,5) = mat(5,0);
//                mat(1,2) = mat(2,1); mat(1,3) = mat(3,1); mat(1,4) = mat(4,1); mat(1,5) = mat(5,1);
//                mat(2,3) = mat(3,2); mat(2,4) = mat(4,2); mat(2,5) = mat(5,2);
//                mat(3,4) = mat(4,3); mat(3,5) = mat(5,3);
//                mat(4,5) = mat(5,4);
                soln = mat.ldlt().solve(rhs);
                q_val = soln(0) + soln(1) * X[0] + soln(2) * X[1] + soln(3) * X[0]*X[0] + soln(4) * X[1]*X[1] + soln(5) * X[0]*X[1];
                done = true;
            }
        }
    }
    q_val = SAMRAI_MPI::sumReduction(q_val);

    for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
        level->deallocatePatchData(data_idx_temp);
    }
    return q_val;
}

double InterpolationUtilities::weightFcn(const std::vector<double>& x, const std::vector<double>& x_j)
{
    double r = sqrt((x[0]-x_j[0])*(x[0]-x_j[0])+(x[1]-x_j[1])*(x[1]-x_j[1]));
    double w = exp(-10.0*r*r);
    w = 1.0;
    return w;
}

} // namespace IBTK
