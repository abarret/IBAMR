#include "InterpolationUtilities.h"
#include "tbox/Pointer.h"
#include <HierarchyCellDataOpsReal.h>
#include <Patch.h>
#include <PatchLevel.h>
//#include <Eigen/src/Core/Matrix.h>
#include <Eigen/Core>
//#include <Eigen/src/QR/ColPivHouseholderQR.h>
#include <Eigen/QR>

using Eigen::Matrix3d;
using Eigen::Vector3d;
// extern "C" {
//     void weno4_(const double*,
//                 const double&,
//                 const double&);
// }

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
                for (CellIterator<NDIM> i(box); i; i++)
                {
                    CellIndex<NDIM> ci = i();
                    for (int d = 0; d < NDIM; ++d) x[d] = x_lower[d] + dx[d] * (ci(d) - ci_l(d) + 0.5);
                    double r = sqrt(x[0] * x[0] + x[1] * x[1]);
                    if (r > 1.0) i_data(ci) = 1;
                }
                // Solve the least squares problem, and we're done!
                Vector3d rhs = Vector3d::Zero();
                Vector3d soln = Vector3d::Zero();
                Matrix3d mat = Matrix3d::Zero();
                for (CellIterator<NDIM> i(box); i; i++)
                {
                    CellIndex<NDIM> ci = i();
                    if (i_data(ci) == 1)
                    {
                        for (int d = 0; d < NDIM; ++d) x[d] = x_lower[d] + dx[d] * (ci(d) - ci_l(d) + 0.5);
                        // Fill in RHS and MATRIX values
                        rhs(0) += (*S_data)(ci, depth);
                        rhs(1) += (*S_data)(ci, depth) * x[0];
                        rhs(2) += (*S_data)(ci, depth) * x[1];
                        mat(0, 0) += 1;
                        mat(0, 1) += x[0];
                        mat(0, 2) += x[1];
                        mat(1, 0) += x[0];
                        mat(1, 1) += x[0] * x[0];
                        mat(1, 2) += x[1] * x[0];
                        mat(2, 0) += x[1];
                        mat(2, 1) += x[0] * x[1];
                        mat(2, 2) += x[1] * x[1];
                    }
                }
                ColPivHouseholderQR<Matrix3d> solver(mat);
                soln = solver.solve(rhs);
                q_val = soln(0) + soln(1) * X[0] + soln(2) * X[1];
                done = true;
            }
        }
    }
    //                 CellIndex<NDIM> ci(idx);
    //                 Pointer<CellData<NDIM, double> > u_c_data = patch->getPatchData(data_idx_temp);
    //                 const double* x_lower = p_geom->getXLower();
    //                 const Index<NDIM> ci_l = patch_box.lower();
    //                 std::vector<double> x_idx(NDIM);
    //                 for(int d = 0; d < NDIM; ++d)
    //                     x_idx[d] = x_lower[d] + dx[d]*(ci(d)-ci_l(d)+0.5);
    //                 std::vector<double> alpha(NDIM);
    //                 for(int d = 0; d < NDIM; ++d)
    //                     alpha[d] = (Q[d]-x_idx[d])/dx[d];
    //                 std::vector<double> q_y(4);
    //                 for (int j = -1; j < 3; ++j)
    //                 {
    //                     std::vector<double> q(4);
    //                     for (int i = -1; i < 3; ++i)
    //                         q[i+1] = (*u_c_data)(ci+IntVector<NDIM>(i,j));
    //                     weno4_(q.data(), alpha[0], q_y[j+1]);
    //                 }
    //                 weno4_(q_y.data(), alpha[1], q_val);
    //             }
    //         }
    //     }
    q_val = SAMRAI_MPI::sumReduction(q_val);

    for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
        level->deallocatePatchData(data_idx_temp);
    }
    return q_val;
}

} // namespace IBTK
