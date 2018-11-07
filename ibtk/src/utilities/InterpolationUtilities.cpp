#include "ibtk/InterpolationUtilities.h"
#include "ibtk/ibtk_utilities.h"
#include "ibtk/namespaces.h"
#include "tbox/Pointer.h"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/QR>
#include <HierarchyCellDataOpsReal.h>
#include <HierarchySideDataOpsReal.h>
#include <Patch.h>
#include <PatchLevel.h>
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
                double r;
                const CellIndex<NDIM> ci_l = patch_box.lower();
                std::vector<double> x(NDIM);
                // Moving least squares
                VectorXd rhs = VectorXd::Zero(6);
                VectorXd soln = VectorXd::Zero(6);
                MatrixXd mat = MatrixXd::Zero(6, 6);
                for (CellIterator<NDIM> i(box); i; i++)
                {
                    CellIndex<NDIM> ci = i();
                    for (int d = 0; d < NDIM; ++d) x[d] = x_lower[d] + dx[d] * (ci(d) - ci_l(d) + 0.5);
                    r = sqrt(x[0] * x[0] + x[1] * x[1]);
                    if (r > 1.0)
                    {
                        // Fill in RHS and MATRIX values
                        double w = weightFcn(X, x);
                        double f = (*S_data)(ci, depth);
                        rhs(0) += w * f;
                        rhs(1) += w * f * x[0];
                        rhs(2) += w * f * x[1];
                        rhs(3) += w * f * x[0] * x[0];
                        rhs(4) += w * f * x[1] * x[1];
                        rhs(5) += w * f * x[0] * x[1];
                        mat(0, 0) += w;
                        mat(1, 0) += x[0] * w;
                        mat(1, 1) += x[0] * x[0] * w;
                        mat(2, 0) += x[1] * w;
                        mat(2, 1) += x[1] * x[0] * w;
                        mat(2, 2) += x[1] * x[1] * w;
                        mat(3, 0) += x[0] * x[0] * w;
                        mat(3, 1) += x[0] * x[0] * x[0] * w;
                        mat(3, 2) += x[0] * x[0] * x[1] * w;
                        mat(3, 3) += x[0] * x[0] * x[0] * x[0] * w;
                        mat(4, 0) += x[1] * x[1] * w;
                        mat(4, 1) += x[1] * x[1] * x[0] * w;
                        mat(4, 2) += x[1] * x[1] * x[1] * w;
                        mat(4, 3) += x[1] * x[1] * x[0] * x[0] * w;
                        mat(4, 4) += x[1] * x[1] * x[1] * x[1] * w;
                        mat(5, 0) += x[0] * x[1] * w;
                        mat(5, 1) += x[0] * x[1] * x[0] * w;
                        mat(5, 2) += x[0] * x[1] * x[1] * w;
                        mat(5, 3) += x[0] * x[1] * x[0] * x[0] * w;
                        mat(5, 4) += x[0] * x[1] * x[1] * x[1] * w;
                        mat(5, 5) += x[0] * x[1] * x[0] * x[1] * w;
                    }
                }
                // mat(0,1) = mat(1,0); mat(0,2) = mat(2,0); mat(0,3) = mat(3,0); mat(0,4) = mat(4,0);
                // mat(0,5) = mat(5,0);
                // mat(1,2) = mat(2,1); mat(1,3) = mat(3,1); mat(1,4) = mat(4,1); mat(1,5) = mat(5,1);
                // mat(2,3) = mat(3,2); mat(2,4) = mat(4,2); mat(2,5) = mat(5,2);
                // mat(3,4) = mat(4,3); mat(3,5) = mat(5,3);
                // mat(4,5) = mat(5,4);
                soln = mat.ldlt().solve(rhs);
                q_val = soln(0) + soln(1) * X[0] + soln(2) * X[1] + soln(3) * X[0] * X[0] + soln(4) * X[1] * X[1] +
                        soln(5) * X[0] * X[1];
                // q_val = soln(1) + 2.0*soln(3)*X[0]+soln(5)*X[1];
                // pout << "relative error is :\n" << (mat*soln-rhs).norm() / rhs.norm() << "\n";
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

double
InterpolationUtilities::weightFcn(const std::vector<double>& x, const std::vector<double>& x_j)
{
    double r = sqrt((x[0] - x_j[0]) * (x[0] - x_j[0]) + (x[1] - x_j[1]) * (x[1] - x_j[1]));
    double w = exp(-10.0 * r * r);
    w = 1.0;
    return w;
}

void
InterpolationUtilities::gradient(const int grad_idx,
                                 const int data_idx,
                                 Pointer<SideVariable<NDIM, double> > Q_var,
                                 Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
                                 const std::vector<RobinBcCoefStrategy<NDIM>*>& bc_coefs,
                                 const double data_time,
                                 const int depth)
{
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    const int data_idx_temp =
        var_db->registerVariableAndContext(Q_var, var_db->getContext("Interpolation"), IntVector<NDIM>(3));
    for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
        level->allocatePatchData(data_idx_temp);
    }
    HierarchySideDataOpsReal<NDIM, double> hier_sc_data_ops(patch_hierarchy);
    hier_sc_data_ops.copyData(data_idx_temp, data_idx);
    typedef HierarchyGhostCellInterpolation::InterpolationTransactionComponent InterpolationTransactionComponent;
    std::vector<InterpolationTransactionComponent> ghost_cell_components(1);
    ghost_cell_components[0] = InterpolationTransactionComponent(
        data_idx_temp, "CONSERVATIVE_LINEAR_REFINE", true, "CONSERVATIVE_COARSEN", "LINEAR", false, bc_coefs, NULL);
    HierarchyGhostCellInterpolation ghost_fill_op;
    ghost_fill_op.initializeOperatorState(ghost_cell_components, patch_hierarchy);
    ghost_fill_op.fillData(data_time);
    for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Pointer<CartesianPatchGeometry<NDIM> > p_geom = patch->getPatchGeometry();
            const double* const dx = p_geom->getDx();
            const double* const x_lower = p_geom->getXLower();
            const Box<NDIM>& patch_box = patch->getBox();
            Pointer<SideData<NDIM, double> > S_data = patch->getPatchData(data_idx_temp);
            Pointer<CellData<NDIM, double> > grad_data = patch->getPatchData(grad_idx);
            const CellIndex<NDIM>& ci_l = patch_box.lower();
            for (CellIterator<NDIM> i(patch_box); i; i++)
            {
                CellIndex<NDIM> ci = i();
                Box<NDIM> box(ci, ci);
                box.grow(2);
                double rc = 0.0, r = 0.0;
                std::vector<double> xc(NDIM), x(NDIM);
                for (int d = 0; d < NDIM; ++d) xc[d] += x_lower[d] + dx[d] * (ci(d) - ci_l(d));
                rc = sqrt(xc[0] * xc[0] + xc[1] * xc[1]);

                VectorXd rhs = VectorXd::Zero(6);
                VectorXd soln = VectorXd::Zero(6);
                MatrixXd mat = MatrixXd::Zero(6, 6);
                for (int axis = 0; axis < NDIM; ++axis)
                {
                    rhs = VectorXd::Zero(6);
                    soln = VectorXd::Zero(6);
                    mat = MatrixXd::Zero(6, 6);
                    for (SideIterator<NDIM> i(box, axis); i; i++)
                    {
                        SideIndex<NDIM> si = i();
                        for (int d = 0; d < NDIM; ++d)
                            x[d] = x_lower[d] + dx[d] * (si(d) - ci_l(d) + (d == axis ? 0.0 : 0.5));
                        r = sqrt(x[0] * x[0] + x[1] * x[1]);
                        if (((r > 1.0) && (rc > 1.0)) || ((r < 1.0) && (rc < 1.0)))
                        {
                            double w = weightFcn(xc, x);
                            double f = (*S_data)(si, depth);
                            rhs(0) += w * f;
                            rhs(1) += w * f * x[0];
                            rhs(2) += w * f * x[1];
                            rhs(3) += w * f * x[0] * x[0];
                            rhs(4) += w * f * x[1] * x[1];
                            rhs(5) += w * f * x[0] * x[1];
                            mat(0, 0) += w;
                            mat(1, 0) += x[0] * w;
                            mat(1, 1) += x[0] * x[0] * w;
                            mat(2, 0) += x[1] * w;
                            mat(2, 1) += x[1] * x[0] * w;
                            mat(2, 2) += x[1] * x[1] * w;
                            mat(3, 0) += x[0] * x[0] * w;
                            mat(3, 1) += x[0] * x[0] * x[0] * w;
                            mat(3, 2) += x[0] * x[0] * x[1] * w;
                            mat(3, 3) += x[0] * x[0] * x[0] * x[0] * w;
                            mat(4, 0) += x[1] * x[1] * w;
                            mat(4, 1) += x[1] * x[1] * x[0] * w;
                            mat(4, 2) += x[1] * x[1] * x[1] * w;
                            mat(4, 3) += x[1] * x[1] * x[0] * x[0] * w;
                            mat(4, 4) += x[1] * x[1] * x[1] * x[1] * w;
                            mat(5, 0) += x[0] * x[1] * w;
                            mat(5, 1) += x[0] * x[1] * x[0] * w;
                            mat(5, 2) += x[0] * x[1] * x[1] * w;
                            mat(5, 3) += x[0] * x[1] * x[0] * x[0] * w;
                            mat(5, 4) += x[0] * x[1] * x[1] * x[1] * w;
                            mat(5, 5) += x[0] * x[1] * x[0] * x[1] * w;
                        }
                    }
                    soln = mat.ldlt().solve(rhs);
                    (*grad_data)(ci, 2 * axis) = soln(1) + 2.0 * soln(3) * xc[0] + soln(5) * xc[1];
                    (*grad_data)(ci, 2 * axis + 1) = soln(2) + 2.0 * soln(4) * xc[1] + soln(5) * xc[0];
                }
            }
        }
    }
    for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
        level->deallocatePatchData(data_idx_temp);
    }

    return;
}

void
InterpolationUtilities::gradient(std::vector<double>& ret,
                                 const std::vector<double>& X,
                                 const int data_idx,
                                 Pointer<SideVariable<NDIM, double> > Q_var,
                                 SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > patch_hierarchy,
                                 const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& bc_coefs,
                                 const double data_time,
                                 const int depth)
{
    ret.resize(NDIM * NDIM);
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    const int data_idx_temp =
        var_db->registerVariableAndContext(Q_var, var_db->getContext("Interpolation"), IntVector<NDIM>(3));
    for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
        level->allocatePatchData(data_idx_temp);
    }
    HierarchySideDataOpsReal<NDIM, double> hier_sc_data_ops(patch_hierarchy);
    hier_sc_data_ops.copyData(data_idx_temp, data_idx);
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
            Pointer<SideData<NDIM, double> > S_data = patch->getPatchData(data_idx_temp);
            if (patch_box.contains(idx))
            {
                // Great. The patch is currently on this level
                // Let's create a box that contains this data
                Box<NDIM> box(idx, idx);
                // Grow it by some number of grid cells
                box.grow(2);
                double r;
                const CellIndex<NDIM> ci_l = patch_box.lower();
                std::vector<double> x(NDIM);
                // Moving least squares
                VectorXd rhs_x = VectorXd::Zero(6), rhs_y = VectorXd::Zero(6);
                VectorXd soln = VectorXd::Zero(6);
                MatrixXd mat_x = MatrixXd::Zero(6, 6), mat_y = MatrixXd::Zero(6, 6);
                for (int axis = 0; axis < NDIM; ++axis)
                {
                    rhs_x = VectorXd::Zero(6);
                    soln = VectorXd::Zero(6);
                    mat_x = MatrixXd::Zero(6, 6);
                    for (SideIterator<NDIM> i(box, axis); i; i++)
                    {
                        SideIndex<NDIM> si = i();
                        for (int d = 0; d < NDIM; ++d)
                            x[d] = x_lower[d] + dx[d] * (si(d) - ci_l(d) + (d == axis ? 0.0 : 0.5));
                        r = sqrt(x[0] * x[0] + x[1] * x[1]);
                        if (r > 1.0)
                        {
                            // Fill in RHS and MATRIX values
                            double w = weightFcn(X, x);
                            double f = (*S_data)(si, depth);
                            rhs_x(0) += w * f;
                            rhs_x(1) += w * f * x[0];
                            rhs_x(2) += w * f * x[1];
                            rhs_x(3) += w * f * x[0] * x[0];
                            rhs_x(4) += w * f * x[1] * x[1];
                            rhs_x(5) += w * f * x[0] * x[1];
                            mat_x(0, 0) += w;
                            mat_x(1, 0) += x[0] * w;
                            mat_x(1, 1) += x[0] * x[0] * w;
                            mat_x(2, 0) += x[1] * w;
                            mat_x(2, 1) += x[1] * x[0] * w;
                            mat_x(2, 2) += x[1] * x[1] * w;
                            mat_x(3, 0) += x[0] * x[0] * w;
                            mat_x(3, 1) += x[0] * x[0] * x[0] * w;
                            mat_x(3, 2) += x[0] * x[0] * x[1] * w;
                            mat_x(3, 3) += x[0] * x[0] * x[0] * x[0] * w;
                            mat_x(4, 0) += x[1] * x[1] * w;
                            mat_x(4, 1) += x[1] * x[1] * x[0] * w;
                            mat_x(4, 2) += x[1] * x[1] * x[1] * w;
                            mat_x(4, 3) += x[1] * x[1] * x[0] * x[0] * w;
                            mat_x(4, 4) += x[1] * x[1] * x[1] * x[1] * w;
                            mat_x(5, 0) += x[0] * x[1] * w;
                            mat_x(5, 1) += x[0] * x[1] * x[0] * w;
                            mat_x(5, 2) += x[0] * x[1] * x[1] * w;
                            mat_x(5, 3) += x[0] * x[1] * x[0] * x[0] * w;
                            mat_x(5, 4) += x[0] * x[1] * x[1] * x[1] * w;
                            mat_x(5, 5) += x[0] * x[1] * x[0] * x[1] * w;
                        }
                    }
                    mat_x(0, 1) = mat_x(1, 0);
                    mat_x(0, 2) = mat_x(2, 0);
                    mat_x(0, 3) = mat_x(3, 0);
                    mat_x(0, 4) = mat_x(4, 0);
                    mat_x(0, 5) = mat_x(5, 0);
                    mat_x(1, 2) = mat_x(2, 1);
                    mat_x(1, 3) = mat_x(3, 1);
                    mat_x(1, 4) = mat_x(4, 1);
                    mat_x(1, 5) = mat_x(5, 1);
                    mat_x(2, 3) = mat_x(3, 2);
                    mat_x(2, 4) = mat_x(4, 2);
                    mat_x(2, 5) = mat_x(5, 2);
                    mat_x(3, 4) = mat_x(4, 3);
                    mat_x(3, 5) = mat_x(5, 3);
                    mat_x(4, 5) = mat_x(5, 4);
                    soln = mat_x.ldlt().solve(rhs_x);
                    ret[2 * axis] = soln(1) + 2.0 * soln(3) * X[0] + soln(5) * X[1];
                    ret[2 * axis + 1] = soln(2) + 2.0 * soln(4) * X[1] + soln(5) * X[0];
                    pout << "relative error is :\n" << (mat_x * soln - rhs_x).norm() / rhs_x.norm() << "\n";
                    pout << rhs_x << "\n";
                    pout << mat_x << "\n";
                    pout << soln << "\n";
                    done = true;
                }
            }
        }
    }
    SAMRAI_MPI::sumReduction(ret.data(), ret.size());

    for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
        level->deallocatePatchData(data_idx_temp);
    }
    return;
}

} // namespace IBTK
