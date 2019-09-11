#include "ibamr/AdvDiffComplexFluidConvectiveOperator.h"
#include <HierarchyDataOpsManager.h>
#include <fstream>
#include <iostream>

extern "C"
{
#if (NDIM == 2)
    void conform_tens_conv_u_s_oper_2d_(const double*,
                                        const double*,
                                        const double*,
                                        const int&,
                                        const double*,
                                        const int&,
                                        const double*,
                                        const int&,
                                        const double*,
                                        const int&,
                                        const double*,
                                        const int&,
                                        const int&,
                                        const int&,
                                        const int&,
                                        const int&);
    void sqrt_tens_conv_u_s_oper_2d_(const double*,
                                     const double*,
                                     const double*,
                                     const int&,
                                     const double*,
                                     const int&,
                                     const double*,
                                     const int&,
                                     const double*,
                                     const int&,
                                     const double*,
                                     const int&,
                                     const int&,
                                     const int&,
                                     const int&,
                                     const int&);
    void log_tens_conv_u_s_oper_2d_(const double*,
                                    const double*,
                                    const double*,
                                    const int&,
                                    const double*,
                                    const int&,
                                    const double*,
                                    const int&,
                                    const double*,
                                    const int&,
                                    const double*,
                                    const int&,
                                    const int&,
                                    const int&,
                                    const int&,
                                    const int&);
    void conform_tens_conv_u_c_oper_2d_(const double*,
                                        const double*,
                                        const int&,
                                        const double*,
                                        const int&,
                                        const double*,
                                        const int&,
                                        const double*,
                                        const int&,
                                        const double*,
                                        const int&,
                                        const int&,
                                        const int&,
                                        const int&,
                                        const int&);
    void sqrt_tens_conv_u_c_oper_2d_(const double*,
                                     const double*,
                                     const int&,
                                     const double*,
                                     const int&,
                                     const double*,
                                     const int&,
                                     const double*,
                                     const int&,
                                     const double*,
                                     const int&,
                                     const int&,
                                     const int&,
                                     const int&,
                                     const int&);
    void log_tens_conv_u_c_oper_2d_(const double*,
                                    const double*,
                                    const int&,
                                    const double*,
                                    const int&,
                                    const double*,
                                    const int&,
                                    const double*,
                                    const int&,
                                    const double*,
                                    const int&,
                                    const int&,
                                    const int&,
                                    const int&,
                                    const int&);
#endif
#if (NDIM == 3)
    void conform_tens_conv_u_s_oper_3d_(const double*,
                                        const double*,
                                        const double*,
                                        const double*,
                                        const int&,
                                        const double*,
                                        const int&,
                                        const double*,
                                        const int&,
                                        const double*,
                                        const int&,
                                        const double*,
                                        const int&,
                                        const int&,
                                        const int&,
                                        const int&,
                                        const int&,
                                        const int&,
                                        const int&);
    void sqrt_tens_conv_u_s_oper_3d_(const double*,
                                     const double*,
                                     const double*,
                                     const double*,
                                     const int&,
                                     const double*,
                                     const int&,
                                     const double*,
                                     const int&,
                                     const double*,
                                     const int&,
                                     const double*,
                                     const int&,
                                     const int&,
                                     const int&,
                                     const int&,
                                     const int&,
                                     const int&,
                                     const int&);
    void log_tens_conv_u_s_oper_3d_(const double*,
                                    const double*,
                                    const double*,
                                    const double*,
                                    const int&,
                                    const double*,
                                    const int&,
                                    const double*,
                                    const int&,
                                    const double*,
                                    const int&,
                                    const double*,
                                    const int&,
                                    const int&,
                                    const int&,
                                    const int&,
                                    const int&,
                                    const int&,
                                    const int&);
    void conform_tens_conv_u_c_oper_3d_(const double*,
                                        const double*,
                                        const int&,
                                        const double*,
                                        const int&,
                                        const double*,
                                        const int&,
                                        const double*,
                                        const int&,
                                        const double*,
                                        const int&,
                                        const int&,
                                        const int&,
                                        const int&,
                                        const int&,
                                        const int&,
                                        const int&);
    void sqrt_tens_conv_u_c_oper_3d_(const double*,
                                     const double*,
                                     const int&,
                                     const double*,
                                     const int&,
                                     const double*,
                                     const int&,
                                     const double*,
                                     const int&,
                                     const double*,
                                     const int&,
                                     const int&,
                                     const int&,
                                     const int&,
                                     const int&,
                                     const int&,
                                     const int&);
    void log_tens_conv_u_c_oper_3d_(const double*,
                                    const double*,
                                    const int&,
                                    const double*,
                                    const int&,
                                    const double*,
                                    const int&,
                                    const double*,
                                    const int&,
                                    const double*,
                                    const int&,
                                    const int&,
                                    const int&,
                                    const int&,
                                    const int&,
                                    const int&,
                                    const int&);
#endif
}

namespace IBAMR
{
static const IntVector<NDIM> ghosts_cc = 1;

// Constructor
AdvDiffComplexFluidConvectiveOperator::AdvDiffComplexFluidConvectiveOperator(
    const std::string& object_name,
    Pointer<CellVariable<NDIM, double> > Q_var,
    Pointer<Database> input_db,
    const CFConvectiveOperatorType difference_form,
    const std::vector<RobinBcCoefStrategy<NDIM>*>& conc_bc_coefs,
    const std::vector<RobinBcCoefStrategy<NDIM>*>& u_bc_coefs)
    : ConvectiveOperator(object_name, UNKNOWN_CONVECTIVE_DIFFERENCING_TYPE),
      d_hierarchy(NULL),
      d_coarsest_ln(-1),
      d_finest_ln(-1),
      d_Q_var(Q_var),
      d_u_adv_var(new SideVariable<NDIM, double>("Complex U var")),
      d_u_scratch_idx(-1),
      d_difference_form(difference_form),
      d_convec_oper(NULL),
      d_Q_convec_idx(-1),
      d_conc_bc_coefs(conc_bc_coefs),
      d_u_bc_coefs(u_bc_coefs),
      d_rhs(NULL),
      d_R_idx(-1),
      d_conform_tens(true),
      d_sqr_root(false),
      d_log_conform(false),
      d_lambda(-1),
      d_eta(-1)
{
    d_sqr_root = input_db->getBoolWithDefault("square_root_evolve", false);
    d_log_conform = input_db->getBoolWithDefault("log_conform_evolve", false);
    // Register some scratch variables
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<VariableContext> context = var_db->getContext(d_object_name + "::CONTEXT");
    d_Q_convec_idx = var_db->registerVariableAndContext(
        d_Q_var, var_db->getContext(d_object_name + "::CONVECTIVE"), IntVector<NDIM>(0));
    d_lambda = input_db->getDouble("relaxation_time");
    Pointer<VariableContext> rhs_ctx = var_db->getContext(d_object_name + "::RHS_CONTEXT");
    d_R_idx = var_db->registerVariableAndContext(d_Q_var, rhs_ctx);
    Pointer<VariableContext> new_cxt = var_db->getContext(d_object_name + "::U_ADV_CXT");
    d_u_scratch_idx = var_db->registerVariableAndContext(d_u_adv_var, new_cxt, /*ghosts_cc*/ IntVector<NDIM>(7));
    switch (d_difference_form)
    {
    case CENTERED:
        d_convec_oper = new AdvDiffCenteredConvectiveOperator(
            d_object_name + "::convec_oper", d_Q_var, input_db, ADVECTIVE, conc_bc_coefs);
        break;
    case WAVE_PROP:
        d_convec_oper = new AdvDiffWavePropConvectiveOperator(
            d_object_name + "::convec_oper", d_Q_var, input_db, ADVECTIVE, conc_bc_coefs);
        break;
    case CF_PPM:
        d_convec_oper = new AdvDiffPPMConvectiveOperator(
            d_object_name + "::convec_oper", d_Q_var, input_db, ADVECTIVE, conc_bc_coefs);
        break;
    case UNKNOWN:
        TBOX_ERROR("\n\n CONVECTIVE OPERATOR: " + enum_to_string<CFConvectiveOperatorType>(d_difference_form) +
                   " UNKNOWN\n\n");
        break;
    }
} // Constructor
AdvDiffComplexFluidConvectiveOperator::~AdvDiffComplexFluidConvectiveOperator()
{
    deallocateOperatorState();
    // intentionally blank
    return;
}
void
AdvDiffComplexFluidConvectiveOperator::applyConvectiveOperator(int Q_idx, int Y_idx)
{
    if (!d_is_initialized)
    {
        TBOX_ERROR("AdvDiffCenteredConvectiveOperator::applyConvectiveOperator():\n"
                   << "  operator must be initialized prior to call to applyConvectiveOperator\n");
    }
    // Set up velocity information:
    d_convec_oper->setAdvectionVelocity(d_u_idx);
    d_convec_oper->setSolutionTime(d_solution_time);
    Pointer<CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();
    // Copy data to side centered velocity field
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            Pointer<FaceData<NDIM, double> > u_f_data = patch->getPatchData(d_u_idx);
            Pointer<SideData<NDIM, double> > u_s_data = patch->getPatchData(d_u_scratch_idx);
            const Box<NDIM> box = patch->getBox();
            for (int axis = 0; axis < NDIM; ++axis)
            {
                for (SideIterator<NDIM> i(box, axis); i; i++)
                {
                    SideIndex<NDIM> si = i();
                    FaceIndex<NDIM> fi(si.toCell(0), axis, 1);
                    (*u_s_data)(si) = (*u_f_data)(fi);
                }
            }
        }
    }
    // Fill boundary conditions for side centered velocity
    typedef HierarchyGhostCellInterpolation::InterpolationTransactionComponent InterpolationTransactionComponent;
    std::vector<InterpolationTransactionComponent> ghost_cell_components(1);
    ghost_cell_components[0] = InterpolationTransactionComponent(
        d_u_scratch_idx, "NONE", true, "CUBIC_COARSEN", "LINEAR", false, d_u_bc_coefs, NULL);
    HierarchyGhostCellInterpolation ghost_fill_op;
    ghost_fill_op.initializeOperatorState(ghost_cell_components, d_hierarchy);
    ghost_fill_op.fillData(d_solution_time);

    d_rhs->setPatchDataIndex(Q_idx);
    d_rhs->setDataOnPatchHierarchy(d_R_idx, d_Q_var, d_hierarchy, d_solution_time);

    d_convec_oper->applyConvectiveOperator(Q_idx, d_Q_convec_idx);

    for (int level_num = d_coarsest_ln; level_num <= d_finest_ln; ++level_num)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(level_num);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            Pointer<SideData<NDIM, double> > u_ADV_s_data = patch->getPatchData(d_u_scratch_idx);
            Pointer<CellData<NDIM, double> > u_ADV_c_data = patch->getPatchData(d_u_scratch_idx);
            if (u_ADV_s_data)
            {
                const Pointer<CartesianPatchGeometry<NDIM> > p_geom = patch->getPatchGeometry();
                const double* dx = p_geom->getDx();
                const Box<NDIM>& patch_box = patch->getBox();
                const IntVector<NDIM> patch_lower = patch_box.lower();
                const IntVector<NDIM> patch_upper = patch_box.upper();
                Pointer<CellData<NDIM, double> > Q_data = patch->getPatchData(Q_idx);
                Pointer<CellData<NDIM, double> > Y_data = patch->getPatchData(Y_idx);
                const IntVector<NDIM> Q_data_gcw = Q_data->getGhostCellWidth();

                const IntVector<NDIM> u_data_gcw = u_ADV_s_data->getGhostCellWidth();
                const IntVector<NDIM> Y_data_gcw = Y_data->getGhostCellWidth();
                Pointer<CellData<NDIM, double> > C_data = patch->getPatchData(d_Q_convec_idx);
                const IntVector<NDIM> C_data_gcw = C_data->getGhostCellWidth();
                Pointer<CellData<NDIM, double> > R_data = patch->getPatchData(d_R_idx);
                const IntVector<NDIM> R_data_gcw = R_data->getGhostCellWidth();
                if (d_sqr_root)
                {
#if (NDIM == 2)
                    sqrt_tens_conv_u_s_oper_2d_(dx,
                                                u_ADV_s_data->getPointer(0),
                                                u_ADV_s_data->getPointer(1),
                                                u_data_gcw.max(),
                                                Q_data->getPointer(0),
                                                Q_data_gcw.max(),
                                                R_data->getPointer(0),
                                                R_data_gcw.max(),
                                                C_data->getPointer(0),
                                                C_data_gcw.max(),
                                                Y_data->getPointer(0),
                                                Y_data_gcw.max(),
                                                patch_lower(0),
                                                patch_upper(0),
                                                patch_lower(1),
                                                patch_upper(1));
#endif
#if (NDIM == 3)
                    sqrt_tens_conv_u_s_oper_3d_(dx,
                                                u_ADV_s_data->getPointer(0),
                                                u_ADV_s_data->getPointer(1),
                                                u_ADV_s_data->getPointer(2),
                                                u_data_gcw.max(),
                                                Q_data->getPointer(0),
                                                Q_data_gcw.max(),
                                                R_data->getPointer(0),
                                                R_data_gcw.max(),
                                                C_data->getPointer(0),
                                                C_data_gcw.max(),
                                                Y_data->getPointer(0),
                                                Y_data_gcw.max(),
                                                patch_lower(0),
                                                patch_upper(0),
                                                patch_lower(1),
                                                patch_upper(1),
                                                patch_lower(2),
                                                patch_upper(2));
#endif
                }
                else if (d_log_conform)
                {
#if (NDIM == 2)
                    log_tens_conv_u_s_oper_2d_(dx,
                                               u_ADV_s_data->getPointer(0),
                                               u_ADV_s_data->getPointer(1),
                                               u_data_gcw.max(),
                                               Q_data->getPointer(0),
                                               Q_data_gcw.max(),
                                               R_data->getPointer(0),
                                               R_data_gcw.max(),
                                               C_data->getPointer(0),
                                               C_data_gcw.max(),
                                               Y_data->getPointer(0),
                                               Y_data_gcw.max(),
                                               patch_lower(0),
                                               patch_upper(0),
                                               patch_lower(1),
                                               patch_upper(1));
#endif
#if (NDIM == 3)
                    log_tens_conv_u_s_oper_3d_(dx,
                                               u_ADV_s_data->getPointer(0),
                                               u_ADV_s_data->getPointer(1),
                                               u_ADV_s_data->getPointer(2),
                                               u_data_gcw.max(),
                                               Q_data->getPointer(0),
                                               Q_data_gcw.max(),
                                               R_data->getPointer(0),
                                               R_data_gcw.max(),
                                               C_data->getPointer(0),
                                               C_data_gcw.max(),
                                               Y_data->getPointer(0),
                                               Y_data_gcw.max(),
                                               patch_lower(0),
                                               patch_upper(0),
                                               patch_lower(1),
                                               patch_upper(1),
                                               patch_lower(2),
                                               patch_upper(2));
#endif
                }
                else
                {
#if (NDIM == 2)
                    conform_tens_conv_u_s_oper_2d_(dx,
                                                   u_ADV_s_data->getPointer(0),
                                                   u_ADV_s_data->getPointer(1),
                                                   u_data_gcw.max(),
                                                   Q_data->getPointer(0),
                                                   Q_data_gcw.max(),
                                                   R_data->getPointer(0),
                                                   R_data_gcw.max(),
                                                   C_data->getPointer(0),
                                                   C_data_gcw.max(),
                                                   Y_data->getPointer(0),
                                                   Y_data_gcw.max(),
                                                   patch_lower(0),
                                                   patch_upper(0),
                                                   patch_lower(1),
                                                   patch_upper(1));
#endif
#if (NDIM == 3)
                    conform_tens_conv_u_s_oper_3d_(dx,
                                                   u_ADV_s_data->getPointer(0),
                                                   u_ADV_s_data->getPointer(1),
                                                   u_ADV_s_data->getPointer(2),
                                                   u_data_gcw.max(),
                                                   Q_data->getPointer(0),
                                                   Q_data_gcw.max(),
                                                   R_data->getPointer(0),
                                                   R_data_gcw.max(),
                                                   C_data->getPointer(0),
                                                   C_data_gcw.max(),
                                                   Y_data->getPointer(0),
                                                   Y_data_gcw.max(),
                                                   patch_lower(0),
                                                   patch_upper(0),
                                                   patch_lower(1),
                                                   patch_upper(1),
                                                   patch_lower(2),
                                                   patch_upper(2));
#endif
                }
            }
            else if (u_ADV_c_data)
            {
                const Pointer<CartesianPatchGeometry<NDIM> > p_geom = patch->getPatchGeometry();
                const double* dx = p_geom->getDx();
                const Box<NDIM>& patch_box = patch->getBox();
                const IntVector<NDIM> patch_lower = patch_box.lower();
                const IntVector<NDIM> patch_upper = patch_box.upper();
                Pointer<CellData<NDIM, double> > Q_data = patch->getPatchData(Q_idx);
                Pointer<CellData<NDIM, double> > Y_data = patch->getPatchData(Y_idx);
                const IntVector<NDIM> Q_data_gcw = Q_data->getGhostCellWidth();
                const IntVector<NDIM> u_data_gcw = u_ADV_c_data->getGhostCellWidth();
                const IntVector<NDIM> Y_data_gcw = Y_data->getGhostCellWidth();
                Pointer<CellData<NDIM, double> > C_data = patch->getPatchData(d_Q_convec_idx);
                const IntVector<NDIM> C_data_gcw = C_data->getGhostCellWidth();
                Pointer<CellData<NDIM, double> > R_data = patch->getPatchData(d_R_idx);
                const IntVector<NDIM> R_data_gcw = R_data->getGhostCellWidth();
                if (d_sqr_root)
                {
#if (NDIM == 2)
                    sqrt_tens_conv_u_c_oper_2d_(dx,
                                                u_ADV_c_data->getPointer(0),
                                                u_data_gcw.max(),
                                                Q_data->getPointer(0),
                                                Q_data_gcw.max(),
                                                R_data->getPointer(0),
                                                R_data_gcw.max(),
                                                C_data->getPointer(0),
                                                C_data_gcw.max(),
                                                Y_data->getPointer(0),
                                                Y_data_gcw.max(),
                                                patch_lower(0),
                                                patch_upper(0),
                                                patch_lower(1),
                                                patch_upper(1));
#endif
#if (NDIM == 3)
                    sqrt_tens_conv_u_c_oper_3d_(dx,
                                                u_ADV_c_data->getPointer(0),
                                                u_data_gcw.max(),
                                                Q_data->getPointer(0),
                                                Q_data_gcw.max(),
                                                R_data->getPointer(0),
                                                R_data_gcw.max(),
                                                C_data->getPointer(0),
                                                C_data_gcw.max(),
                                                Y_data->getPointer(0),
                                                Y_data_gcw.max(),
                                                patch_lower(0),
                                                patch_upper(0),
                                                patch_lower(1),
                                                patch_upper(1),
                                                patch_lower(2),
                                                patch_upper(2));
#endif
                }
                else if (d_log_conform)
                {
#if (NDIM == 2)
                    log_tens_conv_u_c_oper_2d_(dx,
                                               u_ADV_c_data->getPointer(0),
                                               u_data_gcw.max(),
                                               Q_data->getPointer(0),
                                               Q_data_gcw.max(),
                                               R_data->getPointer(0),
                                               R_data_gcw.max(),
                                               C_data->getPointer(0),
                                               C_data_gcw.max(),
                                               Y_data->getPointer(0),
                                               Y_data_gcw.max(),
                                               patch_lower(0),
                                               patch_upper(0),
                                               patch_lower(1),
                                               patch_upper(1));
#endif
#if (NDIM == 3)
                    log_tens_conv_u_c_oper_3d_(dx,
                                               u_ADV_c_data->getPointer(0),
                                               u_data_gcw.max(),
                                               Q_data->getPointer(0),
                                               Q_data_gcw.max(),
                                               R_data->getPointer(0),
                                               R_data_gcw.max(),
                                               C_data->getPointer(0),
                                               C_data_gcw.max(),
                                               Y_data->getPointer(0),
                                               Y_data_gcw.max(),
                                               patch_lower(0),
                                               patch_upper(0),
                                               patch_lower(1),
                                               patch_upper(1),
                                               patch_lower(2),
                                               patch_upper(2));
#endif
                }
                else
                {
#if (NDIM == 2)
                    conform_tens_conv_u_c_oper_2d_(dx,
                                                   u_ADV_c_data->getPointer(0),
                                                   u_data_gcw.max(),
                                                   Q_data->getPointer(0),
                                                   Q_data_gcw.max(),
                                                   R_data->getPointer(0),
                                                   R_data_gcw.max(),
                                                   C_data->getPointer(0),
                                                   C_data_gcw.max(),
                                                   Y_data->getPointer(0),
                                                   Y_data_gcw.max(),
                                                   patch_lower(0),
                                                   patch_upper(0),
                                                   patch_lower(1),
                                                   patch_upper(1));
#endif
#if (NDIM == 3)
                    conform_tens_conv_u_c_oper_3d_(dx,
                                                   u_ADV_c_data->getPointer(0),
                                                   u_data_gcw.max(),
                                                   Q_data->getPointer(0),
                                                   Q_data_gcw.max(),
                                                   R_data->getPointer(0),
                                                   R_data_gcw.max(),
                                                   C_data->getPointer(0),
                                                   C_data_gcw.max(),
                                                   Y_data->getPointer(0),
                                                   Y_data_gcw.max(),
                                                   patch_lower(0),
                                                   patch_upper(0),
                                                   patch_lower(1),
                                                   patch_upper(1),
                                                   patch_lower(2),
                                                   patch_upper(2));
#endif
                }
            }
        } // end Patch loop
    }     // end Level loop
} // end applyConvectiveOperator
void
AdvDiffComplexFluidConvectiveOperator::initializeOperatorState(const SAMRAIVectorReal<NDIM, double>& in,
                                                               const SAMRAIVectorReal<NDIM, double>& out)
{
    if (d_is_initialized) deallocateOperatorState();
    // Get Hierarchy Information
    d_hierarchy = in.getPatchHierarchy();
    d_coarsest_ln = in.getCoarsestLevelNumber();
    d_finest_ln = in.getFinestLevelNumber();
#if !defined(NDEBUG)
    TBOX_ASSERT(d_hierarchy == out.getPatchHierarchy());
    TBOX_ASSERT(d_coarsest_ln == out.getCoarsestLevelNumber());
    TBOX_ASSERT(d_finest_ln == out.getFinestLevelNumber());
#else
    NULL_USE(out);
#endif
    d_convec_oper->initializeOperatorState(in, out);
    d_convec_oper->setAdvectionVelocity(d_u_idx);
    // Allocate Patch Data
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        if (!level->checkAllocated(d_u_scratch_idx)) level->allocatePatchData(d_u_scratch_idx);

        if (!level->checkAllocated(d_Q_convec_idx)) level->allocatePatchData(d_Q_convec_idx);

        if (!level->checkAllocated(d_R_idx)) level->allocatePatchData(d_R_idx);
    }
    d_is_initialized = true;
    return;
}
void
AdvDiffComplexFluidConvectiveOperator::deallocateOperatorState()
{
    if (!d_is_initialized) return;
    d_convec_oper->deallocateOperatorState();
    // Deallocate scratch data
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        if (level->checkAllocated(d_Q_convec_idx)) level->deallocatePatchData(d_Q_convec_idx);
        if (level->checkAllocated(d_u_scratch_idx)) level->deallocatePatchData(d_u_scratch_idx);
        if (level->checkAllocated(d_R_idx)) level->deallocatePatchData(d_R_idx);
    }
    d_is_initialized = false;
    return;
}

void
AdvDiffComplexFluidConvectiveOperator::registerRHSTerm(Pointer<RHS_Operator> rhs)
{
    d_rhs = rhs;
    return;
}
} // namespace IBAMR
