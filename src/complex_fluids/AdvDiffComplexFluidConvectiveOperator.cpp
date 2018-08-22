#include "ibamr/AdvDiffComplexFluidConvectiveOperator.h"
#include <HierarchyDataOpsManager.h>
#include <boost/lexical_cast.hpp>
#include <fstream>
#include <iostream>

extern "C"
{
#if (NDIM == 2)
    void conform_tens_conv_u_f_oper_2d_(const double*,
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
                                        const double&);
    void sqrt_tens_conv_u_f_oper_2d_(const double*,
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
                                     const double&);
    void log_tens_conv_u_f_oper_2d_(const double*,
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
                                    const double&);
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
                                        const int&,
                                        const double&);
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
                                     const int&,
                                     const double&);
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
                                    const int&,
                                    const double&);
#endif
#if (NDIM == 3)
    void conform_tens_conv_u_f_oper_3d_(const double*,
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
                                        const int&,
                                        const double&);
    void sqrt_tens_conv_u_f_oper_3d_(const double*,
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
                                     const int&,
                                     const double&);
    void log_tens_conv_u_f_oper_3d_(const double*,
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
                                    const int&,
                                    const double&);
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
                                        const int&,
                                        const double&);
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
                                     const int&,
                                     const double&);
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
                                    const int&,
                                    const double&);
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
    Pointer<FaceVariable<NDIM, double> > u_var)
    : ConvectiveOperator(object_name, UNKNOWN_CONVECTIVE_DIFFERENCING_TYPE),
      d_ghostfill_alg_Q(NULL),
      d_ghostfill_alg_u(NULL),
      d_ghostfill_scheds_u(),
      d_ghostfill_scheds_Q(),
      d_conc_bc_coefs(conc_bc_coefs),
      d_outflow_bdry_extrap_type("CONSTANT"),
      d_hierarchy(NULL),
      d_coarsest_ln(-1),
      d_finest_ln(-1),
      d_Q_var(Q_var),
      d_Q_data_depth(0),
      d_rhs(NULL),
      d_R_idx(-1),
      d_lambda(-1),
      d_eta(-1),
      d_sqr_root(false),
      d_log_conform(false),
      d_difference_form(difference_form),
      d_u_adv_var(u_var)
{
    d_sqr_root = input_db->getBoolWithDefault("square_root_evolve", false);
    d_log_conform = input_db->getBoolWithDefault("log_conform_evolve", false);
    // Register some scratch variables
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<VariableContext> context = var_db->getContext(d_object_name + "::CONTEXT");
    d_Q_scratch_idx = var_db->registerVariableAndContext(d_Q_var, context, ghosts_cc);
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
    Pointer<CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();
    // Set up refine algorithms for Q and u.
    Pointer<RefineAlgorithm<NDIM> > refine_alg_u = new RefineAlgorithm<NDIM>();
    Pointer<RefineOperator<NDIM> > refine_op_u =
        grid_geom->lookupRefineOperator(d_u_adv_var, "CONSERVATIVE_LINEAR_REFINE");
    refine_alg_u->registerRefine(d_u_scratch_idx, d_u_idx, d_u_scratch_idx, refine_op_u);
    d_ghostfill_alg_u = new RefineAlgorithm<NDIM>();
    d_ghostfill_alg_u->registerRefine(d_u_scratch_idx, d_u_idx, d_u_scratch_idx, refine_op_u);
    d_ghostfill_strategy_u = new CartExtrapPhysBdryOp(d_u_scratch_idx, "LINEAR");
    d_ghostfill_scheds_u.resize(d_finest_ln + 1);
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        d_ghostfill_scheds_u[ln] =
            d_ghostfill_alg_u->createSchedule(level, ln - 1, d_hierarchy, d_ghostfill_strategy_u);
    }
    // Set up coarsen algorithms for Q and u.
    Pointer<CoarsenAlgorithm<NDIM> > coarsen_alg_u = new CoarsenAlgorithm<NDIM>();
    Pointer<CoarsenOperator<NDIM> > coarsen_op_u =
        grid_geom->lookupCoarsenOperator(d_u_adv_var, "CONSERVATIVE_COARSEN");
    coarsen_alg_u->registerCoarsen(d_u_scratch_idx, d_u_scratch_idx, coarsen_op_u);
    // Refine the data for Q and u
    d_ghostfill_scheds_u.resize(d_finest_ln + 1);
    for (int level_num = d_coarsest_ln; level_num <= d_finest_ln; ++level_num)
    {
        refine_alg_u->resetSchedule(d_ghostfill_scheds_u[level_num]);
        d_ghostfill_scheds_u[level_num]->fillData(d_solution_time);
        d_ghostfill_alg_u->resetSchedule(d_ghostfill_scheds_u[level_num]);
    }
    for (int level_num = d_finest_ln; level_num > d_coarsest_ln; --level_num)
    {
        d_coarsen_scheds_u[level_num]->coarsenData();
    }

    d_rhs->setPatchDataIndex(Q_idx);
    d_rhs->setDataOnPatchHierarchy(d_R_idx, d_Q_var, d_hierarchy, d_solution_time);

    d_convec_oper->applyConvectiveOperator(Q_idx, d_Q_convec_idx);

    for (int level_num = d_coarsest_ln; level_num <= d_finest_ln; ++level_num)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(level_num);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            Pointer<FaceData<NDIM, double> > u_ADV_f_data = patch->getPatchData(d_u_scratch_idx);
            Pointer<CellData<NDIM, double> > u_ADV_c_data = patch->getPatchData(d_u_scratch_idx);
            if (u_ADV_f_data)
            {
                const Pointer<CartesianPatchGeometry<NDIM> > p_geom = patch->getPatchGeometry();
                const double* dx = p_geom->getDx();
                const Box<NDIM>& patch_box = patch->getBox();
                const IntVector<NDIM> patch_lower = patch_box.lower();
                const IntVector<NDIM> patch_upper = patch_box.upper();
                Pointer<CellData<NDIM, double> > Q_data = patch->getPatchData(Q_idx);
                Pointer<CellData<NDIM, double> > Y_data = patch->getPatchData(Y_idx);
                const IntVector<NDIM> Q_data_gcw = Q_data->getGhostCellWidth();

                const IntVector<NDIM> u_data_gcw = u_ADV_f_data->getGhostCellWidth();
                const IntVector<NDIM> Y_data_gcw = Y_data->getGhostCellWidth();
                Pointer<CellData<NDIM, double> > C_data = patch->getPatchData(d_Q_convec_idx);
                const IntVector<NDIM> C_data_gcw = C_data->getGhostCellWidth();
                Pointer<CellData<NDIM, double> > R_data = patch->getPatchData(d_R_idx);
                const IntVector<NDIM> R_data_gcw = R_data->getGhostCellWidth();
                if (d_sqr_root)
                {
#if (NDIM == 2)
                    sqrt_tens_conv_u_f_oper_2d_(dx,
                                                u_ADV_f_data->getPointer(0),
                                                u_ADV_f_data->getPointer(1),
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
                                                d_lambda);
#endif
#if (NDIM == 3)
                    sqrt_tens_conv_u_f_oper_3d_(dx,
                                                u_ADV_f_data->getPointer(0),
                                                u_ADV_f_data->getPointer(1),
                                                u_ADV_f_data->getPointer(2),
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
                                                patch_upper(2),
                                                d_lambda);
#endif
                }
                else if (d_log_conform)
                {
#if (NDIM == 2)
                    log_tens_conv_u_f_oper_2d_(dx,
                                               u_ADV_f_data->getPointer(0),
                                               u_ADV_f_data->getPointer(1),
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
                                               d_lambda);
#endif
#if (NDIM == 3)
                    log_tens_conv_u_f_oper_3d_(dx,
                                               u_ADV_f_data->getPointer(0),
                                               u_ADV_f_data->getPointer(1),
                                               u_ADV_f_data->getPointer(2),
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
                                               patch_upper(2),
                                               d_lambda);
#endif
                }
                else
                {
#if (NDIM == 2)
                    conform_tens_conv_u_f_oper_2d_(dx,
                                                   u_ADV_f_data->getPointer(0),
                                                   u_ADV_f_data->getPointer(1),
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
                                                   d_lambda);
#endif
#if (NDIM == 3)
                    conform_tens_conv_u_f_oper_3d_(dx,
                                                   u_ADV_f_data->getPointer(0),
                                                   u_ADV_f_data->getPointer(1),
                                                   u_ADV_f_data->getPointer(2),
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
                                                   patch_upper(2),
                                                   d_lambda);
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
                                                patch_upper(1),
                                                d_lambda);
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
                                                patch_upper(2),
                                                d_lambda);
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
                                               patch_upper(1),
                                               d_lambda);
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
                                               patch_upper(2),
                                               d_lambda);
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
                                                   patch_upper(1),
                                                   d_lambda);
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
                                                   patch_upper(2),
                                                   d_lambda);
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
    /* Set up the coasen operations. These COARSEN the data (i.e. fills data at coarse interfaces)
     * General process:
     * 1) Set up a coarsen algorithm
     * 2) Register a coarsen operator with the algorithm
     * 3) Fill a coarsen schedule with the coarsen algorithm
     * 4) To actually coarsen data, use coarsen schedule -> coarsen data()
     */
    Pointer<CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();
    Pointer<CoarsenOperator<NDIM> > coarsen_op_Q = grid_geom->lookupCoarsenOperator(d_Q_var, "CONSERVATIVE_COARSEN");
    Pointer<CoarsenOperator<NDIM> > coarsen_op_u =
        grid_geom->lookupCoarsenOperator(d_u_adv_var, "CONSERVATIVE_COARSEN");
    // Step 1) and 2)
    d_coarsen_alg_Q = new CoarsenAlgorithm<NDIM>();
    d_coarsen_alg_Q->registerCoarsen(d_Q_scratch_idx, d_Q_scratch_idx, coarsen_op_Q);
    d_coarsen_scheds_Q.resize(d_finest_ln + 1);
    d_coarsen_alg_u = new CoarsenAlgorithm<NDIM>();
    d_coarsen_alg_u->registerCoarsen(d_u_scratch_idx, d_u_scratch_idx, coarsen_op_u);
    d_coarsen_scheds_u.resize(d_finest_ln + 1);
    // Step 3)
    for (int ln = d_coarsest_ln + 1; ln <= d_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        Pointer<PatchLevel<NDIM> > coarser_level = d_hierarchy->getPatchLevel(ln - 1);
        d_coarsen_scheds_u[ln] = d_coarsen_alg_u->createSchedule(coarser_level, level);
        d_coarsen_scheds_Q[ln] = d_coarsen_alg_Q->createSchedule(coarser_level, level);
    }
    /* Set Refine Algorithms. This interpolates data onto finer grid
     * General process:
     * 1) Set up a refine algorithm
     * 2) Register a refine operation with the algorithm
     * 3) Fill a refine schedule with the refine algorithm
     * 4) Invoke fill data() inside refine schedule
     */
    // Note we only set up refine algorithms for Q here because u has not been set yet.
    Pointer<RefineOperator<NDIM> > refine_op_Q = grid_geom->lookupRefineOperator(d_Q_var, "CONSERVATIVE_LINEAR_REFINE");
    d_ghostfill_alg_Q = new RefineAlgorithm<NDIM>();
    d_ghostfill_alg_Q->registerRefine(d_Q_scratch_idx, in.getComponentDescriptorIndex(0), d_Q_scratch_idx, refine_op_Q);
    if (d_outflow_bdry_extrap_type != "NONE")
        d_ghostfill_strategy_Q = new CartExtrapPhysBdryOp(d_Q_scratch_idx, d_outflow_bdry_extrap_type);
    d_ghostfill_scheds_Q.resize(d_finest_ln + 1);
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        d_ghostfill_scheds_Q[ln] =
            d_ghostfill_alg_Q->createSchedule(level, ln - 1, d_hierarchy, d_ghostfill_strategy_Q);
    }
    // Allocate Patch Data
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        if (!level->checkAllocated(d_Q_scratch_idx)) level->allocatePatchData(d_Q_scratch_idx);

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
        if (level->checkAllocated(d_Q_scratch_idx)) level->deallocatePatchData(d_Q_scratch_idx);
        if (level->checkAllocated(d_Q_convec_idx)) level->deallocatePatchData(d_Q_convec_idx);
        if (level->checkAllocated(d_u_scratch_idx)) level->deallocatePatchData(d_u_scratch_idx);
        if (level->checkAllocated(d_R_idx)) level->deallocatePatchData(d_R_idx);
    }
    // Deallocate the refine algorithm, operator, patch strategy, and schedules.
    d_ghostfill_alg_Q.setNull();
    d_ghostfill_alg_u.setNull();
    d_ghostfill_strategy_Q.setNull();
    d_ghostfill_strategy_u.setNull();
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        d_ghostfill_scheds_u[ln].setNull();
        d_ghostfill_scheds_Q[ln].setNull();
    }
    d_ghostfill_scheds_u.clear();
    d_ghostfill_scheds_Q.clear();
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
