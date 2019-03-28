#include "ibamr/ComplexFluidForcing.h"
extern "C"
{
#if (NDIM == 2)
    void div_tensor_s_to_c_2d_(const double* dx,
                               const double* d_data_0,
                               const double* d_data_1,
                               const int& d_gcw,
                               const double* s_data,
                               const int& s_gcw,
                               const int& ilower0,
                               const int& iupper0,
                               const int& ilower1,
                               const int& iupper1,
                               const double& alpha);
    void div_tensor_c_to_c_2d_(const double* dx,
                               const double* d_data,
                               const int& d_gcw,
                               const double* s_data,
                               const int& s_gcw,
                               const int& ilower0,
                               const int& iupper0,
                               const int& ilower1,
                               const int& iupper1,
                               const double& alpha);
    void div_tensor_s_to_c_weno_2d_(const double*,
                                    const double*,
                                    const int&,
                                    const double*,
                                    const int&,
                                    const int&,
                                    const int&,
                                    const int&,
                                    const int&,
                                    const double&,
                                    const double*,
                                    const int&);
#endif
#if (NDIM == 3)
    void div_tensor_s_to_c_3d_(const double* dx,
                               const double* d_data_0,
                               const double* d_data_1,
                               const double* d_data_2,
                               const int& d_gcw,
                               const double* s_data,
                               const int& s_gcw,
                               const int& ilower0,
                               const int& iupper0,
                               const int& ilower1,
                               const int& iupper1,
                               const int& ilower2,
                               const int& iupper2,
                               const double& alpha);
    void div_tensor_c_to_c_3d_(const double* dx,
                               const double* d_data,
                               const int& d_gcw,
                               const double* s_data,
                               const int& s_gcw,
                               const int& ilower0,
                               const int& iupper0,
                               const int& ilower1,
                               const int& iupper1,
                               const int& ilower2,
                               const int& iupper2,
                               const double& alpha);
#endif
}
// Namespace
namespace IBAMR
{
// Constructor
ComplexFluidForcing::ComplexFluidForcing(const std::string& object_name,
                                         Pointer<Database> input_db,
                                         Pointer<CartGridFunction> u_fcn,
                                         Pointer<CartesianGridGeometry<NDIM> > grid_geometry,
                                         Pointer<AdvDiffSemiImplicitHierarchyIntegrator> adv_diff_integrator,
                                         Pointer<VisItDataWriter<NDIM> > visit_data_writer)
    : d_object_name(object_name),
      d_project_conform(false),
      d_context(NULL),
      init_conds(NULL),
      d_adv_diff_integrator(adv_diff_integrator),
      d_convec_oper(NULL),
      d_convec_oper_type(CENTERED),
      d_W_cc_var(NULL),
      d_conform_var_draw(NULL),
      d_stress_var_draw(NULL),
      d_divW_var_draw(NULL),
      d_conform_idx_draw(-1),
      d_stress_idx_draw(-1),
      d_divW_idx_draw(-1),
      d_conform_draw(true),
      d_stress_draw(false),
      d_divW_draw(false),
      d_lambda(-1),
      d_eta(-1),
      d_W_cc_idx(-1),
      d_W_cc_scratch_idx(-1),
      d_conc_bc_coefs(),
      d_fluid_model("OLDROYDB"),
      d_sqr_root(false),
      d_log_conform(false),
      d_hierarchy(NULL),
      d_max_det(std::numeric_limits<double>::max()),
      d_min_det(-1.0),
      d_log_det(false),
      d_positive_def(true),
      d_error_on_spd(false),
      d_min_norm(-1.0),
      d_max_norm(std::numeric_limits<double>::max()),
      d_log_divW(false),
      d_divW_rel_tag(false),
      d_divW_abs_tag(false),
      d_u_fcn(u_fcn),
      d_u_var(new FaceVariable<NDIM, double>("Complex Fluid Velocity"))
{
    // Set up muParserCartGridFunctions
    init_conds = new muParserCartGridFunction(object_name, input_db->getDatabase("InitialConditions"), grid_geometry);
    // Register Variables and variable context objects.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    d_context = var_db->getContext(d_object_name + "::CONTEXT");
    d_W_cc_var = new CellVariable<NDIM, double>(d_object_name + "::W_cc", NDIM * (NDIM + 1) / 2);
    static const IntVector<NDIM> ghosts_cc = 7;
    // Set up Advection Diffusion Integrator
    d_adv_diff_integrator->registerTransportedQuantity(d_W_cc_var);
    d_adv_diff_integrator->setInitialConditions(d_W_cc_var, init_conds);
    d_adv_diff_integrator->registerAdvectionVelocity(d_u_var);
    d_adv_diff_integrator->setAdvectionVelocityFunction(d_u_var, d_u_fcn);
    d_adv_diff_integrator->setAdvectionVelocity(d_W_cc_var, d_u_var);
    d_W_cc_scratch_idx = var_db->registerVariableAndContext(d_W_cc_var, d_context, ghosts_cc);
    d_W_cc_idx = var_db->registerClonedPatchDataIndex(d_W_cc_var, d_W_cc_scratch_idx);
    getFromInput(input_db, visit_data_writer);
    LocationIndexRobinBcCoefs<NDIM> bc_coef(d_object_name + "::bc_coef", Pointer<Database>(NULL));
    for (int d = 0; d < NDIM; ++d)
    {
        bc_coef.setBoundarySlope(2 * d, 0.0);
        bc_coef.setBoundarySlope(2 * d + 1, 0.0);
    }
    std::vector<RobinBcCoefStrategy<NDIM>*> bc_coefs(NDIM * (NDIM + 1) / 2, &bc_coef);
    // Create boundary conditions for advected materials if not periodic
    const IntVector<NDIM>& periodic_shift = grid_geometry->getPeriodicShift();
    if (periodic_shift.min() <= 0)
    {
        vector<RobinBcCoefStrategy<NDIM>*> conc_bc_coefs(NDIM * (NDIM + 1) / 2);
        d_conc_bc_coefs.resize(NDIM * (NDIM + 1) / 2);
        for (int d = 0; d < NDIM * (NDIM + 1) / 2; ++d)
        {
            ostringstream conc_bc_name;
            conc_bc_name << "c_bc_coef_" << d;
            const string bc_coefs_name = conc_bc_name.str();

            ostringstream bc_coefs_db_name_stream;
            bc_coefs_db_name_stream << "ExtraStressBoundaryConditions_" << d;
            const string bc_coefs_db_name = bc_coefs_db_name_stream.str();

            d_conc_bc_coefs[d] =
                new muParserRobinBcCoefs(bc_coefs_name, input_db->getDatabase(bc_coefs_db_name), grid_geometry);
        }
        d_adv_diff_integrator->setPhysicalBcCoefs(d_W_cc_var, d_conc_bc_coefs);
    }
    d_convec_oper = new AdvDiffComplexFluidConvectiveOperator("ComplexFluidConvectiveOperator",
                                                              d_W_cc_var,
                                                              input_db,
                                                              d_convec_oper_type,
                                                              d_conc_bc_coefs,
                                                              std::vector<RobinBcCoefStrategy<NDIM>*>(NULL));
    d_adv_diff_integrator->setConvectiveOperator(d_W_cc_var, d_convec_oper);
    if (d_fluid_model == "OLDROYDB")
    {
        registerRHSTerm(new OldroydB_RHS("OldroydB_RHS", input_db, d_W_cc_var, adv_diff_integrator));
    }
    else if (d_fluid_model == "GIESEKUS")
    {
        registerRHSTerm(new Giesekus_RHS("Giesekus_RHS", input_db, d_W_cc_var, adv_diff_integrator));
    }
    else if (d_fluid_model == "ROLIEPOLY")
    {
        registerRHSTerm(new RoliePoly_RHS("RoliePoly_RHS", input_db, d_W_cc_var, adv_diff_integrator));
    }
    else
    {
        TBOX_ERROR("Fluid model: " + d_fluid_model + " not known.\n\n");
    }
    if (d_divW_abs_tag || d_divW_rel_tag)
    {
        d_adv_diff_integrator->registerApplyGradientDetectorCallback(&applyGradientDetectorCallback, this);
    }
    return;
} // End Constructor
// Constructor
ComplexFluidForcing::ComplexFluidForcing(const std::string& object_name,
                                         Pointer<Database> input_db,
                                         const Pointer<INSHierarchyIntegrator> fluid_solver,
                                         Pointer<CartesianGridGeometry<NDIM> > grid_geometry,
                                         Pointer<AdvDiffSemiImplicitHierarchyIntegrator> adv_diff_integrator,
                                         Pointer<VisItDataWriter<NDIM> > visit_data_writer)
    : d_object_name(object_name),
      d_project_conform(false),
      d_context(NULL),
      init_conds(NULL),
      d_adv_diff_integrator(adv_diff_integrator),
      d_convec_oper(NULL),
      d_convec_oper_type(CENTERED),
      d_W_cc_var(NULL),
      d_conform_var_draw(NULL),
      d_stress_var_draw(NULL),
      d_divW_var_draw(NULL),
      d_conform_idx_draw(-1),
      d_stress_idx_draw(-1),
      d_divW_idx_draw(-1),
      d_conform_draw(true),
      d_stress_draw(false),
      d_divW_draw(false),
      d_lambda(-1),
      d_eta(-1),
      d_W_cc_idx(-1),
      d_W_cc_scratch_idx(-1),
      d_conc_bc_coefs(),
      d_fluid_model("OLDROYDB"),
      d_sqr_root(false),
      d_log_conform(false),
      d_hierarchy(NULL),
      d_max_det(std::numeric_limits<double>::max()),
      d_min_det(-1.0),
      d_log_det(false),
      d_positive_def(true),
      d_error_on_spd(false),
      d_min_norm(-1.0),
      d_max_norm(std::numeric_limits<double>::max()),
      d_log_divW(false),
      d_divW_rel_tag(false),
      d_divW_abs_tag(false),
      d_u_fcn(NULL),
      d_u_var(NULL)
{
    // Set up muParserCartGridFunctions
    init_conds = new muParserCartGridFunction(object_name, input_db->getDatabase("InitialConditions"), grid_geometry);
    // Register Variables and variable context objects.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    d_context = var_db->getContext(d_object_name + "::CONTEXT");
    d_W_cc_var = new CellVariable<NDIM, double>(d_object_name + "::W_cc", NDIM * (NDIM + 1) / 2);
    static const IntVector<NDIM> ghosts_cc = 7;
    // Set up Advection Diffusion Integrator
    d_adv_diff_integrator->registerTransportedQuantity(d_W_cc_var);
    d_adv_diff_integrator->setInitialConditions(d_W_cc_var, init_conds);
    d_adv_diff_integrator->setAdvectionVelocity(d_W_cc_var, fluid_solver->getAdvectionVelocityVariable());
    d_W_cc_scratch_idx = var_db->registerVariableAndContext(d_W_cc_var, d_context, ghosts_cc);
    d_W_cc_idx = var_db->registerClonedPatchDataIndex(d_W_cc_var, d_W_cc_scratch_idx);
    getFromInput(input_db, visit_data_writer);
    LocationIndexRobinBcCoefs<NDIM> bc_coef(d_object_name + "::bc_coef", Pointer<Database>(NULL));
    for (int d = 0; d < NDIM; ++d)
    {
        bc_coef.setBoundarySlope(2 * d, 0.0);
        bc_coef.setBoundarySlope(2 * d + 1, 0.0);
    }
    std::vector<RobinBcCoefStrategy<NDIM>*> bc_coefs(NDIM * (NDIM + 1) / 2, &bc_coef);
    // Create boundary conditions for advected materials if not periodic
    const IntVector<NDIM>& periodic_shift = grid_geometry->getPeriodicShift();
    if (periodic_shift.min() <= 0)
    {
        vector<RobinBcCoefStrategy<NDIM>*> conc_bc_coefs(NDIM * (NDIM + 1) / 2);
        d_conc_bc_coefs.resize(NDIM * (NDIM + 1) / 2);
        for (int d = 0; d < NDIM * (NDIM + 1) / 2; ++d)
        {
            ostringstream conc_bc_name;
            conc_bc_name << "c_bc_coef_" << d;
            const string bc_coefs_name = conc_bc_name.str();

            ostringstream bc_coefs_db_name_stream;
            bc_coefs_db_name_stream << "ExtraStressBoundaryConditions_" << d;
            const string bc_coefs_db_name = bc_coefs_db_name_stream.str();

            d_conc_bc_coefs[d] =
                new muParserRobinBcCoefs(bc_coefs_name, input_db->getDatabase(bc_coefs_db_name), grid_geometry);
        }
        d_adv_diff_integrator->setPhysicalBcCoefs(d_W_cc_var, d_conc_bc_coefs);
    }
    d_convec_oper = new AdvDiffComplexFluidConvectiveOperator("ComplexFluidConvectiveOperator",
                                                              d_W_cc_var,
                                                              input_db,
                                                              d_convec_oper_type,
                                                              d_conc_bc_coefs,
                                                              fluid_solver->getVelocityBoundaryConditions());
    d_adv_diff_integrator->setConvectiveOperator(d_W_cc_var, d_convec_oper);
    if (d_fluid_model == "OLDROYDB")
    {
        registerRHSTerm(new OldroydB_RHS("OldroydB_RHS", input_db, d_W_cc_var, adv_diff_integrator));
    }
    else if (d_fluid_model == "GIESEKUS")
    {
        registerRHSTerm(new Giesekus_RHS("Giesekus_RHS", input_db, d_W_cc_var, adv_diff_integrator));
    }
    else if (d_fluid_model == "ROLIEPOLY")
    {
        registerRHSTerm(new RoliePoly_RHS("RoliePoly_RHS", input_db, d_W_cc_var, adv_diff_integrator));
    }
    else
    {
        TBOX_ERROR("Fluid model: " + d_fluid_model + " not known.\n\n");
    }
    if (d_divW_abs_tag || d_divW_rel_tag)
    {
        d_adv_diff_integrator->registerApplyGradientDetectorCallback(&applyGradientDetectorCallback, this);
    }
    return;
} // End Constructor
// Destructor
ComplexFluidForcing::~ComplexFluidForcing()
{
    // deallocate draw data...
    int finest_ln = d_hierarchy->getFinestLevelNumber();
    for (int ln = 0; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        if (d_conform_draw && level->checkAllocated(d_conform_idx_draw)) level->deallocatePatchData(d_conform_idx_draw);
        if (d_stress_draw && level->checkAllocated(d_stress_idx_draw)) level->deallocatePatchData(d_stress_idx_draw);
        if ((d_divW_idx_draw > -1) && level->checkAllocated(d_divW_idx_draw))
            level->deallocatePatchData(d_divW_idx_draw);
    }
    return;
} // End Destructor

void
ComplexFluidForcing::getFromInput(const Pointer<Database> input_db, Pointer<VisItDataWriter<NDIM> > visit_data_writer)
{
    d_lambda = input_db->getDouble("relaxation_time");
    d_eta = input_db->getDouble("viscosity");
    d_sqr_root = input_db->getBoolWithDefault("square_root_evolve", false);
    d_log_conform = input_db->getBoolWithDefault("log_conform_evolve", false);
    if (input_db->keyExists("D")) d_adv_diff_integrator->setDiffusionCoefficient(d_W_cc_var, input_db->getDouble("D"));
    d_log_det = input_db->getBoolWithDefault("log_determinant", false);
    d_convec_oper_type = IBAMR::string_to_enum<CFConvectiveOperatorType>(input_db->getStringWithDefault(
        "convective_operator_type", IBAMR::enum_to_string<CFConvectiveOperatorType>(CENTERED)));
    d_fluid_model = input_db->getStringWithDefault("fluid_model", "OldroydB");
    boost::to_upper(d_fluid_model);
    if (input_db->keyExists("error_on_spd"))
    {
        d_error_on_spd = input_db->getBool("error_on_spd");
    }
    if (input_db->keyExists("log_divergence"))
    {
        d_log_divW = input_db->getBool("log_divergence");
    }
    if (input_db->keyExists("project_conformation_tensor"))
    {
        d_project_conform = input_db->getBool("project_conformation_tensor");
    }
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    d_conform_draw = input_db->getBoolWithDefault("output_conformation_tensor", false);
    if (d_conform_draw)
    {
        d_conform_var_draw = new CellVariable<NDIM, double>(d_object_name + "::conform_draw", NDIM * NDIM);
        d_conform_idx_draw = var_db->registerVariableAndContext(d_conform_var_draw, d_context, IntVector<NDIM>(0));
        visit_data_writer->registerPlotQuantity("Conformation_Tensor", "TENSOR", d_conform_idx_draw);
    }
    d_stress_draw = input_db->getBoolWithDefault("output_stress_tensor", true);
    if (d_stress_draw)
    {
        d_stress_var_draw = new CellVariable<NDIM, double>(d_object_name + "::stress_draw", NDIM * NDIM);
        d_stress_idx_draw = var_db->registerVariableAndContext(d_stress_var_draw, d_context, IntVector<NDIM>(0));
        visit_data_writer->registerPlotQuantity("Stress_Tensor", "TENSOR", d_stress_idx_draw);
    }
    d_divW_draw = input_db->getBoolWithDefault("output_divergence", false);
    d_divW_rel_tag = input_db->getBoolWithDefault("divergence_rel_tagging", false);
    d_divW_abs_tag = input_db->getBoolWithDefault("divergence_abs_tagging", false);
    if (d_log_divW || d_divW_draw || d_divW_abs_tag || d_divW_rel_tag)
    {
        d_divW_var_draw = new CellVariable<NDIM, double>(d_object_name + "::divW_draw", NDIM);
        d_divW_idx_draw = var_db->registerVariableAndContext(d_divW_var_draw, d_context, IntVector<NDIM>(0));
        if (d_divW_draw) visit_data_writer->registerPlotQuantity("Stress_Divergence", "VECTOR", d_divW_idx_draw);
    }
    if (d_divW_rel_tag) d_divW_rel_thresh = input_db->getDoubleArray("divergence_rel_thresh");
    if (d_divW_abs_tag) d_divW_abs_thresh = input_db->getDoubleArray("divergence_abs_thresh");
    return;
}

// Time Dependent?
bool
ComplexFluidForcing::isTimeDependent() const
{
    return true;
} // End Time Dependent?

// setDataOnPatchHierarchy
void
ComplexFluidForcing::setDataOnPatchHierarchy(const int data_idx,
                                             Pointer<Variable<NDIM> > var,
                                             Pointer<PatchHierarchy<NDIM> > hierarchy,
                                             const double data_time,
                                             const bool initial_time,
                                             const int coarsest_ln_in,
                                             const int finest_ln_in)
{
    d_hierarchy = hierarchy;
    const int coarsest_ln = (coarsest_ln_in == -1 ? 0 : coarsest_ln_in);
    const int finest_ln = (finest_ln_in == -1 ? hierarchy->getFinestLevelNumber() : finest_ln_in);
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    const IntVector<NDIM> ghosts_cc = 1;

    // Allocate Data to store components of the Complex stress tensor
    for (int level_num = coarsest_ln; level_num <= finest_ln; ++level_num)
    {
        Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_num);
        if (!level->checkAllocated(d_W_cc_scratch_idx)) level->allocatePatchData(d_W_cc_scratch_idx);
        if (d_stress_draw && !level->checkAllocated(d_stress_idx_draw)) level->allocatePatchData(d_stress_idx_draw);
        if (d_conform_draw && !level->checkAllocated(d_conform_idx_draw)) level->allocatePatchData(d_conform_idx_draw);
        if ((d_divW_idx_draw > -1) && !level->checkAllocated(d_divW_idx_draw))
            level->allocatePatchData(d_divW_idx_draw);
        if (!level->checkAllocated(d_W_cc_idx)) level->allocatePatchData(d_W_cc_idx);
    }
    if (initial_time)
    {
        init_conds->setDataOnPatchHierarchy(d_W_cc_idx, d_W_cc_var, hierarchy, data_time, coarsest_ln, finest_ln);
    }
    else
    {
        const int W_current_idx =
            var_db->mapVariableAndContextToIndex(d_W_cc_var, d_adv_diff_integrator->getCurrentContext());
        const int W_new_idx = var_db->mapVariableAndContextToIndex(d_W_cc_var, d_adv_diff_integrator->getNewContext());
        const bool W_new_is_allocated = d_adv_diff_integrator->isAllocatedPatchData(W_new_idx);
        HierarchyDataOpsManager<NDIM>* hier_data_ops_manager = HierarchyDataOpsManager<NDIM>::getManager();
        Pointer<HierarchyDataOpsReal<NDIM, double> > hier_cc_data_ops =
            hier_data_ops_manager->getOperationsDouble(d_W_cc_var, hierarchy, true);
        if (d_adv_diff_integrator->getCurrentCycleNumber() == 0 || !W_new_is_allocated)
        {
            hier_cc_data_ops->copyData(d_W_cc_idx, W_current_idx);
        }
        else
        {
            hier_cc_data_ops->linearSum(d_W_cc_idx, 0.5, W_current_idx, 0.5, W_new_idx);
        }
    }

    HierarchyDataOpsManager<NDIM>* hier_data_ops_manager = HierarchyDataOpsManager<NDIM>::getManager();
    Pointer<HierarchyDataOpsReal<NDIM, double> > hier_cc_data_ops =
        hier_data_ops_manager->getOperationsDouble(d_W_cc_var, hierarchy);
    if (d_sqr_root)
    {
        squareMatrix(d_W_cc_scratch_idx, d_W_cc_var, hierarchy, data_time, initial_time, coarsest_ln, finest_ln);
    }
    else if (d_log_conform)
    {
        exponentiateMatrix(d_W_cc_scratch_idx, d_W_cc_var, hierarchy, data_time, initial_time, coarsest_ln, finest_ln);
    }
    else
    {
        if (d_project_conform)
        {
            projectMatrix(d_W_cc_idx, d_W_cc_var, hierarchy, data_time, initial_time, coarsest_ln, finest_ln);
        }
        hier_cc_data_ops->copyData(d_W_cc_scratch_idx, d_W_cc_idx);
    }

    typedef HierarchyGhostCellInterpolation::InterpolationTransactionComponent InterpolationTransactionComponent;
    std::vector<InterpolationTransactionComponent> ghost_cell_components(1);
    ghost_cell_components[0] = InterpolationTransactionComponent(
        d_W_cc_scratch_idx, "CONSERVATIVE_LINEAR_REFINE", true, "CONSERVATIVE_COARSEN", "LINEAR", false, d_conc_bc_coefs, NULL);
    HierarchyGhostCellInterpolation ghost_fill_op;
    ghost_fill_op.initializeOperatorState(ghost_cell_components, hierarchy);
    ghost_fill_op.fillData(data_time);

    Pointer<FaceVariable<NDIM, double> > u_var = d_adv_diff_integrator->getAdvectionVelocity(d_W_cc_var);
    const int u_idx = var_db->mapVariableAndContextToIndex(u_var, d_adv_diff_integrator->getCurrentContext());
/*
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            Pointer<CellData<NDIM, double> > Q_data = patch->getPatchData(d_W_cc_scratch_idx);
            Pointer<FaceData<NDIM, double> > u_adv_data = patch->getPatchData(u_idx);
            AdvDiffPhysicalBoundaryUtilities::setPhysicalBoundaryConditions(
                Q_data, u_adv_data, patch, d_conc_bc_coefs, data_time, "LINEAR", false);
        }
    }
*/
    // Check to ensure conformation tensor is positive definite
    d_positive_def = true;
    for (int level_num = coarsest_ln; level_num <= finest_ln; ++level_num)
    {
        checkPositiveDefiniteOnLevel(
            d_W_cc_scratch_idx, d_W_cc_var, hierarchy->getPatchLevel(level_num), data_time, initial_time);
    }
    int temp = d_positive_def ? 1 : 0;
    int output;
    MPI_Allreduce(&temp, &output, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    plog << "Conformation tensor is " << (output == 1 ? "SPD" : "NOT SPD") << "\n";
    if (d_error_on_spd && output == 0)
        TBOX_ERROR("ERROR: \n "
                   << "CONFORMATION TENSOR IS NOT SPD! \n\n");

    // Check max and min determinant
    if (d_log_det)
    {
        d_max_det = 0.0;
        d_min_det = 1e10;
        for (int level_num = coarsest_ln; level_num <= finest_ln; ++level_num)
        {
            findDeterminant(
                d_W_cc_scratch_idx, d_W_cc_var, hierarchy->getPatchLevel(level_num), data_time, initial_time);
        }
        double temp_max, temp_min;
        MPI_Allreduce(&d_max_det, &temp_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allreduce(&d_min_det, &temp_min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        plog << "Largest det:  " << temp_max << "\n";
        plog << "Smallest det: " << temp_min << "\n";
    }

    // Compute Div W on each patch level
    d_min_norm = 1.0e10;
    d_max_norm = -1.0e10;
    for (int level_num = coarsest_ln; level_num <= finest_ln; ++level_num)
    {
        setDataOnPatchLevel(data_idx, var, hierarchy->getPatchLevel(level_num), data_time, initial_time);
    }
    // Output largest and smallest max norm of Div W
    if (d_log_divW || d_divW_rel_tag)
    {
        double temp_max, temp_min;
        MPI_Allreduce(&d_max_norm, &temp_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allreduce(&d_min_norm, &temp_min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        plog << "Largest max norm of Div W:  " << temp_max << "\n";
        plog << "Smallest max norm of Div W: " << temp_min << "\n";
    }

    for (int level_num = coarsest_ln; level_num <= finest_ln; ++level_num)
    {
        Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_num);
        if (level->checkAllocated(d_W_cc_scratch_idx)) level->deallocatePatchData(d_W_cc_scratch_idx);
    }
    return;
} // End setDataOnPatchHierarchy

void
ComplexFluidForcing::setDataOnPatchLevel(const int data_idx,
                                         Pointer<Variable<NDIM> > var,
                                         Pointer<PatchLevel<NDIM> > level,
                                         const double data_time,
                                         const bool initial_time)
{
    if (initial_time)
    {
        if (!level->checkAllocated(d_W_cc_scratch_idx)) level->allocatePatchData(d_W_cc_scratch_idx);
        if (d_conform_draw && !level->checkAllocated(d_conform_idx_draw)) level->allocatePatchData(d_conform_idx_draw);
        if (d_stress_draw && !level->checkAllocated(d_stress_idx_draw)) level->allocatePatchData(d_stress_idx_draw);
        if ((d_divW_idx_draw > -1) && !level->checkAllocated(d_divW_idx_draw)) level->allocatePatchData(d_divW_idx_draw);
        init_conds->setDataOnPatchLevel(d_W_cc_scratch_idx, d_W_cc_var, level, data_time, initial_time);
    }
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = level->getPatch(p());
        setDataOnPatch(data_idx, var, patch, data_time, initial_time, level);
    }
    return;
}

void
ComplexFluidForcing::setDataOnPatch(const int data_idx,
                                    Pointer<Variable<NDIM> > /*var*/,
                                    Pointer<Patch<NDIM> > patch,
                                    const double /*data_time*/,
                                    const bool initial_time,
                                    Pointer<PatchLevel<NDIM> > /*patch_level*/)
{
    const Box<NDIM>& patch_box = patch->getBox();
    const Pointer<CartesianPatchGeometry<NDIM> > p_geom = patch->getPatchGeometry();
    const double* dx = p_geom->getDx();
    Pointer<CellData<NDIM, double> > W_cc_data = patch->getPatchData(d_W_cc_scratch_idx);
    Pointer<CellData<NDIM, double> > divW_draw_data(NULL);
    if (d_divW_idx_draw > -1) divW_draw_data = patch->getPatchData(d_divW_idx_draw);
    if (d_log_divW || d_divW_draw || d_divW_abs_tag || d_divW_rel_tag) divW_draw_data->fillAll(0.0);
    if (d_conform_draw)
    {
        Pointer<CellData<NDIM, double> > conform_data_draw = patch->getPatchData(d_conform_idx_draw);
#if (NDIM == 2)
        conform_data_draw->copyDepth(0, *W_cc_data, 0);
        conform_data_draw->copyDepth(1, *W_cc_data, 2);
        conform_data_draw->copyDepth(2, *W_cc_data, 2);
        conform_data_draw->copyDepth(3, *W_cc_data, 1);
#endif
#if (NDIM == 3)
        conform_data_draw->copyDepth(0, *W_cc_data, 0);
        conform_data_draw->copyDepth(1, *W_cc_data, 5);
        conform_data_draw->copyDepth(2, *W_cc_data, 4);
        conform_data_draw->copyDepth(3, *W_cc_data, 5);
        conform_data_draw->copyDepth(4, *W_cc_data, 1);
        conform_data_draw->copyDepth(5, *W_cc_data, 3);
        conform_data_draw->copyDepth(6, *W_cc_data, 4);
        conform_data_draw->copyDepth(7, *W_cc_data, 3);
        conform_data_draw->copyDepth(8, *W_cc_data, 2);
#endif
    }
    if (d_stress_draw)
    {
        Pointer<CellData<NDIM, double> > stress_data_draw = patch->getPatchData(d_stress_idx_draw);
        for (CellIterator<NDIM> ci(patch_box); ci; ci++)
        {
            CellIndex<NDIM> idx = *ci;
#if (NDIM == 2)
            (*stress_data_draw)(idx, 0) = d_eta / d_lambda * ((*W_cc_data)(idx, 0) - 1.0);
            (*stress_data_draw)(idx, 1) = d_eta / d_lambda * (*W_cc_data)(idx, 2);
            (*stress_data_draw)(idx, 2) = d_eta / d_lambda * (*W_cc_data)(idx, 2);
            (*stress_data_draw)(idx, 3) = d_eta / d_lambda * ((*W_cc_data)(idx, 1) - 1.0);
#endif
#if (NDIM == 3)
            (*stress_data_draw)(idx, 0) = d_eta / d_lambda * ((*W_cc_data)(idx, 0) - 1.0);
            (*stress_data_draw)(idx, 1) = d_eta / d_lambda * (*W_cc_data)(idx, 5);
            (*stress_data_draw)(idx, 2) = d_eta / d_lambda * (*W_cc_data)(idx, 4);
            (*stress_data_draw)(idx, 3) = d_eta / d_lambda * (*W_cc_data)(idx, 5);
            (*stress_data_draw)(idx, 4) = d_eta / d_lambda * ((*W_cc_data)(idx, 1) - 1.0);
            (*stress_data_draw)(idx, 5) = d_eta / d_lambda * (*W_cc_data)(idx, 3);
            (*stress_data_draw)(idx, 6) = d_eta / d_lambda * (*W_cc_data)(idx, 4);
            (*stress_data_draw)(idx, 7) = d_eta / d_lambda * (*W_cc_data)(idx, 3);
            (*stress_data_draw)(idx, 8) = d_eta / d_lambda * ((*W_cc_data)(idx, 2) - 1.0);
#endif
        }
    }
    Pointer<SideData<NDIM, double> > divW_sc_data = patch->getPatchData(data_idx);
    Pointer<CellData<NDIM, double> > divW_cc_data = patch->getPatchData(data_idx);
    if (divW_sc_data)
    {
        divW_sc_data->fillAll(0.0);
        if (initial_time) return;
        const IntVector<NDIM> W_cc_ghosts = W_cc_data->getGhostCellWidth();
        const IntVector<NDIM> divW_sc_ghosts = divW_sc_data->getGhostCellWidth();
        // Compute prefactor of divergence of stress.
        // If we are computing the conformation tensor, we need to convert to stress
        // NOTE: This term depends on the particular fluid model (Currently works for OldroydB and Giesekus, maybe
        // RoliePoly?) Perhaps a better way would be to convert to stress, then take divergence (takes care of nonlinear
        // models).
        double alpha = d_eta / d_lambda;
        const IntVector<NDIM>& patch_lower = patch_box.lower();
        const IntVector<NDIM>& patch_upper = patch_box.upper();
#if (NDIM == 2)
        //     div_tensor_s_to_c_weno_2d_(divW_sc_data->getPointer(0), divW_sc_data->getPointer(1),
        //     divW_sc_ghosts.max(),
        //                                W_cc_data->getPointer(0), W_cc_ghosts.min(), patch_lower(0), patch_upper(0),
        //                                patch_lower(1), patch_upper(1), alpha, dx, 3);
        div_tensor_s_to_c_2d_(dx,
                              divW_sc_data->getPointer(0),
                              divW_sc_data->getPointer(1),
                              divW_sc_ghosts.min(),
                              W_cc_data->getPointer(0),
                              W_cc_ghosts.min(),
                              patch_lower(0),
                              patch_upper(0),
                              patch_lower(1),
                              patch_upper(1),
                              alpha);
#endif
#if (NDIM == 3)
        div_tensor_s_to_c_3d_(dx,
                              divW_sc_data->getPointer(0),
                              divW_sc_data->getPointer(1),
                              divW_sc_data->getPointer(2),
                              divW_sc_ghosts.min(),
                              W_cc_data->getPointer(0),
                              W_cc_ghosts.min(),
                              patch_lower(0),
                              patch_upper(0),
                              patch_lower(1),
                              patch_upper(1),
                              patch_lower(2),
                              patch_upper(2),
                              alpha);
#endif
        if (d_log_divW || d_divW_draw || d_divW_abs_tag || d_divW_rel_tag)
        {
            for (CellIterator<NDIM> ci(patch_box); ci; ci++)
            {
                CellIndex<NDIM> idx = *ci;
                SideIndex<NDIM> f_x_n = SideIndex<NDIM>(idx, 0, 0);
                SideIndex<NDIM> f_x_p = SideIndex<NDIM>(idx, 0, 1);
                SideIndex<NDIM> f_y_n = SideIndex<NDIM>(idx, 1, 0);
                SideIndex<NDIM> f_y_p = SideIndex<NDIM>(idx, 1, 1);
                double max_norm;
#if (NDIM == 2)
                max_norm = std::max<double>(std::fabs(0.5 * ((*divW_sc_data)(f_x_n) + (*divW_sc_data)(f_x_p))),
                                            std::fabs(0.5 * ((*divW_sc_data)(f_y_n) + (*divW_sc_data)(f_y_p))));
                (*divW_draw_data)(idx, 0) = 0.5 * ((*divW_sc_data)(f_x_n) + (*divW_sc_data)(f_x_p));
                (*divW_draw_data)(idx, 1) = 0.5 * ((*divW_sc_data)(f_y_n) + (*divW_sc_data)(f_y_p));
#endif
#if (NDIM == 3)
                SideIndex<NDIM> f_z_n = SideIndex<NDIM>(idx, 2, 0);
                SideIndex<NDIM> f_z_p = SideIndex<NDIM>(idx, 2, 1);
                max_norm = std::max<double>(
                    std::fabs(0.5 * ((*divW_sc_data)(f_x_n) + (*divW_sc_data)(f_x_p))),
                    std::max<double>(std::fabs(0.5 * ((*divW_sc_data)(f_y_n) + (*divW_sc_data)(f_y_p))),
                                     std::fabs(0.5 * ((*divW_sc_data)(f_z_n) + (*divW_sc_data)(f_z_p)))));
                (*divW_draw_data)(idx, 0) = 0.5 * ((*divW_sc_data)(f_x_n) + (*divW_sc_data)(f_x_p));
                (*divW_draw_data)(idx, 1) = 0.5 * ((*divW_sc_data)(f_y_n) + (*divW_sc_data)(f_y_p));
                (*divW_draw_data)(idx, 2) = 0.5 * ((*divW_sc_data)(f_z_n) + (*divW_sc_data)(f_z_p));
#endif
                if (d_log_divW || d_divW_rel_tag)
                {
                    d_min_norm = (d_min_norm > max_norm) ? max_norm : d_min_norm;
                    d_max_norm = (d_max_norm > max_norm) ? d_max_norm : max_norm;
                }
            }
        }
    }
    else if (divW_cc_data)
    {
        divW_cc_data->fillAll(0.0);
        if (initial_time) return;
        const IntVector<NDIM> W_cc_ghosts = W_cc_data->getGhostCellWidth();
        const IntVector<NDIM> divW_cc_ghosts = divW_cc_data->getGhostCellWidth();
        // Compute prefactor of divergence of stress.
        // If we are computing the conformation tensor, we need to convert to stress
        // NOTE: This term depends on the particular fluid model (Currently works for OldroydB and Giesekus, maybe
        // RoliePoly?) Perhaps a better way would be to convert to stress, then take divergence (takes care of nonlinear
        // models).
        double alpha = d_eta / d_lambda;
        const IntVector<NDIM>& patch_lower = patch_box.lower();
        const IntVector<NDIM>& patch_upper = patch_box.upper();
#if (NDIM == 2)
        div_tensor_c_to_c_2d_(dx,
                              divW_cc_data->getPointer(0),
                              divW_cc_ghosts.min(),
                              W_cc_data->getPointer(0),
                              W_cc_ghosts.min(),
                              patch_lower(0),
                              patch_upper(0),
                              patch_lower(1),
                              patch_upper(1),
                              alpha);
#endif
#if (NDIM == 3)
        div_tensor_c_to_c_3d_(dx,
                              divW_cc_data->getPointer(0),
                              divW_cc_ghosts.min(),
                              W_cc_data->getPointer(0),
                              W_cc_ghosts.min(),
                              patch_lower(0),
                              patch_upper(0),
                              patch_lower(1),
                              patch_upper(1),
                              patch_lower(2),
                              patch_upper(2),
                              alpha);
#endif
        for (CellIterator<NDIM> ci(patch_box); ci; ci++)
        {
            CellIndex<NDIM> idx = *ci;
            double max_norm;
#if (NDIM == 2)
            max_norm = std::max<double>(std::fabs((*divW_cc_data)(idx, 0)), std::fabs((*divW_cc_data)(idx, 1)));
#endif
#if (NDIM == 3)
            max_norm = std::max<double>(
                std::fabs((*divW_cc_data)(idx, 0)),
                std::max<double>(std::fabs((*divW_cc_data)(idx, 1)), std::fabs((*divW_cc_data)(idx, 2))));
#endif
            if (d_log_divW || d_divW_rel_tag)
            {
                d_min_norm = (d_min_norm > max_norm) ? max_norm : d_min_norm;
                d_max_norm = (d_max_norm > max_norm) ? d_max_norm : max_norm;
            }
        }
        if (d_divW_draw || d_divW_abs_tag || d_divW_rel_tag) divW_draw_data->copy(*divW_cc_data);
    }
} // end setDataOnPatch

void
ComplexFluidForcing::registerRHSTerm(Pointer<RHS_Operator> rhs)
{
    d_convec_oper->registerRHSTerm(rhs);
    return;
}

void
ComplexFluidForcing::checkPositiveDefiniteOnLevel(const int data_idx,
                                                  const Pointer<Variable<NDIM> > var,
                                                  const Pointer<PatchLevel<NDIM> > level,
                                                  const double data_time,
                                                  const bool initial_time)
{
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        checkPositiveDefiniteOnPatch(data_idx, var, level->getPatch(p()), data_time, initial_time);
    }
    return;
}

void
ComplexFluidForcing::checkPositiveDefiniteOnPatch(const int data_idx,
                                                  const Pointer<Variable<NDIM> > /*var*/,
                                                  const Pointer<Patch<NDIM> > patch,
                                                  const double /*data_time*/,
                                                  const bool initial_time)
{
    const Box<NDIM>& box = patch->getBox();
    const Pointer<PatchGeometry<NDIM> > p_geom = patch->getPatchGeometry();
    if (initial_time) return;
    Pointer<CellData<NDIM, double> > s_data = patch->getPatchData(data_idx);
    for (CellIterator<NDIM> it(box); it; it++)
    {
        CellIndex<NDIM> ci = *it;
#if (NDIM == 2)
        Eigen::Matrix2d tens(Eigen::Matrix2d::Zero());
        tens(0, 0) = (*s_data)(ci, 0);
        tens(0, 1) = tens(1, 0) = (*s_data)(ci, 2);
        tens(1, 1) = (*s_data)(ci, 1);
        Eigen::LLT<Eigen::Matrix2d> llt;
#endif
#if (NDIM == 3)
        Eigen::Matrix3d tens(Eigen::Matrix3d::Zero());
        tens(0, 0) = (*s_data)(ci, 0);
        tens(1, 1) = (*s_data)(ci, 1);
        tens(2, 2) = (*s_data)(ci, 2);
        tens(0, 1) = tens(1, 0) = (*s_data)(ci, 5);
        tens(0, 2) = tens(2, 0) = (*s_data)(ci, 4);
        tens(1, 2) = tens(2, 1) = (*s_data)(ci, 3);
        Eigen::LLT<Eigen::Matrix3d> llt;
#endif
        llt.compute(tens);
        if (llt.info() == Eigen::NumericalIssue)
        {
            d_positive_def = false;
        }
    }
    return;
}

void
ComplexFluidForcing::squareMatrix(const int data_idx,
                                  const Pointer<Variable<NDIM> > /*var*/,
                                  const Pointer<PatchHierarchy<NDIM> > hierarchy,
                                  const double /*data_time*/,
                                  const bool initial_time,
                                  const int coarsest_ln,
                                  const int finest_ln)
{
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Box<NDIM>& box = patch->getBox();
            const Pointer<PatchGeometry<NDIM> > p_geom = patch->getPatchGeometry();
            if (initial_time) return;
            Pointer<CellData<NDIM, double> > half_data = patch->getPatchData(d_W_cc_idx);
            Pointer<CellData<NDIM, double> > tens_data = patch->getPatchData(data_idx);
            for (CellIterator<NDIM> it(box); it; it++)
            {
                CellIndex<NDIM> i = *it;
#if (NDIM == 2)
                (*tens_data)(i, 0) = (*half_data)(i, 0) * (*half_data)(i, 0) + (*half_data)(i, 2) * (*half_data)(i, 2);
                (*tens_data)(i, 1) = (*half_data)(i, 2) * (*half_data)(i, 2) + (*half_data)(i, 1) * (*half_data)(i, 1);
                (*tens_data)(i, 2) = (*half_data)(i, 0) * (*half_data)(i, 2) + (*half_data)(i, 1) * (*half_data)(i, 2);
#endif
#if (NDIM == 3)
                (*tens_data)(i, 0) = (*half_data)(i, 0) * (*half_data)(i, 0) + (*half_data)(i, 4) * (*half_data)(i, 4) +
                                     (*half_data)(i, 5) * (*half_data)(i, 5);
                (*tens_data)(i, 1) = (*half_data)(i, 1) * (*half_data)(i, 1) + (*half_data)(i, 3) * (*half_data)(i, 3) +
                                     (*half_data)(i, 5) * (*half_data)(i, 5);
                (*tens_data)(i, 2) = (*half_data)(i, 2) * (*half_data)(i, 2) + (*half_data)(i, 3) * (*half_data)(i, 3) +
                                     (*half_data)(i, 4) * (*half_data)(i, 4);
                (*tens_data)(i, 3) = (*half_data)(i, 1) * (*half_data)(i, 3) + (*half_data)(i, 2) * (*half_data)(i, 3) +
                                     (*half_data)(i, 4) * (*half_data)(i, 5);
                (*tens_data)(i, 4) = (*half_data)(i, 0) * (*half_data)(i, 4) + (*half_data)(i, 2) * (*half_data)(i, 4) +
                                     (*half_data)(i, 3) * (*half_data)(i, 5);
                (*tens_data)(i, 5) = (*half_data)(i, 3) * (*half_data)(i, 4) + (*half_data)(i, 0) * (*half_data)(i, 5) +
                                     (*half_data)(i, 1) * (*half_data)(i, 5);
#endif
            }
        }
    }
    return;
}

void
ComplexFluidForcing::findDeterminant(const int data_idx,
                                     const Pointer<Variable<NDIM> > /*var*/,
                                     const Pointer<PatchLevel<NDIM> > level,
                                     const double /*data_time*/,
                                     const bool /*initial_time*/)
{
    double det;
    for (PatchLevel<NDIM>::Iterator i(level); i; i++)
    {
        Pointer<Patch<NDIM> > patch = level->getPatch(i());
        const Box<NDIM>& box = patch->getBox();
        Pointer<CellData<NDIM, double> > data = patch->getPatchData(data_idx);
        for (CellIterator<NDIM> it(box); it; it++)
        {
            CellIndex<NDIM> i = *it;
#if (NDIM == 2)
            det = (*data)(i, 0) * (*data)(i, 1) - (*data)(i, 2) * (*data)(i, 2);
#endif
#if (NDIM == 3)
            det = (*data)(i, 0) * (*data)(i, 1) * (*data)(i, 2) - (*data)(i, 0) * (*data)(i, 3) * (*data)(i, 3) -
                  (*data)(i, 1) * (*data)(i, 4) * (*data)(i, 4) + 2.0 * (*data)(i, 3) * (*data)(i, 4) * (*data)(i, 5) -
                  (*data)(i, 2) * (*data)(i, 5) * (*data)(i, 5);
#endif
            d_max_det = det > d_max_det ? det : d_max_det;
            d_min_det = det < d_min_det ? det : d_min_det;
        }
    }
    return;
}

void
ComplexFluidForcing::exponentiateMatrix(const int data_idx,
                                        const Pointer<Variable<NDIM> > /*var*/,
                                        const Pointer<PatchHierarchy<NDIM> > hierarchy,
                                        const double /*data_time*/,
                                        const bool /*initial_time*/,
                                        const int coarsest_ln,
                                        const int finest_ln)
{
    for (int ln = coarsest_ln; ln <= finest_ln; ln++)
    {
        Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Box<NDIM>& box = patch->getBox();
            Pointer<CellData<NDIM, double> > data = patch->getPatchData(data_idx);
            Pointer<CellData<NDIM, double> > log_data = patch->getPatchData(d_W_cc_idx);
            for (CellIterator<NDIM> it(box); it; it++)
            {
                CellIndex<NDIM> i = *it;
#if (NDIM == 2)
                Eigen::Matrix2d mat;
                mat(0, 0) = (*log_data)(i, 0);
                mat(1, 0) = (*log_data)(i, 2);
                mat(0, 1) = (*log_data)(i, 2);
                mat(1, 1) = (*log_data)(i, 1);
                mat = mat.exp();
                (*data)(i, 0) = mat(0, 0);
                (*data)(i, 1) = mat(1, 1);
                (*data)(i, 2) = mat(0, 1);
#endif
#if (NDIM == 3)
                Eigen::Matrix3d mat;
                mat(0, 0) = (*log_data)(i, 0);
                mat(1, 1) = (*log_data)(i, 1);
                mat(2, 2) = (*log_data)(i, 2);
                mat(1, 2) = mat(2, 1) = (*log_data)(i, 3);
                mat(0, 2) = mat(2, 0) = (*log_data)(i, 4);
                mat(0, 1) = mat(1, 0) = (*log_data)(i, 5);
                mat = mat.exp();
                (*data)(i, 0) = mat(0, 0);
                (*data)(i, 1) = mat(1, 1);
                (*data)(i, 2) = mat(2, 2);
                (*data)(i, 3) = mat(1, 2);
                (*data)(i, 4) = mat(0, 2);
                (*data)(i, 5) = mat(0, 1);
#endif
            }
        }
    }
    return;
}

void
ComplexFluidForcing::projectMatrix(const int data_idx,
                                   const Pointer<Variable<NDIM> > /*var*/,
                                   const Pointer<PatchHierarchy<NDIM> > hierarchy,
                                   const double /*data_time*/,
                                   const bool initial_time,
                                   const int coarsest_ln,
                                   const int finest_ln)
{
    for (int ln = coarsest_ln; ln <= finest_ln; ln++)
    {
        Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Box<NDIM>& box = patch->getBox();
            if (initial_time) return;
            Pointer<CellData<NDIM, double> > data = patch->getPatchData(data_idx);
            for (CellIterator<NDIM> it(box); it; it++)
            {
                CellIndex<NDIM> i = *it;
#if (NDIM == 2)
                Eigen::Matrix2d c;
                Eigen::SelfAdjointEigenSolver<Matrix2d> eigs;
                c(0, 0) = (*data)(i, 0);
                c(1, 1) = (*data)(i, 1);
                c(0, 1) = c(1, 0) = (*data)(i, 2);
                eigs.computeDirect(c);
                Eigen::Matrix2d eig_vals(Eigen::Matrix2d::Zero());
                eig_vals(0, 0) = eigs.eigenvalues()(0);
                eig_vals(1, 1) = eigs.eigenvalues()(1);
                eig_vals(0, 0) = max(eig_vals(0, 0), 0.0);
                eig_vals(1, 1) = max(eig_vals(1, 1), 0.0);
                Eigen::Matrix2d eig_vecs = eigs.eigenvectors();
#endif
#if (NDIM == 3)
                Eigen::Matrix3d c;
                Eigen::SelfAdjointEigenSolver<Matrix3d> eigs;
                c(0, 0) = (*data)(i, 0);
                c(1, 1) = (*data)(i, 1);
                c(2, 2) = (*data)(i, 2);
                c(2, 1) = c(1, 2) = (*data)(i, 3);
                c(0, 2) = c(2, 0) = (*data)(i, 4);
                c(0, 1) = c(1, 0) = (*data)(i, 5);
                eigs.computeDirect(c);
                Eigen::Matrix3d eig_vals(Eigen::Matrix3d::Zero());
                eig_vals(0, 0) = eigs.eigenvalues()(0);
                eig_vals(1, 1) = eigs.eigenvalues()(1);
                eig_vals(2, 2) = eigs.eigenvalues()(2);
                eig_vals(0, 0) = max(eig_vals(0, 0), 0.0);
                eig_vals(1, 1) = max(eig_vals(1, 1), 0.0);
                eig_vals(2, 2) = max(eig_vals(2, 2), 0.0);
                Eigen::Matrix3d eig_vecs = eigs.eigenvectors();
#endif
                c = eig_vecs * eig_vals * eig_vecs.transpose();
#if (NDIM == 2)
                (*data)(i, 0) = c(0, 0);
                (*data)(i, 1) = c(1, 1);
                (*data)(i, 2) = c(0, 1);
#endif
#if (NDIM == 3)
                (*data)(i, 0) = c(0, 0);
                (*data)(i, 1) = c(1, 1);
                (*data)(i, 2) = c(2, 2);
                (*data)(i, 3) = c(1, 2);
                (*data)(i, 4) = c(0, 2);
                (*data)(i, 5) = c(0, 1);
#endif
            }
        }
    }
    return;
}

void
ComplexFluidForcing::applyGradientDetector(Pointer<BasePatchHierarchy<NDIM> > hierarchy,
                                           int level_number,
                                           double /*error_data_time*/,
                                           int tag_index,
                                           bool initial_time,
                                           bool /*richardson_extrapolation_too*/)
{
    if (initial_time) return;
    Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_number);
    double divW_rel_thresh = 0.0;
    if (d_divW_rel_thresh.size() > 0)
        divW_rel_thresh = d_divW_rel_thresh[std::max(std::min(level_number, d_divW_rel_thresh.size() - 1), 0)];
    double divW_abs_thresh = 0.0;
    if (d_divW_abs_thresh.size() > 0)
        divW_abs_thresh = d_divW_abs_thresh[std::max(std::min(level_number, d_divW_abs_thresh.size() - 1), 0)];
    if (divW_rel_thresh > 0.0 || divW_abs_thresh > 0.0)
    {
        double thresh = std::numeric_limits<double>::max();
        if (divW_abs_thresh > 0.0) thresh = std::min(thresh, divW_abs_thresh);
        if (divW_rel_thresh > 0.0) thresh = std::min(thresh, divW_rel_thresh * d_max_norm);
        thresh += sqrt(std::numeric_limits<double>::epsilon());
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            Pointer<CellData<NDIM, double> > W_data = patch->getPatchData(d_divW_idx_draw);
            if (!W_data) continue;
            Pointer<CellData<NDIM, int> > tag_data = patch->getPatchData(tag_index);
            const Box<NDIM>& box = patch->getBox();
            double norm;
            for (CellIterator<NDIM> ic(box); ic; ic++)
            {
                const CellIndex<NDIM>& i = ic();
                norm = 0.0;
                for (int d = 0; d < NDIM; ++d) norm += (*W_data)(i, d) * (*W_data)(i, d);
                norm = sqrt(norm);
                if (norm > thresh) (*tag_data)(i) = 1;
            }
        }
    }
    return;
}
void
ComplexFluidForcing::applyGradientDetectorCallback(Pointer<BasePatchHierarchy<NDIM> > hierarchy,
                                                   int level_number,
                                                   double error_data_time,
                                                   int tag_index,
                                                   bool initial_time,
                                                   bool richardson_extrapolation_too,
                                                   void* ctx)
{
    ComplexFluidForcing* object = static_cast<ComplexFluidForcing*>(ctx);
    object->applyGradientDetector(
        hierarchy, level_number, error_data_time, tag_index, initial_time, richardson_extrapolation_too);
    return;
}
} // namespace IBAMR
