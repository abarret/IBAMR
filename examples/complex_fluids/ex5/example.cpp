// ---------------------------------------------------------------------
//
// Copyright (c) 2019 - 2022 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#include <ibamr/CFIIMethod.h>
#include <ibamr/CFINSForcing.h>
#include <ibamr/CFOneSidedUpperConvectiveOperator.h>
#include <ibamr/IBExplicitHierarchyIntegrator.h>
#include <ibamr/IIMethod.h>
#include <ibamr/INSCollocatedHierarchyIntegrator.h>
#include <ibamr/INSStaggeredHierarchyIntegrator.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/LEInteractor.h>
#include <ibtk/libmesh_utilities.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

#include "libmesh/face_tri3_subdivision.h"
#include "libmesh/mesh_modification.h"
#include "libmesh/mesh_refinement.h"
#include "libmesh/mesh_subdivision_support.h"
#include "libmesh/mesh_tools.h"
#include <libmesh/boundary_info.h>
#include <libmesh/boundary_mesh.h>
#include <libmesh/equation_systems.h>
#include <libmesh/exact_solution.h>
#include <libmesh/exodusII_io.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/mesh_triangle_interface.h>

#include <petscsys.h>

#include <boost/multi_array.hpp>

#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <LoadBalancer.h>
#include <SAMRAI_config.h>
#include <StandardTagAndInitialize.h>

#include <ibamr/app_namespaces.h>

// Local includes
#include "ExtraStressExact.h"
#include "VelocityInit.h"

// Elasticity model data.
namespace ModelData
{
// Tether (penalty) force functions.
static double kappa_s = 1.0e6;
static double eta_s = 0.0;
static double rho = 1.0;
static double dx = -1.0;
static bool ERROR_ON_MOVE = false;
void
tether_force_function(VectorValue<double>& F,
                      const libMesh::VectorValue<double>& /*n*/,
                      const libMesh::VectorValue<double>& /*N*/,
                      const TensorValue<double>& /*FF*/,
                      const libMesh::Point& x,
                      const libMesh::Point& X,
                      Elem* const /*elem*/,
                      unsigned short int /*side*/,
                      const vector<const vector<double>*>& var_data,
                      const vector<const vector<VectorValue<double> >*>& /*grad_var_data*/,
                      double /*time*/,
                      void* /*ctx*/)
{
    const std::vector<double>& U = *var_data[0];
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        F(d) = kappa_s * (X(d) - x(d)) - eta_s * U[d];
    }
    std::vector<double> d(NDIM);
    d[0] = std::abs(x(0) - X(0));
    d[1] = std::abs(X(1) - x(1));
    if (ERROR_ON_MOVE && ((d[0] > 0.25 * dx) || (d[1] > 0.25 * dx)))
    {
        TBOX_ERROR("Structure has moved too much.\n");
    }

    return;
} // tether_force_function
} // namespace ModelData
using namespace ModelData;

static std::string SIG_IN_ERR_SYS_NAME = "SIG_IN_ERR";
static std::string SIG_OUT_ERR_SYS_NAME = "SIG_OUT_ERR";

// Function prototypes
void setInsideOfChannel(int ls_idx,
                        Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
                        const double y_low,
                        const double y_up,
                        const double theta,
                        double data_time);

void computeErrorAndPlotEulerian(Pointer<INSHierarchyIntegrator> time_integrator,
                                 Pointer<AdvDiffHierarchyIntegrator> adv_diff_integrator,
                                 Pointer<CFINSForcing> cf_forcing,
                                 Pointer<VisItDataWriter<NDIM> > visit_data_writer,
                                 Pointer<PatchHierarchy<NDIM> > hierarchy,
                                 Pointer<CartGridFunction> u_exact_fcn,
                                 Pointer<CartGridFunction> sig_exact_fcn,
                                 const int u_err_idx,
                                 Pointer<SideVariable<NDIM, double> > u_err_var,
                                 const int sig_xx_err_idx,
                                 const int sig_yy_err_idx,
                                 const int sig_xy_err_idx,
                                 const int u_draw_err_idx,
                                 const int sig_xx_exact_idx,
                                 const int sig_yy_exact_idx,
                                 const int sig_xy_exact_idx,
                                 const int u_draw_exact_idx,
                                 Pointer<CellVariable<NDIM, double> > u_draw_err_var,
                                 const int iteration_num,
                                 const double loop_time);

void computeErrorAndPlotLagrangian(const std::unique_ptr<ExodusII_IO>& lower_io,
                                   const std::string& lower_filename,
                                   EquationSystems* lower_eq_sys,
                                   const std::unique_ptr<ExodusII_IO>& upper_io,
                                   const std::string& upper_filename,
                                   EquationSystems* upper_eq_sys,
                                   Pointer<ExtraStressExact> sig_in_fcn,
                                   Pointer<ExtraStressExact> sig_out_fcn,
                                   Pointer<VelocityInit> vel_fcn,
                                   const double loop_time,
                                   const int iter_num,
                                   const int which_sys);

/*******************************************************************************
 * For each run, the input filename and restart information (if needed) must   *
 * be given on the command line.  For non-restarted case, command line is:     *
 *                                                                             *
 *    executable <input file name>                                             *
 *                                                                             *
 * For restarted run, command line is:                                         *
 *                                                                             *
 *    executable <input file name> <restart directory> <restart number>        *
 *                                                                             *
 *******************************************************************************/
int
main(int argc, char* argv[])
{
    // Initialize IBAMR and libraries. Deinitialization is handled by this object as well.
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);
    const LibMeshInit& init = ibtk_init.getLibMeshInit();

    { // cleanup dynamically allocated objects prior to shutdown

        // Parse command line options, set some standard options from the input
        // file, initialize the restart database (if this is a restarted run),
        // and enable file logging.
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "IB.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();

        // Get various standard options set in the input file.
        const bool dump_viz_data = app_initializer->dumpVizData();
        const int viz_dump_interval = app_initializer->getVizDumpInterval();
        const bool uses_visit = dump_viz_data && app_initializer->getVisItDataWriter();
#ifdef LIBMESH_HAVE_EXODUS_API
        const bool uses_exodus = dump_viz_data && !app_initializer->getExodusIIFilename().empty();
#else
        const bool uses_exodus = false;
        if (!app_initializer->getExodusIIFilename().empty())
        {
            plog << "WARNING: libMesh was compiled without Exodus support, so no "
                 << "Exodus output will be written in this program.\n";
        }
#endif
        const string lower_filename = app_initializer->getExodusIIFilename("lower");
        const string upper_filename = app_initializer->getExodusIIFilename("upper");

        const bool dump_restart_data = app_initializer->dumpRestartData();
        const int restart_dump_interval = app_initializer->getRestartDumpInterval();
        const string restart_dump_dirname = app_initializer->getRestartDumpDirectory();

        const bool dump_postproc_data = app_initializer->dumpPostProcessingData();
        const int postproc_data_dump_interval = app_initializer->getPostProcessingDataDumpInterval();
        const string postproc_data_dump_dirname = app_initializer->getPostProcessingDataDumpDirectory();
        if (dump_postproc_data && (postproc_data_dump_interval > 0) && !postproc_data_dump_dirname.empty())
        {
            Utilities::recursiveMkdir(postproc_data_dump_dirname);
        }

        const bool dump_timer_data = app_initializer->dumpTimerData();
        const int timer_dump_interval = app_initializer->getTimerDumpInterval();

        // Create a simple FE mesh.
        Mesh solid_mesh(init.comm(), NDIM);
        dx = input_db->getDouble("DX");
        ERROR_ON_MOVE = input_db->getBool("ERROR_ON_MOVE");
        const double theta = input_db->getDouble("THETA"); // Channel angle.
        const double L =
            input_db->getDouble("LX"); // Channel Length. Note this is different than domain length if slope != 0.
        const double y_low = input_db->getDouble("CHANNEL_LOW");
        const double y_up = input_db->getDouble("CHANNEL_UP");
        const double x_low = input_db->getDouble("X_LOW");
        const double x_up = input_db->getDouble("X_UP");
        const double mfac = input_db->getDouble("MFAC");
        const double ds = mfac * dx;
        string elem_type = input_db->getString("ELEM_TYPE");
        ReplicatedMesh lower_mesh(init.comm(), NDIM), upper_mesh(init.comm(), NDIM);
        MeshTools::Generation::build_line(lower_mesh,
                                          static_cast<int>(ceil(L / ds)),
                                          0.0 + 0.0 * dx,
                                          L - 0.0 * dx,
                                          Utility::string_to_enum<ElemType>(elem_type));
        MeshTools::Generation::build_line(upper_mesh,
                                          static_cast<int>(ceil(L / ds)),
                                          0.0 + 0.0 * dx,
                                          L - 0.0 * dx,
                                          Utility::string_to_enum<ElemType>(elem_type));

        for (MeshBase::node_iterator it = lower_mesh.nodes_begin(); it != lower_mesh.nodes_end(); ++it)
        {
            Node* n = *it;
            libMesh::Point& X = *n;
            X(1) = y_low + std::tan(theta) * X(0);
        }

        for (MeshBase::node_iterator it = upper_mesh.nodes_begin(); it != upper_mesh.nodes_end(); ++it)
        {
            Node* n = *it;
            libMesh::Point& X = *n;
            X(1) = y_up + std::tan(theta) * X(0);
        }

        for (MeshBase::element_iterator it = lower_mesh.elements_begin(); it != lower_mesh.elements_end(); ++it)
        {
            Elem* elem = *it;
            bool elem_near_bdry = false;
            for (unsigned int n_id = 0; n_id < elem->n_nodes(); ++n_id)
            {
                Node& node = elem->node_ref(n_id);
                if (node(0) < (x_low + 0.4) || node(0) > (x_up - 0.4)) elem_near_bdry = true;
            }
            if (elem_near_bdry) elem->subdomain_id() = 2;
        }

        for (MeshBase::element_iterator it = upper_mesh.elements_begin(); it != upper_mesh.elements_end(); ++it)
        {
            Elem* elem = *it;
            bool elem_near_bdry = false;
            for (unsigned int n_id = 0; n_id < elem->n_nodes(); ++n_id)
            {
                Node& node = elem->node_ref(n_id);
                if (node(0) < (x_low + 0.4) || node(0) > (x_up - 0.4)) elem_near_bdry = true;
            }
            if (elem_near_bdry) elem->subdomain_id() = 2;
        }

        lower_mesh.prepare_for_use();
        upper_mesh.prepare_for_use();

        std::vector<MeshBase*> meshes = { &lower_mesh, &upper_mesh };

        kappa_s = input_db->getDouble("KAPPA_S");
        eta_s = input_db->getDouble("ETA_S");
        rho = input_db->getDouble("RHO");

        // Create major algorithm and data objects that comprise the
        // application.  These objects are configured from the input database
        // and, if this is a restarted run, from the restart database.
        Pointer<INSHierarchyIntegrator> navier_stokes_integrator;
        const string solver_type = app_initializer->getComponentDatabase("Main")->getString("solver_type");
        if (solver_type == "STAGGERED")
        {
            navier_stokes_integrator = new INSStaggeredHierarchyIntegrator(
                "INSStaggeredHierarchyIntegrator",
                app_initializer->getComponentDatabase("INSStaggeredHierarchyIntegrator"));
        }
        else if (solver_type == "COLLOCATED")
        {
            navier_stokes_integrator = new INSCollocatedHierarchyIntegrator(
                "INSCollocatedHierarchyIntegrator",
                app_initializer->getComponentDatabase("INSCollocatedHierarchyIntegrator"));
        }
        else
        {
            TBOX_ERROR("Unsupported solver type: " << solver_type << "\n"
                                                   << "Valid options are: COLLOCATED, STAGGERED");
        }
        Pointer<AdvDiffSemiImplicitHierarchyIntegrator> adv_diff_integrator =
            new AdvDiffSemiImplicitHierarchyIntegrator(
                "AdvDiffSemiImplicitHierarchyIntegrator",
                app_initializer->getComponentDatabase("AdvDiffSemiImplicitHierarchyIntegrator"));
        navier_stokes_integrator->registerAdvDiffHierarchyIntegrator(adv_diff_integrator);
        bool use_cfiim = input_db->getBool("USE_CFIIM");
        Pointer<IIMethod> ib_method_ops;
        if (use_cfiim)
        {
            ib_method_ops =
                new CFIIMethod("CFIIMethod",
                               app_initializer->getComponentDatabase("IIMethod"),
                               meshes,
                               app_initializer->getComponentDatabase("GriddingAlgorithm")->getInteger("max_levels"));
        }
        else
        {
            ib_method_ops =
                new IIMethod("IIMethod",
                             app_initializer->getComponentDatabase("IIMethod"),
                             meshes,
                             app_initializer->getComponentDatabase("GriddingAlgorithm")->getInteger("max_levels"));
        }
        Pointer<IBHierarchyIntegrator> time_integrator =
            new IBExplicitHierarchyIntegrator("IBHierarchyIntegrator",
                                              app_initializer->getComponentDatabase("IBHierarchyIntegrator"),
                                              ib_method_ops,
                                              navier_stokes_integrator);
        Pointer<CartesianGridGeometry<NDIM> > grid_geometry = new CartesianGridGeometry<NDIM>(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        Pointer<PatchHierarchy<NDIM> > patch_hierarchy = new PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);
        Pointer<StandardTagAndInitialize<NDIM> > error_detector =
            new StandardTagAndInitialize<NDIM>("StandardTagAndInitialize",
                                               time_integrator,
                                               app_initializer->getComponentDatabase("StandardTagAndInitialize"));
        Pointer<BergerRigoutsos<NDIM> > box_generator = new BergerRigoutsos<NDIM>();
        Pointer<LoadBalancer<NDIM> > load_balancer =
            new LoadBalancer<NDIM>("LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
        Pointer<GriddingAlgorithm<NDIM> > gridding_algorithm =
            new GriddingAlgorithm<NDIM>("GriddingAlgorithm",
                                        app_initializer->getComponentDatabase("GriddingAlgorithm"),
                                        error_detector,
                                        box_generator,
                                        load_balancer);

        // Configure the II solver.
        ib_method_ops->initializeFEEquationSystems();
        std::vector<int> vars(NDIM);
        for (unsigned int d = 0; d < NDIM; ++d) vars[d] = d;
        std::string vel_sys_name;
        vel_sys_name = IIMethod::VELOCITY_SYSTEM_NAME;
        vector<SystemData> sys_data(1, SystemData(vel_sys_name, vars));
        IIMethod::LagSurfaceForceFcnData body_fcn_data(tether_force_function, sys_data);
        ib_method_ops->registerLagSurfaceForceFunction(body_fcn_data, 0);
        ib_method_ops->registerLagSurfaceForceFunction(body_fcn_data, 1);
        EquationSystems* lower_eq_sys = ib_method_ops->getFEDataManager(0)->getEquationSystems();
        EquationSystems* upper_eq_sys = ib_method_ops->getFEDataManager(1)->getEquationSystems();

        // Create Eulerian initial condition specification objects.
        Pointer<VelocityInit> u_fcn = new VelocityInit("Velocity", app_initializer->getComponentDatabase("Params"), 0);
        Pointer<VelocityInit> v_fcn = new VelocityInit("Velocity", app_initializer->getComponentDatabase("Params"), 1);
        navier_stokes_integrator->registerPhysicalBoundaryConditions({ u_fcn.getPointer(), v_fcn.getPointer() });
        navier_stokes_integrator->registerVelocityInitialConditions(u_fcn);

        if (input_db->keyExists("PressureInitialConditions"))
        {
            Pointer<CartGridFunction> p_init = new muParserCartGridFunction(
                "p_init", app_initializer->getComponentDatabase("PressureInitialConditions"), grid_geometry);
            navier_stokes_integrator->registerPressureInitialConditions(p_init);
        }

        // Set up visualization plot file writers.
        Pointer<VisItDataWriter<NDIM> > visit_data_writer = app_initializer->getVisItDataWriter();
        if (uses_visit)
        {
            time_integrator->registerVisItDataWriter(visit_data_writer);
        }

        // Create level set object
        Pointer<CellVariable<NDIM, double> > ls_var = new CellVariable<NDIM, double>("ls");
        auto var_db = VariableDatabase<NDIM>::getDatabase();
        const int ls_idx = var_db->registerVariableAndContext(ls_var, var_db->getContext("LS"), 5 /*ghosts*/);

        // Create Eulerian body force function specification objects.
        Pointer<CFINSForcing> polymericStressForcing;
        if (input_db->keyExists("ForcingFunction"))
        {
            Pointer<CartGridFunction> f_fcn = new muParserCartGridFunction(
                "f_fcn", app_initializer->getComponentDatabase("ForcingFunction"), grid_geometry);
            time_integrator->registerBodyForceFunction(f_fcn);
        }
        if (input_db->keyExists("ComplexFluid"))
        {
            polymericStressForcing = new CFINSForcing("PolymericStressForcing",
                                                      app_initializer->getComponentDatabase("ComplexFluid"),
                                                      navier_stokes_integrator,
                                                      grid_geometry,
                                                      adv_diff_integrator,
                                                      visit_data_writer);
            time_integrator->registerBodyForceFunction(polymericStressForcing);

            Pointer<CFOneSidedUpperConvectiveOperator> cf_upper_convec_op =
                polymericStressForcing->getUpperConvectiveOperator();
            if (cf_upper_convec_op) cf_upper_convec_op->setLSIdx(ls_idx);
        }

        Pointer<ExtraStressExact> sig_fcn =
            new ExtraStressExact("ExtraStress", app_initializer->getComponentDatabase("Params"), 0, true);
        Pointer<ExtraStressExact> sig_yy_fcn =
            new ExtraStressExact("ExtraStress", app_initializer->getComponentDatabase("Params"), 1, true);
        Pointer<ExtraStressExact> sig_xy_fcn =
            new ExtraStressExact("ExtraStress", app_initializer->getComponentDatabase("Params"), 2, true);
        Pointer<ExtraStressExact> sig_out_fcn =
            new ExtraStressExact("ExtraStress", app_initializer->getComponentDatabase("Params"), 0, false);
        if (polymericStressForcing)
        {
            Pointer<hier::Variable<NDIM> > sig_var = polymericStressForcing->getVariable();
            adv_diff_integrator->setInitialConditions(sig_var, sig_fcn);
            adv_diff_integrator->setPhysicalBcCoefs(
                sig_var, { sig_fcn.getPointer(), sig_yy_fcn.getPointer(), sig_xy_fcn.getPointer() });
        }

        Pointer<CartGridFunction> p_exact_fcn =
            new muParserCartGridFunction("p_exact", app_initializer->getComponentDatabase("PExact"), grid_geometry);
        Pointer<SideVariable<NDIM, double> > u_err_var = new SideVariable<NDIM, double>("U_err");
        Pointer<CellVariable<NDIM, double> > sig_xx_err_var = new CellVariable<NDIM, double>("SIG_XX_ERR"),
                                             sig_yy_err_var = new CellVariable<NDIM, double>("SIG_YY_ERR"),
                                             sig_xy_err_var = new CellVariable<NDIM, double>("SIG_XY_ERR"),
                                             p_err_var = new CellVariable<NDIM, double>("P_ERR"),
                                             u_err_draw_var = new CellVariable<NDIM, double>("U_DRAW_ERR", NDIM);
        Pointer<VariableContext> err_ctx = var_db->getContext("Error");
        Pointer<VariableContext> exa_ctx = var_db->getContext("Exact");
        int u_err_idx = var_db->registerVariableAndContext(u_err_var, err_ctx, 0),
            sig_xx_err_idx = var_db->registerVariableAndContext(sig_xx_err_var, err_ctx, 0),
            sig_yy_err_idx = var_db->registerVariableAndContext(sig_yy_err_var, err_ctx, 0),
            sig_xy_err_idx = var_db->registerVariableAndContext(sig_xy_err_var, err_ctx, 0),
            u_draw_err_idx = var_db->registerVariableAndContext(u_err_draw_var, err_ctx, 0);
        int sig_xx_exact_idx = var_db->registerVariableAndContext(sig_xx_err_var, exa_ctx, 0),
            sig_yy_exact_idx = var_db->registerVariableAndContext(sig_yy_err_var, exa_ctx, 0),
            sig_xy_exact_idx = var_db->registerVariableAndContext(sig_xy_err_var, exa_ctx, 0),
            u_draw_exact_idx = var_db->registerVariableAndContext(u_err_draw_var, exa_ctx, 0);
        visit_data_writer->registerPlotQuantity("U_err", "VECTOR", u_draw_err_idx);
        visit_data_writer->registerPlotQuantity("sig_XX_err", "SCALAR", sig_xx_err_idx);
        visit_data_writer->registerPlotQuantity("sig_YY_err", "SCALAR", sig_yy_err_idx);
        visit_data_writer->registerPlotQuantity("sig_XY_err", "SCALAR", sig_xy_err_idx);
        visit_data_writer->registerPlotQuantity("U_exact", "VECTOR", u_draw_exact_idx);
        visit_data_writer->registerPlotQuantity("sig_XX_exact", "SCALAR", sig_xx_exact_idx);
        visit_data_writer->registerPlotQuantity("sig_YY_exact", "SCALAR", sig_yy_exact_idx);
        visit_data_writer->registerPlotQuantity("sig_XY_exact", "SCALAR", sig_xy_exact_idx);
        visit_data_writer->registerPlotQuantity("LS", "SCALAR", ls_idx);
        // Generate errors for sigma.
        {
            System& sig_err_in_sys = lower_eq_sys->add_system("Explicit", SIG_IN_ERR_SYS_NAME);
            sig_err_in_sys.add_variable("sig_err_in_xx");
            sig_err_in_sys.add_variable("sig_err_in_yy");
            sig_err_in_sys.add_variable("sig_err_in_xy");
            System& sig_err_out_sys = lower_eq_sys->add_system("Explicit", SIG_OUT_ERR_SYS_NAME);
            sig_err_out_sys.add_variable("sig_err_out_xx");
            sig_err_out_sys.add_variable("sig_err_out_yy");
            sig_err_out_sys.add_variable("sig_err_out_xy");
        }
        {
            System& sig_err_in_sys = upper_eq_sys->add_system("Explicit", SIG_IN_ERR_SYS_NAME);
            sig_err_in_sys.add_variable("sig_err_in_xx");
            sig_err_in_sys.add_variable("sig_err_in_yy");
            sig_err_in_sys.add_variable("sig_err_in_xy");
            System& sig_err_out_sys = upper_eq_sys->add_system("Explicit", SIG_OUT_ERR_SYS_NAME);
            sig_err_out_sys.add_variable("sig_err_out_xx");
            sig_err_out_sys.add_variable("sig_err_out_yy");
            sig_err_out_sys.add_variable("sig_err_out_xy");
        }

        // Note for regrids, we need to tell the integrator to call setInsideCylinder
        auto ls_data = std::make_tuple(ls_idx, y_low, y_up, theta);
        auto hierarchy_callback =
            [](Pointer<BasePatchHierarchy<NDIM> > hierarchy, double data_time, bool initial_time, void* ctx) {
                if (initial_time) return;
                const auto& ls_data = *static_cast<std::tuple<int, double, double, double>*>(ctx);
                setInsideOfChannel(std::get<0>(ls_data),
                                   hierarchy,
                                   std::get<1>(ls_data),
                                   std::get<2>(ls_data),
                                   std::get<3>(ls_data),
                                   data_time);
            };
        time_integrator->registerRegridHierarchyCallback(hierarchy_callback, static_cast<void*>(&ls_data));

        // Create an extra sigma index with large ghost width
        const int sig_scr_idx =
            var_db->registerVariableAndContext(polymericStressForcing->getVariable(), var_db->getContext("Ctx"), 5);

        std::unique_ptr<ExodusII_IO> lower_io(uses_exodus ? new ExodusII_IO(lower_mesh) : NULL);
        std::unique_ptr<ExodusII_IO> upper_io(uses_exodus ? new ExodusII_IO(upper_mesh) : NULL);

        // Initialize hierarchy configuration and data on all patches.
        ib_method_ops->initializeFEData();
        time_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);

        // Deallocate initialization objects.
        app_initializer.setNull();

        // Print the input database contents to the log file.
        plog << "Input database:\n";
        input_db->printClassData(plog);

        // Write out initial visualization data.
        int iteration_num = time_integrator->getIntegratorStep();
        double loop_time = time_integrator->getIntegratorTime();
        if (dump_viz_data)
        {
            pout << "\n\nWriting visualization files...\n\n";
            setInsideOfChannel(ls_idx, patch_hierarchy, y_low, y_up, theta, loop_time);
            if (uses_visit)
            {
                // Compute error and setup for plotting
                computeErrorAndPlotEulerian(navier_stokes_integrator,
                                            adv_diff_integrator,
                                            polymericStressForcing,
                                            visit_data_writer,
                                            patch_hierarchy,
                                            u_fcn,
                                            sig_fcn,
                                            u_err_idx,
                                            u_err_var,
                                            sig_xx_err_idx,
                                            sig_yy_err_idx,
                                            sig_xy_err_idx,
                                            u_draw_err_idx,
                                            sig_xx_exact_idx,
                                            sig_yy_exact_idx,
                                            sig_xy_exact_idx,
                                            u_draw_exact_idx,
                                            u_err_draw_var,
                                            iteration_num,
                                            loop_time);
            }
            if (uses_exodus)
            {
                computeErrorAndPlotLagrangian(lower_io,
                                              lower_filename,
                                              lower_eq_sys,
                                              upper_io,
                                              upper_filename,
                                              upper_eq_sys,
                                              sig_fcn,
                                              sig_out_fcn,
                                              u_fcn,
                                              loop_time,
                                              iteration_num / viz_dump_interval + 1,
                                              0);
            }
        }

        // Main time step loop.
        double loop_time_end = time_integrator->getEndTime();
        double dt = 0.0;
        while (!IBTK::rel_equal_eps(loop_time, loop_time_end) && time_integrator->stepsRemaining())
        {
            iteration_num = time_integrator->getIntegratorStep();
            loop_time = time_integrator->getIntegratorTime();

            // Set level set data for grad(u)
            setInsideOfChannel(ls_idx, patch_hierarchy, y_low, y_up, theta, loop_time);
            // Compute the stress jumps
            const int sig_idx = var_db->mapVariableAndContextToIndex(polymericStressForcing->getVariable(),
                                                                     adv_diff_integrator->getCurrentContext());
            for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
            {
                Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
                level->allocatePatchData(sig_scr_idx, loop_time);
            }
            // Fill in ghost cells
            using ITC = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
            ITC ghost_cell_comp(sig_scr_idx,
                                sig_idx,
                                "CONSERVATIVE_LINEAR_REFINE",
                                false,
                                "CONSERVATIVE_COARSEN",
                                "LINEAR",
                                true,
                                adv_diff_integrator->getPhysicalBcCoefs(polymericStressForcing->getVariable()));
            HierarchyGhostCellInterpolation ghost_cell_fill;
            ghost_cell_fill.initializeOperatorState(
                ghost_cell_comp, patch_hierarchy, 0, patch_hierarchy->getFinestLevelNumber());
            ghost_cell_fill.fillData(loop_time);
            Pointer<CFIIMethod> cf_ii_method = ib_method_ops;
            if (cf_ii_method)
            {
                cf_ii_method->computeStressJumps(sig_scr_idx, ls_idx, loop_time, 0);
                cf_ii_method->computeStressJumps(sig_scr_idx, ls_idx, loop_time, 1);
            }
            for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
            {
                Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
                level->deallocatePatchData(sig_scr_idx);
            }

            pout << "\n";
            pout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
            pout << "At beginning of timestep # " << iteration_num << "\n";
            pout << "Simulation time is " << loop_time << "\n";

            dt = time_integrator->getMaximumTimeStepSize();
            time_integrator->advanceHierarchy(dt);
            loop_time += dt;

            pout << "\n";
            pout << "At end       of timestep # " << iteration_num << "\n";
            pout << "Simulation time is " << loop_time << "\n";
            pout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
            pout << "\n";

            // At specified intervals, write visualization and restart files,
            // print out timer data, and store hierarchy data for post
            // processing.
            iteration_num += 1;
            const bool last_step = !time_integrator->stepsRemaining();
            if (dump_viz_data && (iteration_num % viz_dump_interval == 0 || last_step))
            {
                pout << "\nWriting visualization files...\n\n";
                if (uses_visit)
                {
                    // Compute error and setup for plotting
                    computeErrorAndPlotEulerian(navier_stokes_integrator,
                                                adv_diff_integrator,
                                                polymericStressForcing,
                                                visit_data_writer,
                                                patch_hierarchy,
                                                u_fcn,
                                                sig_fcn,
                                                u_err_idx,
                                                u_err_var,
                                                sig_xx_err_idx,
                                                sig_yy_err_idx,
                                                sig_xy_err_idx,
                                                u_draw_err_idx,
                                                sig_xx_exact_idx,
                                                sig_yy_exact_idx,
                                                sig_xy_exact_idx,
                                                u_draw_exact_idx,
                                                u_err_draw_var,
                                                iteration_num,
                                                loop_time);
                }
                if (uses_exodus)
                {
                    computeErrorAndPlotLagrangian(lower_io,
                                                  lower_filename,
                                                  lower_eq_sys,
                                                  upper_io,
                                                  upper_filename,
                                                  upper_eq_sys,
                                                  sig_fcn,
                                                  sig_out_fcn,
                                                  u_fcn,
                                                  loop_time,
                                                  iteration_num / viz_dump_interval + 1,
                                                  0);
                }
            }
            if (dump_restart_data && (iteration_num % restart_dump_interval == 0 || last_step))
            {
                pout << "\nWriting restart files...\n\n";
                RestartManager::getManager()->writeRestartFile(restart_dump_dirname, iteration_num);
            }
            if (dump_timer_data && (iteration_num % timer_dump_interval == 0 || last_step))
            {
                pout << "\nWriting timer data...\n\n";
                TimerManager::getManager()->print(plog);
            }
        }
    } // cleanup dynamically allocated objects prior to shutdown
} // main

void
computeErrorAndPlotEulerian(Pointer<INSHierarchyIntegrator> time_integrator,
                            Pointer<AdvDiffHierarchyIntegrator> adv_diff_integrator,
                            Pointer<CFINSForcing> cf_forcing,
                            Pointer<VisItDataWriter<NDIM> > visit_data_writer,
                            Pointer<PatchHierarchy<NDIM> > hierarchy,
                            Pointer<CartGridFunction> u_exact_fcn,
                            Pointer<CartGridFunction> sig_exact_fcn,
                            const int u_err_idx,
                            Pointer<SideVariable<NDIM, double> > u_err_var,
                            const int sig_xx_err_idx,
                            const int sig_yy_err_idx,
                            const int sig_xy_err_idx,
                            const int u_draw_err_idx,
                            const int sig_xx_exact_idx,
                            const int sig_yy_exact_idx,
                            const int sig_xy_exact_idx,
                            const int u_draw_exact_idx,
                            Pointer<CellVariable<NDIM, double> > u_draw_err_var,
                            const int iteration_num,
                            const double loop_time)
{
    int coarsest_ln = 0;
    int finest_ln = hierarchy->getFinestLevelNumber();
    HierarchyMathOps hier_math_ops("HierMathOps", hierarchy, coarsest_ln, finest_ln);
    HierarchySideDataOpsReal<NDIM, double> hier_sc_data_ops(hierarchy, coarsest_ln, finest_ln);
    HierarchyCellDataOpsReal<NDIM, double> hier_cc_data_ops(hierarchy, coarsest_ln, finest_ln);

    auto var_db = VariableDatabase<NDIM>::getDatabase();
    const int u_cur_idx = var_db->mapVariableAndContextToIndex(time_integrator->getVelocityVariable(),
                                                               time_integrator->getCurrentContext());
    const int sig_cur_idx =
        var_db->mapVariableAndContextToIndex(cf_forcing->getVariable(), adv_diff_integrator->getCurrentContext());
    const int sig_scr_idx =
        var_db->mapVariableAndContextToIndex(cf_forcing->getVariable(), adv_diff_integrator->getScratchContext());
    bool deallocate_sig_scr_after = !adv_diff_integrator->isAllocatedPatchData(sig_scr_idx, coarsest_ln, finest_ln);

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
        if (!level->checkAllocated(u_err_idx)) level->allocatePatchData(u_err_idx, loop_time);
        if (!level->checkAllocated(sig_xx_err_idx)) level->allocatePatchData(sig_xx_err_idx, loop_time);
        if (!level->checkAllocated(sig_yy_err_idx)) level->allocatePatchData(sig_yy_err_idx, loop_time);
        if (!level->checkAllocated(sig_xy_err_idx)) level->allocatePatchData(sig_xy_err_idx, loop_time);
        if (!level->checkAllocated(u_draw_err_idx)) level->allocatePatchData(u_draw_err_idx, loop_time);
        if (!level->checkAllocated(sig_xx_exact_idx)) level->allocatePatchData(sig_xx_exact_idx, loop_time);
        if (!level->checkAllocated(sig_yy_exact_idx)) level->allocatePatchData(sig_yy_exact_idx, loop_time);
        if (!level->checkAllocated(sig_xy_exact_idx)) level->allocatePatchData(sig_xy_exact_idx, loop_time);
        if (!level->checkAllocated(u_draw_exact_idx)) level->allocatePatchData(u_draw_exact_idx, loop_time);
        if (!level->checkAllocated(sig_scr_idx)) level->allocatePatchData(sig_scr_idx, loop_time);
    }

    // Fill in exact answer
    u_exact_fcn->setDataOnPatchHierarchy(
        u_err_idx, time_integrator->getVelocityVariable(), hierarchy, loop_time, false, coarsest_ln, finest_ln);
    sig_exact_fcn->setDataOnPatchHierarchy(
        sig_scr_idx, cf_forcing->getVariable(), hierarchy, loop_time, false, coarsest_ln, finest_ln);
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            Pointer<CellData<NDIM, double> > sig_exact_data = patch->getPatchData(sig_scr_idx);
            Pointer<CellData<NDIM, double> > sig_xx_data = patch->getPatchData(sig_xx_exact_idx);
            Pointer<CellData<NDIM, double> > sig_yy_data = patch->getPatchData(sig_yy_exact_idx);
            Pointer<CellData<NDIM, double> > sig_xy_data = patch->getPatchData(sig_xy_exact_idx);
            sig_xx_data->copyDepth(0, *sig_exact_data, 0);
            sig_yy_data->copyDepth(0, *sig_exact_data, 1);
            sig_xy_data->copyDepth(0, *sig_exact_data, 2);
        }
    }
    hier_math_ops.interp(u_draw_exact_idx, u_draw_err_var, u_err_idx, u_err_var, nullptr, loop_time, false);

    // Now compute the error
    hier_sc_data_ops.subtract(u_err_idx, u_err_idx, u_cur_idx);
    hier_cc_data_ops.subtract(sig_scr_idx, sig_scr_idx, sig_cur_idx);
    // Copy the sig error into the individual patch indices
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            Pointer<CellData<NDIM, double> > sig_err_data = patch->getPatchData(sig_scr_idx);
            Pointer<CellData<NDIM, double> > sig_xx_data = patch->getPatchData(sig_xx_err_idx);
            Pointer<CellData<NDIM, double> > sig_yy_data = patch->getPatchData(sig_yy_err_idx);
            Pointer<CellData<NDIM, double> > sig_xy_data = patch->getPatchData(sig_xy_err_idx);
            sig_xx_data->copyDepth(0, *sig_err_data, 0);
            sig_yy_data->copyDepth(0, *sig_err_data, 1);
            sig_xy_data->copyDepth(0, *sig_err_data, 2);
        }
    }

    // Now compute error norms
    const int wgt_cc_idx = hier_math_ops.getCellWeightPatchDescriptorIndex();
    const int wgt_sc_idx = hier_math_ops.getSideWeightPatchDescriptorIndex();
    pout << "  Computing error norms at time " << loop_time << "\n";
    pout << "  U L1 error:  " << hier_sc_data_ops.L1Norm(u_err_idx, wgt_sc_idx) << "\n";
    pout << "  U L2 error:  " << hier_sc_data_ops.L2Norm(u_err_idx, wgt_sc_idx) << "\n";
    pout << "  U max error: " << hier_sc_data_ops.maxNorm(u_err_idx, wgt_sc_idx) << "\n";
    pout << "\n";
    pout << "  Stress_xx L1 error:  " << hier_cc_data_ops.L1Norm(sig_xx_err_idx, wgt_cc_idx) << "\n";
    pout << "  Stress_xx L2 error:  " << hier_cc_data_ops.L2Norm(sig_xx_err_idx, wgt_cc_idx) << "\n";
    pout << "  Stress_xx max error: " << hier_cc_data_ops.maxNorm(sig_xx_err_idx, wgt_cc_idx) << "\n";
    pout << "\n";
    pout << "  Stress_yy L1 error:  " << hier_cc_data_ops.L1Norm(sig_yy_err_idx, wgt_cc_idx) << "\n";
    pout << "  Stress_yy L2 error:  " << hier_cc_data_ops.L2Norm(sig_yy_err_idx, wgt_cc_idx) << "\n";
    pout << "  Stress_yy max error: " << hier_cc_data_ops.maxNorm(sig_yy_err_idx, wgt_cc_idx) << "\n";
    pout << "\n";
    pout << "  Stress_xy L1 error:  " << hier_cc_data_ops.L1Norm(sig_xy_err_idx, wgt_cc_idx) << "\n";
    pout << "  Stress_xy L2 error:  " << hier_cc_data_ops.L2Norm(sig_xy_err_idx, wgt_cc_idx) << "\n";
    pout << "  Stress_xy max error: " << hier_cc_data_ops.maxNorm(sig_xy_err_idx, wgt_cc_idx) << "\n";

    // Now draw the solutions
    time_integrator->setupPlotData();
    hier_math_ops.interp(u_draw_err_idx, u_draw_err_var, u_err_idx, u_err_var, nullptr, loop_time, false);
    visit_data_writer->writePlotData(hierarchy, iteration_num, loop_time);

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
        level->deallocatePatchData(u_err_idx);
        level->deallocatePatchData(sig_xx_err_idx);
        level->deallocatePatchData(sig_yy_err_idx);
        level->deallocatePatchData(sig_xy_err_idx);
        level->deallocatePatchData(u_draw_err_idx);
        level->deallocatePatchData(sig_xx_exact_idx);
        level->deallocatePatchData(sig_yy_exact_idx);
        level->deallocatePatchData(sig_xy_exact_idx);
        level->deallocatePatchData(u_draw_exact_idx);
        if (deallocate_sig_scr_after) level->deallocatePatchData(sig_scr_idx);
    }
}

void
computeErrorAndPlotLagrangian(const std::unique_ptr<ExodusII_IO>& lower_io,
                              const std::string& lower_filename,
                              EquationSystems* lower_eq_sys,
                              const std::unique_ptr<ExodusII_IO>& upper_io,
                              const std::string& upper_filename,
                              EquationSystems* upper_eq_sys,
                              Pointer<ExtraStressExact> sig_in_fcn,
                              Pointer<ExtraStressExact> sig_out_fcn,
                              Pointer<VelocityInit> vel_fcn,
                              const double loop_time,
                              const int iter_num,
                              const int which_sys)
{
    std::vector<EquationSystems*> eq_systems;
    switch (which_sys)
    {
    case 0:
        eq_systems.push_back(lower_eq_sys);
        break;
    case 1:
        eq_systems.push_back(upper_eq_sys);
        break;
    case 2:
        eq_systems.push_back(lower_eq_sys);
        eq_systems.push_back(upper_eq_sys);
        break;
    default:
        break;
    }

    double sxx_in_max_err = 0.0, syy_in_max_err = 0.0, sxy_in_max_err = 0.0;
    double sxx_out_max_err = 0.0, syy_out_max_err = 0.0, sxy_out_max_err = 0.0;
    double vel_max_err = 0.0;
    for (auto& eq_sys : eq_systems)
    {
        const MeshBase& mesh = eq_sys->get_mesh();
        // Need to copy data to the component versions of the stress
        System& sig_in_sys = eq_sys->get_system(CFIIMethod::EXTRA_STRESS_IN_SYSTEM_NAME);
        System& sig_out_sys = eq_sys->get_system(CFIIMethod::EXTRA_STRESS_OUT_SYSTEM_NAME);
        const DofMap& sig_in_dof_map = sig_in_sys.get_dof_map();
        const DofMap& sig_out_dof_map = sig_out_sys.get_dof_map();
        NumericVector<double>* sig_in_vec = sig_in_sys.current_local_solution.get();
        NumericVector<double>* sig_out_vec = sig_out_sys.current_local_solution.get();
        System& sig_in_err_sys = eq_sys->get_system(SIG_IN_ERR_SYS_NAME);
        const DofMap& sig_in_err_dof_map = sig_in_err_sys.get_dof_map();
        NumericVector<double>* sig_in_err_vec = sig_in_err_sys.solution.get();
        System& sig_out_err_sys = eq_sys->get_system(SIG_OUT_ERR_SYS_NAME);
        const DofMap& sig_out_err_dof_map = sig_out_err_sys.get_dof_map();
        NumericVector<double>* sig_out_err_vec = sig_out_err_sys.solution.get();

        sig_in_err_vec->zero();
        sig_in_err_vec->add(*sig_in_vec);
        sig_in_err_vec->close();
        sig_in_err_sys.update();
        sig_out_err_vec->zero();
        sig_out_err_vec->add(*sig_out_vec);
        sig_out_err_vec->close();
        sig_out_err_sys.update();

        int sig_in_err_sys_num, sig_out_err_sys_num, vel_sys_num;
        for (int sys_num = 0; sys_num < eq_sys->n_systems(); ++sys_num)
        {
            const std::string& sys_name = eq_sys->get_system(sys_num).name();
            if (sys_name == CFIIMethod::VELOCITY_SYSTEM_NAME) vel_sys_num = sys_num;
            if (sys_name == SIG_IN_ERR_SYS_NAME) sig_in_err_sys_num = sys_num;
            if (sys_name == SIG_OUT_ERR_SYS_NAME) sig_out_err_sys_num = sys_num;
        }
        ExactSolution err_est(*eq_sys);
        err_est.set_excluded_subdomains({ 2 });
        err_est.attach_exact_value(vel_sys_num, vel_fcn.getPointer());
        err_est.attach_exact_value(sig_in_err_sys_num, sig_in_fcn.getPointer());
        err_est.attach_exact_value(sig_out_err_sys_num, sig_out_fcn.getPointer());

        err_est.compute_error(CFIIMethod::VELOCITY_SYSTEM_NAME, "U_0");
        vel_max_err = std::max(vel_max_err, err_est.l_inf_error(CFIIMethod::VELOCITY_SYSTEM_NAME, "U_0"));
        err_est.compute_error(CFIIMethod::VELOCITY_SYSTEM_NAME, "U_1");
        vel_max_err = std::max(vel_max_err, err_est.l_inf_error(CFIIMethod::VELOCITY_SYSTEM_NAME, "U_1"));

        err_est.compute_error(SIG_IN_ERR_SYS_NAME, "sig_err_in_xx");
        sxx_in_max_err = std::max(sxx_in_max_err, err_est.l_inf_error(SIG_IN_ERR_SYS_NAME, "sig_err_in_xx"));
        err_est.compute_error(SIG_IN_ERR_SYS_NAME, "sig_err_in_yy");
        syy_in_max_err = std::max(syy_in_max_err, err_est.l_inf_error(SIG_IN_ERR_SYS_NAME, "sig_err_in_yy"));
        err_est.compute_error(SIG_IN_ERR_SYS_NAME, "sig_err_in_xy");
        sxy_in_max_err = std::max(sxy_in_max_err, err_est.l_inf_error(SIG_IN_ERR_SYS_NAME, "sig_err_in_xy"));

        err_est.compute_error(SIG_OUT_ERR_SYS_NAME, "sig_err_out_xx");
        sxx_out_max_err = std::max(sxx_out_max_err, err_est.l_inf_error(SIG_OUT_ERR_SYS_NAME, "sig_err_out_xx"));
        err_est.compute_error(SIG_OUT_ERR_SYS_NAME, "sig_err_out_yy");
        syy_out_max_err = std::max(syy_out_max_err, err_est.l_inf_error(SIG_OUT_ERR_SYS_NAME, "sig_err_out_yy"));
        err_est.compute_error(SIG_OUT_ERR_SYS_NAME, "sig_err_out_xy");
        sxy_out_max_err = std::max(sxy_out_max_err, err_est.l_inf_error(SIG_OUT_ERR_SYS_NAME, "sig_err_out_xy"));
    }

    pout << "\nErrors on structure at time " << loop_time << "\n";
    pout << "   U       " << vel_max_err << "\n";
    pout << "   sxx in  " << sxx_in_max_err << "\n";
    pout << "   syy in  " << syy_in_max_err << "\n";
    pout << "   sxy in  " << sxy_in_max_err << "\n";
    pout << "   sxx out " << sxx_out_max_err << "\n";
    pout << "   syy out " << syy_out_max_err << "\n";
    pout << "   sxy out " << sxy_out_max_err << "\n";

    lower_io->write_timestep(lower_filename, *lower_eq_sys, iter_num, loop_time);
    upper_io->write_timestep(upper_filename, *upper_eq_sys, iter_num, loop_time);
}

void
setInsideOfChannel(const int ls_idx,
                   Pointer<PatchHierarchy<NDIM> > hierarchy,
                   const double y_low,
                   const double y_up,
                   const double theta,
                   const double data_time)
{
    // allocate and set level set data.
    for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
        if (!level->checkAllocated(ls_idx)) level->allocatePatchData(ls_idx, data_time);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
            const double* const dx = pgeom->getDx();
            const double* const xlow = pgeom->getXLower();
            const hier::Index<NDIM>& idx_low = patch->getBox().lower();
            Pointer<CellData<NDIM, double> > ls_data = patch->getPatchData(ls_idx);
            for (CellIterator<NDIM> ci(patch->getBox()); ci; ci++)
            {
                const CellIndex<NDIM>& idx = ci();
                VectorNd x;
                for (int d = 0; d < NDIM; ++d)
                    x[d] = xlow[d] + dx[d] * (static_cast<double>(idx(d) - idx_low(d)) + 0.5);
                const double y_ref = x[1] - std::tan(theta) * x[0];
                (*ls_data)(idx) = std::min(y_up - y_ref, y_ref - y_low);
            }
        }
    }
    // Now fill ghost cells
    using ITC = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
    std::vector<ITC> ghost_cell_comp(1);
    ghost_cell_comp[0] =
        ITC(ls_idx, "CONSERVATIVE_LINEAR_REFINE", false, "NONE", "LINEAR", false, nullptr, nullptr, "LINEAR");
    HierarchyGhostCellInterpolation hier_ghost_fill;
    hier_ghost_fill.initializeOperatorState(ghost_cell_comp, hierarchy, 0, hierarchy->getFinestLevelNumber());
    hier_ghost_fill.fillData(data_time);
}
