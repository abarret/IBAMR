// Copyright (c) 2002-2014, Boyce Griffith
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//    * Redistributions of source code must retain the above copyright notice,
//      this list of conditions and the following disclaimer.
//
//    * Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions and the following disclaimer in the
//      documentation and/or other materials provided with the distribution.
//
//    * Neither the name of The University of North Carolina nor the names of
//      its contributors may be used to endorse or promote products derived from
//      this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

// Config files
#include <IBAMR_config.h>
#include <IBTK_config.h>
#include <SAMRAI_config.h>

// Headers for basic PETSc functions
#include <petscsys.h>

// Headers for basic SAMRAI objects
#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <LoadBalancer.h>
#include <StandardTagAndInitialize.h>

// Headers for basic libMesh objects
#include "libmesh/mesh_modification.h"
#include "libmesh/mesh_refinement.h"
#include "libmesh/mesh_tools.h"
#include <libmesh/boundary_info.h>
#include <libmesh/boundary_mesh.h>
#include <libmesh/equation_systems.h>
#include <libmesh/exodusII_io.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/mesh_triangle_interface.h>

#include "libmesh/face_tri3_subdivision.h"
#include "libmesh/mesh_subdivision_support.h"

// Headers for application-specific algorithm/data structure objects
#include "ibamr/AdvDiffComplexFluidConvectiveOperator.h"
#include "ibamr/ComplexFluidForcing.h"
#include <boost/multi_array.hpp>
#include <ibamr/IBExplicitHierarchyIntegrator.h>
#include <ibamr/IBFEMethod.h>
#include <ibamr/IBFESurfaceMethod.h>
#include <ibamr/INSCollocatedHierarchyIntegrator.h>
#include <ibamr/INSStaggeredHierarchyIntegrator.h>
#include <ibtk/AppInitializer.h>
#include <ibtk/LEInteractor.h>
#include <ibtk/libmesh_utilities.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

// Set up application namespace declarations
#include <ibamr/app_namespaces.h>

// Elasticity model data.
namespace ModelData
{
// Tether (penalty) force functions.
static double kappa_s = 1.0e6;
static double eta_s = 0.0;
static double rho = 1.0;
static double DX = -1.0;
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
    if (ERROR_ON_MOVE && ((d[0] > (0.25 * DX)) || (d[1] > (0.25 * DX))))
    {
        TBOX_ERROR("Structure has moved too much.\n");
    }
    return;
} // tether_force_function
} // namespace ModelData
using namespace ModelData;

// Function prototypes
void postprocess_data(Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
                      Pointer<INSHierarchyIntegrator> navier_stokes_integrator,
                      Pointer<AdvDiffSemiImplicitHierarchyIntegrator> adv_diff_integrator,
                      Pointer<ComplexFluidForcing> complex_fluid,
                      Mesh& mesh,
                      EquationSystems* equation_systems,
                      const int iteration_num,
                      const double loop_time,
                      const string& data_dump_dirname);
void
setInsideOfCylinder(const int d_idx, Pointer<PatchHierarchy<NDIM> > patch_hierarchy, const double time, const double R);

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
void
run_example(int argc, char* argv[])
{
    // Initialize libMesh, PETSc, MPI, and SAMRAI.
    LibMeshInit init(argc, argv);
    SAMRAI_MPI::setCommunicator(PETSC_COMM_WORLD);
    SAMRAI_MPI::setCallAbortInSerialInsteadOfExit();
    SAMRAIManager::setMaxNumberPatchDataEntries(512);
    SAMRAIManager::startup();

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
        const bool uses_exodus = dump_viz_data && !app_initializer->getExodusIIFilename().empty();
        const string exodus_filename = app_initializer->getExodusIIFilename();

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
        Mesh upper_mesh(init.comm(), NDIM - 1), lower_mesh(init.comm(), NDIM - 1);
        const double dx = input_db->getDouble("DX");
        const double mfac = input_db->getDouble("MFAC");
        const double ds = mfac * dx;
        const double L = input_db->getDouble("L");
        std::cout << dx << " " << mfac << "\n";
        std::cout << ds << "\n";
        string elem_type = input_db->getString("ELEM_TYPE");
        const int N = static_cast<int>(L / ds);
        if (NDIM == 2)
        {
            MeshTools::Generation::build_line(upper_mesh, N, 0.0, L, Utility::string_to_enum<ElemType>(elem_type));
            MeshTools::Modification::translate(upper_mesh, 0.0, 2.0, 0.0);
            MeshTools::Generation::build_line(lower_mesh, N, 0.0, L, Utility::string_to_enum<ElemType>(elem_type));
            MeshTools::Modification::translate(lower_mesh, 0.0, -2.0, 0.0);
        }

        upper_mesh.prepare_for_use();
        upper_mesh.print_info();
        lower_mesh.prepare_for_use();
        lower_mesh.print_info();

        std::vector<Mesh*> meshes(2);
        meshes[0] = &upper_mesh;
        meshes[1] = &lower_mesh;

        // MeshRefinement mesh_refinement (mesh);
        // mesh_refinement.uniformly_refine (3);
        // MeshTools::Modification::flatten (mesh);

        kappa_s = input_db->getDouble("KAPPA_S");
        eta_s = input_db->getDouble("ETA_S");
        rho = input_db->getDouble("RHO");
        bool ignore_center = input_db->getBoolWithDefault("IGNORE_CENTER", false);
        ERROR_ON_MOVE = input_db->getBoolWithDefault("ERROR_ON_MOVE", false);
        DX = input_db->getDouble("DX");

        // MeshTools::Subdivision::prepare_subdivision_mesh (mesh, false);

        // Print information about the subdivision mesh to the screen.
        // mesh.print_info();

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
        Pointer<AdvDiffSemiImplicitHierarchyIntegrator> adv_diff_integrator;
        adv_diff_integrator = new AdvDiffSemiImplicitHierarchyIntegrator(
            "AdvDiffSemiImplicitHierarchyIntegrator",
            app_initializer->getComponentDatabase("AdvDiffSemiImplicitHierarchyIntegrator"));
        navier_stokes_integrator->registerAdvDiffHierarchyIntegrator(adv_diff_integrator);
        Pointer<IBFESurfaceMethod> ib_method_ops =
            new IBFESurfaceMethod("IBFEMethod",
                                  app_initializer->getComponentDatabase("IBFEMethod"),
                                  meshes,
                                  app_initializer->getComponentDatabase("GriddingAlgorithm")->getInteger("max_levels"));
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

        // Configure the IBFE solver.
        ib_method_ops->initializeFEEquationSystems();
        std::vector<int> vars(NDIM);
        for (unsigned int d = 0; d < NDIM; ++d) vars[d] = d;
        vector<SystemData> sys_data(1, SystemData(IBFEMethod::VELOCITY_SYSTEM_NAME, vars));
        IBFESurfaceMethod::LagSurfaceForceFcnData body_fcn_data(tether_force_function, sys_data);
        ib_method_ops->registerLagSurfaceForceFunction(body_fcn_data, 0);
        ib_method_ops->registerLagSurfaceForceFunction(body_fcn_data, 1);
        EquationSystems* equation_systems_u = ib_method_ops->getFEDataManager(0)->getEquationSystems();
        EquationSystems* equation_systems_l = ib_method_ops->getFEDataManager(1)->getEquationSystems();

        // Create Eulerian initial condition specification objects.
        if (input_db->keyExists("VelocityInitialConditions"))
        {
            Pointer<CartGridFunction> u_init = new muParserCartGridFunction(
                "u_init", app_initializer->getComponentDatabase("VelocityInitialConditions"), grid_geometry);
            navier_stokes_integrator->registerVelocityInitialConditions(u_init);
        }

        if (input_db->keyExists("PressureInitialConditions"))
        {
            Pointer<CartGridFunction> p_init = new muParserCartGridFunction(
                "p_init", app_initializer->getComponentDatabase("PressureInitialConditions"), grid_geometry);
            navier_stokes_integrator->registerPressureInitialConditions(p_init);
        }

        // Create Eulerian boundary condition specification objects (when necessary).
        const IntVector<NDIM>& periodic_shift = grid_geometry->getPeriodicShift();
        vector<RobinBcCoefStrategy<NDIM>*> u_bc_coefs(NDIM);
        if (periodic_shift.min() > 0)
        {
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                u_bc_coefs[d] = NULL;
            }
        }
        else
        {
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                ostringstream bc_coefs_name_stream;
                bc_coefs_name_stream << "u_bc_coefs_" << d;
                const string bc_coefs_name = bc_coefs_name_stream.str();

                ostringstream bc_coefs_db_name_stream;
                bc_coefs_db_name_stream << "VelocityBcCoefs_" << d;
                const string bc_coefs_db_name = bc_coefs_db_name_stream.str();

                u_bc_coefs[d] = new muParserRobinBcCoefs(
                    bc_coefs_name, app_initializer->getComponentDatabase(bc_coefs_db_name), grid_geometry);
            }
            navier_stokes_integrator->registerPhysicalBoundaryConditions(u_bc_coefs);
        }

        // Set up visualization plot file writers.
        Pointer<VisItDataWriter<NDIM> > visit_data_writer = app_initializer->getVisItDataWriter();
        if (uses_visit)
        {
            time_integrator->registerVisItDataWriter(visit_data_writer);
        }

        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();

        // Create Eulerian body force function specification objects.
        Pointer<ComplexFluidForcing> complex_fluid;
        if (input_db->keyExists("ForcingFunction"))
        {
            Pointer<CartGridFunction> f_fcn = new muParserCartGridFunction(
                "f_fcn", app_initializer->getComponentDatabase("ForcingFunction"), grid_geometry);
            time_integrator->registerBodyForceFunction(f_fcn);
        }
        if (input_db->keyExists("ComplexFluid"))
        {
            complex_fluid = new ComplexFluidForcing("ComplexFluidForcing",
                                                    app_initializer->getComponentDatabase("ComplexFluid"),
                                                    navier_stokes_integrator,
                                                    grid_geometry,
                                                    adv_diff_integrator,
                                                    visit_data_writer);
            time_integrator->registerBodyForceFunction(complex_fluid);
        }

        libMesh::UniquePtr<ExodusII_IO> exodus_io_l(uses_exodus ? new ExodusII_IO(lower_mesh) : NULL);
        libMesh::UniquePtr<ExodusII_IO> exodus_io_u(uses_exodus ? new ExodusII_IO(upper_mesh) : NULL);
        std::string exodus_u_fn = app_initializer->getExodusIIFilename("UpperMesh");
        std::string exodus_l_fn = app_initializer->getExodusIIFilename("LowerMesh");

        Pointer<muParserCartGridFunction> s_exact = new muParserCartGridFunction("S_exact", app_initializer->getComponentDatabase("S_exact"), grid_geometry);
        Pointer<muParserCartGridFunction> u_exact = new muParserCartGridFunction("U_exact", app_initializer->getComponentDatabase("U_exact"), grid_geometry);

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
            if (uses_visit)
            {
                time_integrator->setupPlotData();
                visit_data_writer->writePlotData(patch_hierarchy, iteration_num, loop_time);
            }
            if (uses_exodus)
            {
                exodus_io_u->write_timestep(
                    exodus_u_fn, *equation_systems_u, iteration_num / viz_dump_interval + 1, loop_time);
                exodus_io_l->write_timestep(
                    exodus_l_fn, *equation_systems_l, iteration_num / viz_dump_interval + 1, loop_time);
            }
        }

        // Main time step loop.
        double loop_time_end = time_integrator->getEndTime();
        double dt = 0.0;
        while (!MathUtilities<double>::equalEps(loop_time, loop_time_end) && time_integrator->stepsRemaining())
        {
            iteration_num = time_integrator->getIntegratorStep();
            loop_time = time_integrator->getIntegratorTime();

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
                    time_integrator->setupPlotData();
                    visit_data_writer->writePlotData(patch_hierarchy, iteration_num, loop_time);
                }
                if (uses_exodus)
                {
                    exodus_io_u->write_timestep(
                        exodus_u_fn, *equation_systems_u, iteration_num / viz_dump_interval + 1, loop_time);
                    exodus_io_l->write_timestep(
                        exodus_l_fn, *equation_systems_l, iteration_num / viz_dump_interval + 1, loop_time);
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
            if (dump_viz_data && (iteration_num % postproc_data_dump_interval == 0 || last_step))
            {
                VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
                const Pointer<hier::Variable<NDIM> > p_var = time_integrator->getPressureVariable();
                const Pointer<VariableContext> p_ctx = time_integrator->getCurrentContext();
                const int p_idx = var_db->mapVariableAndContextToIndex(p_var, p_ctx);

                /*postprocess_data(patch_hierarchy,
                                 navier_stokes_integrator,
                                 adv_diff_integrator,
                                 complex_fluid,
                                 upper_mesh,
                                 equation_systems,
                                 iteration_num,
                                 loop_time,
                                 postproc_data_dump_dirname);*/
            }
        }

        Pointer<CellVariable<NDIM, double> > s_var, sxx_var, syy_var, sxy_var;
        int s_idx, s_cloned_idx, sxx_idx, syy_idx, sxy_idx;

        if (complex_fluid)
        {
            s_var = complex_fluid->getVariable();
            sxx_var = new CellVariable<NDIM, double>("Sxx");
            syy_var = new CellVariable<NDIM, double>("Syy");
            sxy_var = new CellVariable<NDIM, double>("Sxy");
            const Pointer<VariableContext> s_ctx = adv_diff_integrator->getCurrentContext();

            s_idx = var_db->mapVariableAndContextToIndex(s_var, s_ctx);
            s_cloned_idx = var_db->registerClonedPatchDataIndex(s_var, s_idx);
            sxx_idx = var_db->registerVariableAndContext(sxx_var, s_ctx);
            syy_idx = var_db->registerVariableAndContext(syy_var, s_ctx);
            sxy_idx = var_db->registerVariableAndContext(sxy_var, s_ctx);
        }

        const Pointer<SAMRAI::hier::Variable<NDIM> > u_var = navier_stokes_integrator->getVelocityVariable();
        const Pointer<VariableContext> u_ctx = navier_stokes_integrator->getCurrentContext();

        const int u_idx = var_db->mapVariableAndContextToIndex(u_var, u_ctx);
        const int u_cloned_idx = var_db->registerClonedPatchDataIndex(u_var, u_idx);

        const int coarsest_ln = 0;
        const int finest_ln = patch_hierarchy->getFinestLevelNumber();
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
            if (complex_fluid)
            {
                level->allocatePatchData(sxx_idx, loop_time);
                level->allocatePatchData(syy_idx, loop_time);
                level->allocatePatchData(sxy_idx, loop_time);
                level->allocatePatchData(s_cloned_idx, loop_time);
            }
            level->allocatePatchData(u_cloned_idx, loop_time);
        }

        if(complex_fluid) s_exact->setDataOnPatchHierarchy(s_cloned_idx, s_var, patch_hierarchy, loop_time);
        u_exact->setDataOnPatchHierarchy(u_cloned_idx, u_var, patch_hierarchy, loop_time);

        HierarchyMathOps hier_math_ops("HierarchyMathOps", patch_hierarchy);
        hier_math_ops.setPatchHierarchy(patch_hierarchy);
        hier_math_ops.resetLevels(coarsest_ln, finest_ln);
        const int wgt_cc_idx = hier_math_ops.getCellWeightPatchDescriptorIndex();
        const int wgt_sc_idx = hier_math_ops.getSideWeightPatchDescriptorIndex();
        Pointer<CellVariable<NDIM, int> > ind_var = new CellVariable<NDIM, int>("Indicator");
        const int ind_idx = var_db->registerVariableAndContext(ind_var, u_ctx);
        visit_data_writer->registerPlotQuantity("Indicator", "SCALAR", ind_idx);
        if (ignore_center)
        {
            for (int ln = 0; ln <= finest_ln; ++ln)
            {
                Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
                level->allocatePatchData(ind_idx, loop_time);
                for (PatchLevel<NDIM>::Iterator p(level); p; p++)
                {
                    Pointer<Patch<NDIM> > patch = level->getPatch(p());
                    Pointer<CellData<NDIM, double> > wgt_cc_data = patch->getPatchData(wgt_cc_idx);
                    Pointer<CellData<NDIM, int> > ind_data = patch->getPatchData(ind_idx);
                    const Box<NDIM>& box = patch->getBox();
                    Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
                    const double* const dx = pgeom->getDx();
                    const double* const x_low = pgeom->getXLower();
                    const SAMRAI::hier::Index<NDIM>& idx_low = box.lower();
                    std::vector<double> x(NDIM);
                    ind_data->fillAll(1);
                    for (CellIterator<NDIM> ci(box); ci; ci++)
                    {
                        const CellIndex<NDIM>& idx = *ci;
                        for (int d = 0; d < NDIM; ++d) x[d] = x_low[d] + dx[d] * (idx(d) - idx_low(d) + 0.5);
                        if (abs(x[1] - 2.0) < (3.0 * dx[1]) || abs(x[1] + 2.0) < (3.0 * dx[1]))
                        {
                            (*wgt_cc_data)(idx) = patch->getPatchData(wgt_cc_idx);
                            (*ind_data)(idx) = 0;
                        }
                    }
                }
            }
        }
        else
        {
            for (int ln = 0; ln <= finest_ln; ++ln)
            {
                Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
                level->allocatePatchData(ind_idx, loop_time);
                for (PatchLevel<NDIM>::Iterator p(level); p; p++)
                {
                    Pointer<Patch<NDIM> > patch = level->getPatch(p());
                    Pointer<CellData<NDIM, int> > ind_data = patch->getPatchData(ind_idx);
                    ind_data->fillAll(1);
                }
            }
        }

        HierarchyCellDataOpsReal<NDIM, double> hier_cc_data_ops(patch_hierarchy, coarsest_ln, finest_ln);
        HierarchySideDataOpsReal<NDIM, double> hier_sc_data_ops(patch_hierarchy, coarsest_ln, finest_ln);
        hier_sc_data_ops.subtract(u_idx, u_idx, u_cloned_idx);
        if (complex_fluid) hier_cc_data_ops.subtract(s_idx, s_idx, s_cloned_idx);

        if (complex_fluid)
        {
            for (int ln = 0; ln <= finest_ln; ++ln)
            {
                Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
                for (PatchLevel<NDIM>::Iterator p(level); p; p++)
                {
                    Pointer<Patch<NDIM> > patch = level->getPatch(p());
                    Pointer<CellData<NDIM, double> > s_data = patch->getPatchData(s_idx);
                    Pointer<CellData<NDIM, double> > sxx_data = patch->getPatchData(sxx_idx);
                    Pointer<CellData<NDIM, double> > syy_data = patch->getPatchData(syy_idx);
                    Pointer<CellData<NDIM, double> > sxy_data = patch->getPatchData(sxy_idx);

                    sxx_data->copyDepth(0, *s_data, 0);
                    syy_data->copyDepth(0, *s_data, 1);
                    sxy_data->copyDepth(0, *s_data, 2);
                }
            }
        }

        pout << "Error in u at time " << loop_time << ":\n"
             << "  L1-norm:  " << std::setprecision(10) << hier_sc_data_ops.L1Norm(u_idx, wgt_sc_idx) << "\n"
             << "  L2-norm:  " << std::setprecision(10) << hier_sc_data_ops.L2Norm(u_idx, wgt_sc_idx) << "\n"
             << "  max-norm: " << std::setprecision(10) << hier_sc_data_ops.maxNorm(u_idx, wgt_sc_idx) << "\n";

        if(complex_fluid)
        {
            pout << "Error in sxx at time " << loop_time << ":\n"
                 << "  L1-norm:  " << std::setprecision(10) << hier_cc_data_ops.L1Norm(sxx_idx, wgt_cc_idx) << "\n"
                 << "  L2-norm:  " << std::setprecision(10) << hier_cc_data_ops.L2Norm(sxx_idx, wgt_cc_idx) << "\n"
                 << "  max-norm: " << std::setprecision(10) << hier_cc_data_ops.maxNorm(sxx_idx, wgt_cc_idx) << "\n";

            pout << "Error in syy at time " << loop_time << ":\n"
                 << "  L1-norm:  " << std::setprecision(10) << hier_cc_data_ops.L1Norm(syy_idx, wgt_cc_idx) << "\n"
                 << "  L2-norm:  " << std::setprecision(10) << hier_cc_data_ops.L2Norm(syy_idx, wgt_cc_idx) << "\n"
                 << "  max-norm: " << std::setprecision(10) << hier_cc_data_ops.maxNorm(syy_idx, wgt_cc_idx) << "\n";

            pout << "Error in sxy at time " << loop_time << ":\n"
                 << "  L1-norm:  " << std::setprecision(10) << hier_cc_data_ops.L1Norm(sxy_idx, wgt_cc_idx) << "\n"
                 << "  L2-norm:  " << std::setprecision(10) << hier_cc_data_ops.L2Norm(sxy_idx, wgt_cc_idx) << "\n"
                 << "  max-norm: " << std::setprecision(10) << hier_cc_data_ops.maxNorm(sxy_idx, wgt_cc_idx) << "\n";
        }

        time_integrator->setupPlotData();
        if(complex_fluid)
        {
            visit_data_writer->registerPlotQuantity("SXX_Err", "SCALAR", sxx_idx);
            visit_data_writer->registerPlotQuantity("SYY_Err", "SCALAR", syy_idx);
            visit_data_writer->registerPlotQuantity("SXY_Err", "SCALAR", sxy_idx);
        }
        visit_data_writer->writePlotData(patch_hierarchy, iteration_num + 1, loop_time);

        for (int ln = 0; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
            level->deallocatePatchData(u_cloned_idx);
            level->deallocatePatchData(s_cloned_idx);
            if (complex_fluid)
            {
                level->deallocatePatchData(sxx_idx);
                level->deallocatePatchData(syy_idx);
                level->deallocatePatchData(sxy_idx);
            }
            level->deallocatePatchData(ind_idx);
        }


        // Cleanup Eulerian boundary condition specification objects (when
        // necessary).
        for (unsigned int d = 0; d < NDIM; ++d) delete u_bc_coefs[d];

    } // cleanup dynamically allocated objects prior to shutdown

    SAMRAIManager::shutdown();
} // main

void
postprocess_data(Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
                 Pointer<INSHierarchyIntegrator> navier_stokes_integrator,
                 Pointer<AdvDiffSemiImplicitHierarchyIntegrator> adv_diff_integrator,
                 Pointer<ComplexFluidForcing> complex_fluid,
                 Mesh& mesh,
                 EquationSystems* equation_systems,
                 const int iteration_num,
                 const double loop_time,
                 const string& data_dump_dirname)
{
    pout << "Outputting data at time: " << loop_time << "\n And iteration num : " << iteration_num << "\n";
    {
        // Output files
        string file_name = data_dump_dirname + "/hier_data.";
        char temp_buf[128];
        sprintf(temp_buf, "%05d.samrai.%05d", iteration_num, SAMRAI_MPI::getRank());
        file_name += temp_buf;
        Pointer<HDFDatabase> hier_db = new HDFDatabase("hier_db");
        hier_db->create(file_name);
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        ComponentSelector hier_data;
        hier_data.setFlag(var_db->mapVariableAndContextToIndex(complex_fluid->getVariable(),
                                                               adv_diff_integrator->getCurrentContext()));
        pout << "current ctx " << adv_diff_integrator->getCurrentContext()->getName() << "\n";
        hier_data.setFlag(var_db->mapVariableAndContextToIndex(navier_stokes_integrator->getVelocityVariable(),
                                                               navier_stokes_integrator->getCurrentContext()));
        hier_data.setFlag(var_db->mapVariableAndContextToIndex(navier_stokes_integrator->getPressureVariable(),
                                                               navier_stokes_integrator->getCurrentContext()));
        patch_hierarchy->putToDatabase(hier_db->putDatabase("PatchHierarchy"), hier_data);
        hier_db->putDouble("loop_time", loop_time);
        hier_db->putInteger("iteration_num", iteration_num);
        hier_db->close();
    }
    const unsigned int dim = mesh.mesh_dimension();
    double F_integral[NDIM];
    double F_int = 0.0;

    for (unsigned int d = 0; d < NDIM; ++d) F_integral[d] = 0.0;

    System& x_system = equation_systems->get_system(IBFEMethod::COORDS_SYSTEM_NAME);
    System& U_system = equation_systems->get_system(IBFEMethod::VELOCITY_SYSTEM_NAME);

    NumericVector<double>* x_vec = x_system.solution.get();
    NumericVector<double>* x_ghost_vec = x_system.current_local_solution.get();
    x_vec->localize(*x_ghost_vec);
    NumericVector<double>* U_vec = U_system.solution.get();
    NumericVector<double>* U_ghost_vec = U_system.current_local_solution.get();
    U_vec->localize(*U_ghost_vec);
    const DofMap& dof_map = x_system.get_dof_map();
    std::vector<std::vector<unsigned int> > dof_indices(NDIM);

    libMesh::UniquePtr<FEBase> fe(FEBase::build(dim, dof_map.variable_type(0)));
    libMesh::UniquePtr<QBase> qrule = QBase::build(QGAUSS, dim, SEVENTH);
    fe->attach_quadrature_rule(qrule.get());
    const vector<double>& JxW = fe->get_JxW();
    const vector<libMesh::Point>& q_point = fe->get_xyz();
    const vector<vector<double> >& phi = fe->get_phi();
    const vector<vector<VectorValue<double> > >& dphi = fe->get_dphi();

    libMesh::UniquePtr<FEBase> fe_face(FEBase::build(dim, dof_map.variable_type(0)));
    libMesh::UniquePtr<QBase> qrule_face(QBase::build(QGAUSS, dim - 1, SEVENTH));
    fe_face->attach_quadrature_rule(qrule_face.get());
    const vector<double>& JxW_face = fe_face->get_JxW();
    const vector<libMesh::Point>& q_point_face = fe_face->get_xyz();
    const vector<libMesh::Point>& normal_face = fe_face->get_normals();
    const vector<vector<double> >& phi_face = fe_face->get_phi();
    const vector<vector<VectorValue<double> > >& dphi_face = fe_face->get_dphi();

    std::vector<double> U_qp_vec(NDIM);
    std::vector<const std::vector<double>*> var_data(1);
    var_data[0] = &U_qp_vec;
    std::vector<const std::vector<libMesh::VectorValue<double> >*> grad_var_data;
    void* force_fcn_ctx = NULL;

    TensorValue<double> FF, FF_inv_trans;
    boost::multi_array<double, 2> x_node, U_node;
    VectorValue<double> F, N, U, n, x;

    const MeshBase::const_element_iterator el_begin = mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator el_end = mesh.active_local_elements_end();
    for (MeshBase::const_element_iterator el_it = el_begin; el_it != el_end; ++el_it)
    {
        Elem* const elem = *el_it;
        fe->reinit(elem);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            dof_map.dof_indices(elem, dof_indices[d], d);
        }
        get_values_for_interpolation(x_node, *x_ghost_vec, dof_indices);
        get_values_for_interpolation(U_node, *U_ghost_vec, dof_indices);

        const unsigned int n_qp = qrule->n_points();
        for (unsigned int qp = 0; qp < n_qp; ++qp)
        {
            interpolate(x, qp, x_node, phi);
            jacobian(FF, qp, x_node, dphi);
            interpolate(U, qp, U_node, phi);
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                U_qp_vec[d] = U(d);
            }
            tether_force_function(
                F, n, N, FF, x, q_point[qp], elem, 0, var_data, grad_var_data, loop_time, force_fcn_ctx);
            F_int += F(0) * JxW[qp];
            for (int d = 0; d < NDIM; ++d)
            {
                F_integral[d] += F(d) * JxW[qp];
            }
        }

        for (unsigned short int side = 0; side < elem->n_sides(); ++side)
        {
            if (elem->neighbor(side)) continue;
            fe_face->reinit(elem, side);
            const unsigned int n_qp_face = qrule_face->n_points();
            for (unsigned int qp = 0; qp < n_qp_face; ++qp)
            {
                interpolate(x, qp, x_node, phi_face);
                jacobian(FF, qp, x_node, dphi_face);
                interpolate(U, qp, U_node, phi_face);
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    U_qp_vec[d] = U(d);
                }
                N = normal_face[qp];
                tensor_inverse_transpose(FF_inv_trans, FF, NDIM);
                tether_force_function(
                    F, n, N, FF, x, q_point_face[qp], elem, side, var_data, grad_var_data, loop_time, force_fcn_ctx);
                for (int d = 0; d < NDIM; ++d)
                {
                    F_integral[d] += F(d) * JxW_face[qp];
                }
            }
        }
    }
    SAMRAI_MPI::sumReduction(F_integral, NDIM);

    return;
} // postprocess_data
