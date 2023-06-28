// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2022 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

/////////////////////////////// INCLUDES /////////////////////////////////////
#include "ibtk/IBTK_CHKERRQ.h"
#include "ibtk/IBTK_MPI.h"
#include "ibtk/LData.h"
#include "ibtk/LSiloDataWriter.h"

#include "IntVector.h"
#include "PatchHierarchy.h"
#include "tbox/Database.h"
#include "tbox/Pointer.h"
#include "tbox/RestartManager.h"
#include "tbox/Utilities.h"

#include "petscao.h"
#include "petscis.h"
#include "petscistypes.h"
#include "petscsys.h"
#include "petscvec.h"

#include <mpi.h>

#include <algorithm>
#include <cstdio>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "ibtk/namespaces.h" // IWYU pragma: keep

#if defined(IBTK_HAVE_SILO)
#include "silo.h"
#endif

/////////////////////////////// STATIC ///////////////////////////////////////
namespace IBTK
{
namespace
{
// The rank of the root MPI process and the MPI tag number.
static const int SILO_MPI_ROOT = 0;
static const int SILO_MPI_TAG = 0;

// The name of the Silo dumps and database filenames.
static const int SILO_NAME_BUFSIZE = 128;
static const std::string VISIT_DUMPS_FILENAME = "lag_data.visit";
static const std::string SILO_DUMP_DIR_PREFIX = "lag_data.cycle_";
static const std::string SILO_SUMMARY_FILE_PREFIX = "lag_data.cycle_";
static const std::string SILO_SUMMARY_FILE_POSTFIX = ".summary.silo";
static const std::string SILO_PROCESSOR_FILE_PREFIX = "lag_data.proc_";
static const std::string SILO_PROCESSOR_FILE_POSTFIX = ".silo";

// Version of LSiloDataWriter restart file data.
static const int LAG_SILO_DATA_WRITER_VERSION = 1;

#if defined(IBTK_HAVE_SILO)
/*!
 * \brief Build a local mesh database entry corresponding to a cloud of marker
 * points.
 */
void
build_local_marker_cloud(DBfile* dbfile,
                         std::string& dirname,
                         const int nmarks,
                         const double* const X,
                         const int nvars,
                         const std::vector<std::string>& varnames,
                         const std::vector<int>& varstartdepths,
                         const std::vector<int>& varplotdepths,
                         const std::vector<int>& vardepths,
                         const std::vector<const double*>& varvals,
                         const int time_step,
                         const double simulation_time)
{
    std::vector<float> block_X(NDIM * nmarks);
    std::vector<std::vector<float> > block_varvals(nvars);
    for (int v = 0; v < nvars; ++v)
    {
        const int varplotdepth = varplotdepths[v];
        block_varvals[v].resize(varplotdepth * nmarks);
    }

    for (int i = 0; i < nmarks; ++i)
    {
        // Get the coordinate data.
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            block_X[d * nmarks + i] = static_cast<float>(X[NDIM * i + d]);
        }

        // Get the variable data.
        for (int v = 0; v < nvars; ++v)
        {
            const int varstartdepth = varstartdepths[v];
            const int varplotdepth = varplotdepths[v];
            const int vardepth = vardepths[v];
            for (int d = 0; d < varplotdepth; ++d)
            {
                block_varvals[v][d * nmarks + i] = static_cast<float>(varvals[v][vardepth * i + varstartdepth + d]);
            }
        }
    }

    // Set the working directory in the Silo database.
    if (DBSetDir(dbfile, dirname.c_str()) == -1)
    {
        TBOX_ERROR("LSiloDataWriter::build_local_marker_cloud()\n"
                   << "  Could not set directory " << dirname << std::endl);
    }

    // Write out the variables.
    int cycle = time_step;
    auto time = static_cast<float>(simulation_time);
    double dtime = simulation_time;

    static const int MAX_OPTS = 3;
    DBoptlist* optlist = DBMakeOptlist(MAX_OPTS);
    DBAddOption(optlist, DBOPT_CYCLE, &cycle);
    DBAddOption(optlist, DBOPT_TIME, &time);
    DBAddOption(optlist, DBOPT_DTIME, &dtime);

    const char* meshname = "mesh";
    std::vector<float*> coords(NDIM);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        coords[d] = nmarks > 0 ? &block_X[d * nmarks] : nullptr;
    }

    int ndims = NDIM;

    DBPutPointmesh(dbfile, meshname, ndims, &coords[0], nmarks, DB_FLOAT, optlist);

    for (int v = 0; v < nvars; ++v)
    {
        const char* varname = varnames[v].c_str();
        const int varplotdepth = varplotdepths[v];

        std::vector<float*> vars(varplotdepth);
        for (int d = 0; d < varplotdepth; ++d)
        {
            vars[d] = nmarks > 0 ? &block_varvals[v][d * nmarks] : nullptr;
        }

        if (varplotdepth == 1)
        {
            DBPutPointvar1(dbfile, varname, meshname, vars[0], nmarks, DB_FLOAT, optlist);
        }
        else
        {
            DBPutPointvar(dbfile, varname, meshname, varplotdepth, &vars[0], nmarks, DB_FLOAT, optlist);
        }
    }

    DBFreeOptlist(optlist);

    // Reset the working directory in the Silo database.
    if (DBSetDir(dbfile, "..") == -1)
    {
        TBOX_ERROR("LSiloDataWriter::build_local_marker_cloud()\n"
                   << "  Could not return to the base directory from subdirectory " << dirname << std::endl);
    }
    return;
} // build_local_marker_cloud

/*!
 * \brief Build a local mesh database entry corresponding to a quadrilateral
 * curvilinear block.
 */
void
build_local_curv_block(DBfile* dbfile,
                       std::string& dirname,
                       const IntVector<NDIM>& nelem_in,
                       const IntVector<NDIM>& periodic,
                       const double* const X,
                       const int nvars,
                       const std::vector<std::string>& varnames,
                       const std::vector<int>& varstartdepths,
                       const std::vector<int>& varplotdepths,
                       const std::vector<int>& vardepths,
                       const std::vector<const double*>& varvals,
                       const int time_step,
                       const double simulation_time)
{
    // Check for co-dimension 1 or 2 data.
    IntVector<NDIM> nelem, degenerate;
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        if (nelem_in(d) == 1)
        {
            nelem(d) = 2;
            degenerate(d) = 1;
        }
        else
        {
            nelem(d) = nelem_in(d);
            degenerate(d) = 0;
        }
    }

    // Rearrange the data into the format required by Silo.
    const int ntot = 1 * (periodic(0) ? nelem(0) + 1 : nelem(0))
#if (NDIM > 1)
                     * (periodic(1) ? nelem(1) + 1 : nelem(1))
#if (NDIM > 2)
                     * (periodic(2) ? nelem(2) + 1 : nelem(2))
#endif
#endif
        ;

    std::vector<float> block_X(NDIM * ntot);
    std::vector<std::vector<float> > block_varvals(nvars);
    for (int v = 0; v < nvars; ++v)
    {
        const int varplotdepth = varplotdepths[v];
        block_varvals[v].resize(varplotdepth * ntot);
    }

    int offset = 0;
#if (NDIM > 2)
    for (int k = 0; k < nelem(2) + (periodic(2) ? 1 : 0); ++k)
    {
#endif
#if (NDIM > 1)
        for (int j = 0; j < nelem(1) + (periodic(1) ? 1 : 0); ++j)
        {
#endif
            for (int i = 0; i < nelem(0) + (periodic(0) ? 1 : 0); ++i)
            {
                const int idx = +(degenerate(0) ? 0 : (i % nelem(0)))
#if (NDIM > 1)
                                + (degenerate(1) ? 0 : (j % nelem(1)) * nelem(0))
#if (NDIM > 2)
                                + (degenerate(2) ? 0 : (k % nelem(2)) * nelem(1) * nelem(0))
#endif
#endif
                    ;

                // Get the coordinate data.
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    block_X[d * ntot + offset] = static_cast<float>(X[NDIM * idx + d]);
                }

                // Get the variable data.
                for (int v = 0; v < nvars; ++v)
                {
                    const int varstartdepth = varstartdepths[v];
                    const int varplotdepth = varplotdepths[v];
                    const int vardepth = vardepths[v];
                    for (int d = 0; d < varplotdepth; ++d)
                    {
                        block_varvals[v][d * ntot + offset] =
                            static_cast<float>(varvals[v][vardepth * idx + varstartdepth + d]);
                    }
                }

                // Increment the counter.
                ++offset;
            }
#if (NDIM > 1)
        }
#endif
#if (NDIM > 2)
    }
#endif

    // Set the working directory in the Silo database.
    if (DBSetDir(dbfile, dirname.c_str()) == -1)
    {
        TBOX_ERROR("LSiloDataWriter::build_local_curv_block()\n"
                   << "  Could not set directory " << dirname << std::endl);
    }

    // Write out the variables.
    int cycle = time_step;
    auto time = static_cast<float>(simulation_time);
    double dtime = simulation_time;

    static const int MAX_OPTS = 3;
    DBoptlist* optlist = DBMakeOptlist(MAX_OPTS);
    DBAddOption(optlist, DBOPT_CYCLE, &cycle);
    DBAddOption(optlist, DBOPT_TIME, &time);
    DBAddOption(optlist, DBOPT_DTIME, &dtime);

    const char* meshname = "mesh";
    const char* coordnames[3] = { "xcoords", "ycoords", "zcoords" };
    std::vector<float*> coords(NDIM);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        coords[d] = ntot > 0 ? &block_X[d * ntot] : nullptr;
    }

    int ndims = NDIM;
    std::vector<int> dims(NDIM);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        dims[d] = nelem(d) + (periodic(d) ? 1 : 0);
    }

    DBPutQuadmesh(dbfile,
                  meshname,
                  const_cast<char**>(coordnames),
                  &coords[0],
                  &dims[0],
                  ndims,
                  DB_FLOAT,
                  DB_NONCOLLINEAR,
                  optlist);

    for (int v = 0; v < nvars; ++v)
    {
        const char* varname = varnames[v].c_str();
        const int varplotdepth = varplotdepths[v];
        std::vector<std::string> compnames;
        for (int d = 0; d < varplotdepth; ++d)
        {
            compnames.push_back(varnames[v] + "_" + std::to_string(d));
        }
        std::vector<const char*> compnames_ptrs;
        for (int d = 0; d < varplotdepth; ++d)
        {
            compnames_ptrs.push_back(compnames[d].c_str());
        }

        std::vector<float*> vars(varplotdepth);
        for (int d = 0; d < varplotdepth; ++d)
        {
            vars[d] = ntot > 0 ? &block_varvals[v][d * ntot] : nullptr;
        }

        if (varplotdepth == 1)
        {
            DBPutQuadvar1(
                dbfile, varname, meshname, vars[0], &dims[0], ndims, nullptr, 0, DB_FLOAT, DB_NODECENT, optlist);
        }
        else
        {
            DBPutQuadvar(dbfile,
                         varname,
                         meshname,
                         varplotdepth,
                         &compnames_ptrs[0],
                         &vars[0],
                         &dims[0],
                         ndims,
                         nullptr,
                         0,
                         DB_FLOAT,
                         DB_NODECENT,
                         optlist);
        }
    }

    DBFreeOptlist(optlist);

    // Reset the working directory in the Silo database.
    if (DBSetDir(dbfile, "..") == -1)
    {
        TBOX_ERROR("LSiloDataWriter::build_local_curv_block()\n"
                   << "  Could not return to the base directory from subdirectory " << dirname << std::endl);
    }
    return;
} // build_local_curv_block

/*!
 * \brief Build a local mesh database entry corresponding to an unstructured
 * mesh.
 */
void
build_local_ucd_mesh(DBfile* dbfile,
                     std::string& dirname,
                     const std::set<int>& vertices,
                     const std::multimap<int, std::pair<int, int> >& edge_map,
                     const double* const X,
                     const int nvars,
                     const std::vector<std::string>& varnames,
                     const std::vector<int>& varstartdepths,
                     const std::vector<int>& varplotdepths,
                     const std::vector<int>& vardepths,
                     const std::vector<const double*>& varvals,
                     const int time_step,
                     const double simulation_time)
{
    // Rearrange the data into the format required by Silo.
    const int ntot = static_cast<int>(vertices.size());

    std::vector<float> block_X(NDIM * ntot);
    std::vector<std::vector<float> > block_varvals(nvars);
    for (int v = 0; v < nvars; ++v)
    {
        const int varplotdepth = varplotdepths[v];
        block_varvals[v].resize(varplotdepth * ntot);
    }

    int offset = 0;
    std::map<int, int> local_vertex_map;
    for (int idx : vertices)
    {
        local_vertex_map[idx] = offset;

        // Get the coordinate data.
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            block_X[d * ntot + offset] = static_cast<float>(X[NDIM * offset + d]);
        }

        // Get the variable data.
        for (int v = 0; v < nvars; ++v)
        {
            const int varstartdepth = varstartdepths[v];
            const int varplotdepth = varplotdepths[v];
            const int vardepth = vardepths[v];
            for (int d = 0; d < varplotdepth; ++d)
            {
                block_varvals[v][d * ntot + offset] =
                    static_cast<float>(varvals[v][vardepth * offset + varstartdepth + d]);
            }
        }

        // Increment the counter.
        ++offset;
    }

    // Prune duplicate edges.
    std::set<std::pair<int, int> > local_edge_set;
    for (const auto& edge_pair : edge_map)
    {
        std::pair<int, int> e = edge_pair.second;
#if !defined(NDEBUG)
        TBOX_ASSERT(vertices.count(e.first) == 1);
        TBOX_ASSERT(vertices.count(e.second) == 1);
#endif
        if (e.first > e.second)
        {
            std::swap<int>(e.first, e.second);
        }
        local_edge_set.insert(e);
    }

    // Create an edge map corresponding to the pruned edge list.
    std::multimap<int, int> local_edge_map;
    for (const auto& edge_pair : local_edge_set)
    {
        const int e1 = edge_pair.first;
        const int e2 = edge_pair.second;
        local_edge_map.insert(std::make_pair(local_vertex_map[e1], local_vertex_map[e2]));
    }

    // Set the working directory in the Silo database.
    if (DBSetDir(dbfile, dirname.c_str()) == -1)
    {
        TBOX_ERROR("LSiloDataWriter::build_local_ucd_mesh()\n"
                   << "  Could not set directory " << dirname << std::endl);
    }

    // Node coordinates.
    int ndims = NDIM;

    int cycle = time_step;
    auto time = static_cast<float>(simulation_time);
    double dtime = simulation_time;

    static const int MAX_OPTS = 3;
    DBoptlist* optlist = DBMakeOptlist(MAX_OPTS);
    DBAddOption(optlist, DBOPT_CYCLE, &cycle);
    DBAddOption(optlist, DBOPT_TIME, &time);
    DBAddOption(optlist, DBOPT_DTIME, &dtime);

    const char* meshname = "mesh";
    const char* coordnames[3] = { "xcoords", "ycoords", "zcoords" };
    std::vector<float*> coords(NDIM);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        coords[d] = ntot > 0 ? &block_X[d * ntot] : nullptr;
    }
    const int nnodes = ntot;

    // Connectivity.
    std::vector<int> nodelist;
    nodelist.reserve(2 * local_edge_map.size());

    for (const auto& edge_pair : local_edge_map)
    {
        nodelist.push_back(edge_pair.first);
        nodelist.push_back(edge_pair.second);
    }
    int lnodelist = static_cast<int>(nodelist.size());
    int nshapetypes = 1;
    int shapecnt[] = { static_cast<int>(local_edge_map.size()) };
    int shapesize[] = { 2 };
    int shapetype[] = { DB_ZONETYPE_BEAM };
    int nzones = static_cast<int>(local_edge_map.size());

    // Write out connectivity information.
    const int origin = 0;
    const int lo_offset = 0;
    const int hi_offset = 0;

    // Write out connectivity information.
    DBPutZonelist2(dbfile,
                   "zonelist",
                   nzones,
                   ndims,
                   lnodelist > 0 ? &nodelist[0] : nullptr,
                   lnodelist,
                   origin,
                   lo_offset,
                   hi_offset,
                   shapetype,
                   shapesize,
                   shapecnt,
                   nshapetypes,
                   optlist);

    // Write an unstructured mesh.
    DBPutUcdmesh(dbfile,
                 meshname,
                 ndims,
                 const_cast<char**>(coordnames),
                 &coords[0],
                 nnodes,
                 nzones,
                 "zonelist",
                 nullptr,
                 DB_FLOAT,
                 nullptr);

    // Write the variables defined on the unstructured mesh.
    for (int v = 0; v < nvars; ++v)
    {
        const char* varname = varnames[v].c_str();
        const int varplotdepth = varplotdepths[v];
        std::vector<std::string> compnames;
        for (int d = 0; d < varplotdepth; ++d)
        {
            compnames.push_back(varnames[v] + "_" + std::to_string(d));
        }
        std::vector<const char*> compnames_ptrs;
        for (int d = 0; d < varplotdepth; ++d)
        {
            compnames_ptrs.push_back(compnames[d].c_str());
        }

        std::vector<float*> vars(varplotdepth);
        for (int d = 0; d < varplotdepth; ++d)
        {
            vars[d] = ntot > 0 ? &block_varvals[v][d * ntot] : nullptr;
        }

        if (varplotdepth == 1)
        {
            DBPutUcdvar1(dbfile, varname, meshname, vars[0], nnodes, nullptr, 0, DB_FLOAT, DB_NODECENT, optlist);
        }
        else
        {
            DBPutUcdvar(dbfile,
                        varname,
                        meshname,
                        varplotdepth,
                        &compnames_ptrs[0],
                        &vars[0],
                        nnodes,
                        nullptr,
                        0,
                        DB_FLOAT,
                        DB_NODECENT,
                        optlist);
        }
    }

    DBFreeOptlist(optlist);

    // Reset the working directory in the Silo database.
    if (DBSetDir(dbfile, "..") == -1)
    {
        TBOX_ERROR("LSiloDataWriter::build_local_ucd_mesh()\n"
                   << "  Could not return to the base directory from subdirectory " << dirname << std::endl);
    }
    return;
} // build_local_ucd_mesh
#endif // if defined(IBTK_HAVE_SILO)
} // namespace

/////////////////////////////// PUBLIC ///////////////////////////////////////

LSiloDataWriter::LSiloDataWriter(std::string object_name, std::string dump_directory_name, bool register_for_restart)
    : d_object_name(std::move(object_name)),
      d_registered_for_restart(register_for_restart),
      d_dump_directory_name(std::move(dump_directory_name)),
      d_nclouds(d_num_objs, 0),
      d_cloud_names(d_num_objs),
      d_cloud_nmarks(d_num_objs),
      d_cloud_first_lag_idx(d_num_objs),
      d_nblocks(d_num_objs, 0),
      d_block_names(d_num_objs),
      d_block_nelems(d_num_objs),
      d_block_periodic(d_num_objs),
      d_block_first_lag_idx(d_num_objs),
      d_nmbs(d_num_objs, 0),
      d_mb_names(d_num_objs),
      d_mb_nblocks(d_num_objs),
      d_mb_nelems(d_num_objs),
      d_mb_periodic(d_num_objs),
      d_mb_first_lag_idx(d_num_objs),
      d_nucd_meshes(d_num_objs, 0),
      d_ucd_mesh_names(d_num_objs),
      d_ucd_mesh_vertices(d_num_objs),
      d_ucd_mesh_edge_maps(d_num_objs),
      d_coords_data(d_num_objs, Pointer<LData>(nullptr)),
      d_nvars(d_num_objs, 0),
      d_var_names(d_num_objs),
      d_var_start_depths(d_num_objs),
      d_var_plot_depths(d_num_objs),
      d_var_depths(d_num_objs),
      d_var_data(d_num_objs),
      d_ao(d_num_objs),
      d_build_vec_scatters(d_num_objs),
      d_src_vec(d_num_objs),
      d_dst_vec(d_num_objs),
      d_vec_scatter(d_num_objs)
{
#if defined(IBTK_HAVE_SILO)
// intentionally blank
#else
    TBOX_WARNING("LSiloDataWriter::LSiloDataWriter(): SILO is not installed; cannot write data." << std::endl);
#endif
    if (d_registered_for_restart)
    {
        RestartManager::getManager()->registerRestartItem(d_object_name, this);
    }

    // Initialize object with data read from the restart database.
    const bool from_restart = RestartManager::getManager()->isFromRestart();
    if (from_restart)
    {
        getFromRestart();
    }
    return;
} // LSiloDataWriter

LSiloDataWriter::~LSiloDataWriter()
{
    if (d_registered_for_restart)
    {
        RestartManager::getManager()->unregisterRestartItem(d_object_name);
    }

    // Destroy any remaining PETSc objects.
    int ierr;
    for (int obj_num = 0; obj_num < d_num_objs; ++obj_num)
    {
        for (auto& vec : d_dst_vec[obj_num])
        {
            Vec& v = vec.second;
            if (v)
            {
                ierr = VecDestroy(&v);
                IBTK_CHKERRQ(ierr);
            }
        }
        for (auto& vec : d_vec_scatter[obj_num])
        {
            VecScatter& vs = vec.second;
            if (vs)
            {
                ierr = VecScatterDestroy(&vs);
                IBTK_CHKERRQ(ierr);
            }
        }
    }
    return;
} // ~LSiloDataWriter

void
LSiloDataWriter::resetNumObjs(const int num_objs)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(num_objs > 0);
#endif
    d_num_objs = std::max(num_objs, d_num_objs);
    // Resize some arrays.
    d_nclouds.resize(d_num_objs, 0);
    d_cloud_names.resize(d_num_objs);
    d_cloud_nmarks.resize(d_num_objs);
    d_cloud_first_lag_idx.resize(d_num_objs);

    d_nblocks.resize(d_num_objs, 0);
    d_block_names.resize(d_num_objs);
    d_block_nelems.resize(d_num_objs);
    d_block_periodic.resize(d_num_objs);
    d_block_first_lag_idx.resize(d_num_objs);

    d_nmbs.resize(d_num_objs, 0);
    d_mb_nblocks.resize(d_num_objs);
    d_mb_names.resize(d_num_objs);
    d_mb_nelems.resize(d_num_objs);
    d_mb_periodic.resize(d_num_objs);
    d_mb_first_lag_idx.resize(d_num_objs);

    d_nucd_meshes.resize(d_num_objs, 0);
    d_ucd_mesh_names.resize(d_num_objs);
    d_ucd_mesh_vertices.resize(d_num_objs);
    d_ucd_mesh_edge_maps.resize(d_num_objs);

    d_coords_data.resize(d_num_objs, nullptr);
    d_nvars.resize(d_num_objs, 0);
    d_var_names.resize(d_num_objs);
    d_var_start_depths.resize(d_num_objs);
    d_var_plot_depths.resize(d_num_objs);
    d_var_depths.resize(d_num_objs);
    d_var_data.resize(d_num_objs);

    d_ao.resize(d_num_objs);
    d_build_vec_scatters.resize(d_num_objs);
    d_src_vec.resize(d_num_objs);
    d_dst_vec.resize(d_num_objs);
    d_vec_scatter.resize(d_num_objs);
    return;
} // resetNumObjs

void
LSiloDataWriter::registerMarkerCloud(const std::string& name,
                                     const int nmarks,
                                     const int first_lag_idx,
                                     const int obj_num)
{
    if (obj_num >= d_num_objs) resetNumObjs(obj_num + 1);

#if !defined(NDEBUG)
    TBOX_ASSERT(nmarks > 0);
    TBOX_ASSERT(obj_num < d_num_objs);
#endif

    // Check to see if the cloud name has already been registered.
    if (find(d_block_names[obj_num].begin(), d_block_names[obj_num].end(), name) != d_block_names[obj_num].end())
    {
        TBOX_ERROR(d_object_name << "::registerMarkerCloud()\n"
                                 << "  marker clouds must have unique names.\n"
                                 << "  a Cartesian block named ``" << name << "'' has already been registered.\n");
    }

    if (find(d_mb_names[obj_num].begin(), d_mb_names[obj_num].end(), name) != d_mb_names[obj_num].end())
    {
        TBOX_ERROR(d_object_name << "::registerMarkerCloud()\n"
                                 << "  marker clouds must have unique names.\n"
                                 << "  a Cartesian multiblock named ``" << name << "'' has already been registered.\n");
    }

    if (find(d_ucd_mesh_names[obj_num].begin(), d_ucd_mesh_names[obj_num].end(), name) !=
        d_ucd_mesh_names[obj_num].end())
    {
        TBOX_ERROR(d_object_name << "::registerMarkerCloud()\n"
                                 << "  marker clouds must have unique names.\n"
                                 << "  an unstructured mesh named ``" << name << "'' has already been registered.\n");
    }

    // Check to see if we are updating a previously registered cloud.
    for (int k = 0; k < d_nclouds[obj_num]; ++k)
    {
        if (d_cloud_names[obj_num][k] == name)
        {
            d_cloud_nmarks[obj_num][k] = nmarks;
            d_cloud_first_lag_idx[obj_num][k] = first_lag_idx;
            return;
        }
    }

    // Record the layout of the marker cloud.
    ++d_nclouds[obj_num];
    d_cloud_names[obj_num].push_back(name);
    d_cloud_nmarks[obj_num].push_back(nmarks);
    d_cloud_first_lag_idx[obj_num].push_back(first_lag_idx);
    return;
} // registerMarkerCloud

void
LSiloDataWriter::registerLogicallyCartesianBlock(const std::string& name,
                                                 const IntVector<NDIM>& nelem,
                                                 const IntVector<NDIM>& periodic,
                                                 const int first_lag_idx,
                                                 const int obj_num)
{
    if (obj_num >= d_num_objs) resetNumObjs(obj_num + 1);

#if !defined(NDEBUG)
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        TBOX_ASSERT(nelem(d) > 0);
        TBOX_ASSERT(periodic(d) == 0 || periodic(d) == 1);
    }
    TBOX_ASSERT(obj_num < d_num_objs);
#endif

    // Check to see if the block name has already been registered.
    if (find(d_cloud_names[obj_num].begin(), d_cloud_names[obj_num].end(), name) != d_cloud_names[obj_num].end())
    {
        TBOX_ERROR(d_object_name << "::registerLogicallyCartesianBlock()\n"
                                 << "  Cartesian blocks must have unique names.\n"
                                 << "  a marker cloud named ``" << name << "'' has already been registered.\n");
    }

    if (find(d_mb_names[obj_num].begin(), d_mb_names[obj_num].end(), name) != d_mb_names[obj_num].end())
    {
        TBOX_ERROR(d_object_name << "::registerLogicallyCartesianBlock()\n"
                                 << "  Cartesian blocks must have unique names.\n"
                                 << "  a Cartesian multiblock named ``" << name << "'' has already been registered.\n");
    }

    if (find(d_ucd_mesh_names[obj_num].begin(), d_ucd_mesh_names[obj_num].end(), name) !=
        d_ucd_mesh_names[obj_num].end())
    {
        TBOX_ERROR(d_object_name << "::registerLogicallyCartesianBlock()\n"
                                 << "  Cartesian blocks must have unique names.\n"
                                 << "  an unstructured mesh named ``" << name << "'' has already been registered.\n");
    }

    // Check to see if we are updating a previously registered block.
    for (int k = 0; k < d_nblocks[obj_num]; ++k)
    {
        if (d_block_names[obj_num][k] == name)
        {
            d_block_nelems[obj_num][k] = nelem;
            d_block_periodic[obj_num][k] = periodic;
            d_block_first_lag_idx[obj_num][k] = first_lag_idx;
            return;
        }
    }

    // Record the layout of the logically Cartesian block.
    ++d_nblocks[obj_num];
    d_block_names[obj_num].push_back(name);
    d_block_nelems[obj_num].push_back(nelem);
    d_block_periodic[obj_num].push_back(periodic);
    d_block_first_lag_idx[obj_num].push_back(first_lag_idx);
    return;
} // registerLogicallyCartesianBlock

void
LSiloDataWriter::registerLogicallyCartesianMultiblock(const std::string& name,
                                                      const std::vector<IntVector<NDIM> >& nelem,
                                                      const std::vector<IntVector<NDIM> >& periodic,
                                                      const std::vector<int>& first_lag_idx,
                                                      const int obj_num)
{
    if (obj_num >= d_num_objs) resetNumObjs(obj_num + 1);

#if !defined(NDEBUG)
    TBOX_ASSERT(periodic.size() == nelem.size());
    TBOX_ASSERT(first_lag_idx.size() == nelem.size());
    size_t sz = nelem.size();
    for (size_t i = 0; i < sz; ++i)
    {
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            TBOX_ASSERT(nelem[i](d) > 0);
            TBOX_ASSERT(periodic[i](d) == 0 || periodic[i](d) == 1);
        }
    }
    TBOX_ASSERT(obj_num < d_num_objs);
#endif

    // Check to see if the multiblock name has already been registered.
    if (find(d_cloud_names[obj_num].begin(), d_cloud_names[obj_num].end(), name) != d_cloud_names[obj_num].end())
    {
        TBOX_ERROR(d_object_name << "::registerLogicallyCartesianMultiblock()\n"
                                 << "  Cartesian multiblocks must have unique names.\n"
                                 << "  a marker cloud named ``" << name << "'' has already been registered.\n");
    }

    if (find(d_block_names[obj_num].begin(), d_block_names[obj_num].end(), name) != d_block_names[obj_num].end())
    {
        TBOX_ERROR(d_object_name << "::registerLogicallyCartesianMultiblock()\n"
                                 << "  Cartesian multiblocks must have unique names.\n"
                                 << "  a Cartesian block named ``" << name << "'' has already been registered.\n");
    }

    if (find(d_ucd_mesh_names[obj_num].begin(), d_ucd_mesh_names[obj_num].end(), name) !=
        d_ucd_mesh_names[obj_num].end())
    {
        TBOX_ERROR(d_object_name << "::registerLogicallyCartesianMultiblock()\n"
                                 << "  Cartesian multiblocks must have unique names.\n"
                                 << "  an unstructured mesh named ``" << name << "'' has already been registered.\n");
    }

    // Check to see if we are updating a previously registered multiblock.
    for (int k = 0; k < d_nmbs[obj_num]; ++k)
    {
        if (d_mb_names[obj_num][k] == name)
        {
            d_mb_nblocks[obj_num][k] = static_cast<int>(nelem.size());
            d_mb_nelems[obj_num][k] = nelem;
            d_mb_periodic[obj_num][k] = periodic;
            d_mb_first_lag_idx[obj_num][k] = first_lag_idx;
            return;
        }
    }

    // Record the layout of the logically Cartesian multiblock.
    ++d_nmbs[obj_num];
    d_mb_names[obj_num].push_back(name);
    d_mb_nblocks[obj_num].push_back(static_cast<int>(nelem.size()));
    d_mb_nelems[obj_num].push_back(nelem);
    d_mb_periodic[obj_num].push_back(periodic);
    d_mb_first_lag_idx[obj_num].push_back(first_lag_idx);
    return;
} // registerLogicallyCartesianMultiblock

void
LSiloDataWriter::registerUnstructuredMesh(const std::string& name,
                                          const std::multimap<int, std::pair<int, int> >& edge_map,
                                          const int obj_num)
{
    if (obj_num >= d_num_objs) resetNumObjs(obj_num + 1);

#if !defined(NDEBUG)
    TBOX_ASSERT(obj_num < d_num_objs);
#endif

    // Check to see if the unstructured mesh name has already been registered.
    if (find(d_cloud_names[obj_num].begin(), d_cloud_names[obj_num].end(), name) != d_cloud_names[obj_num].end())
    {
        TBOX_ERROR(d_object_name << "::registerUnstructuredMesh()\n"
                                 << "  unstructured meshes must have unique names.\n"
                                 << "  a marker cloud named ``" << name << "'' has already been registered.\n");
    }

    if (find(d_block_names[obj_num].begin(), d_block_names[obj_num].end(), name) != d_block_names[obj_num].end())
    {
        TBOX_ERROR(d_object_name << "::registerUnstructuredMesh()\n"
                                 << "  unstructured meshes must have unique names.\n"
                                 << "  a Cartesian block named ``" << name << "'' has already been registered.\n");
    }

    if (find(d_mb_names[obj_num].begin(), d_mb_names[obj_num].end(), name) != d_mb_names[obj_num].end())
    {
        TBOX_ERROR(d_object_name << "::registerUnstructuredMesh()\n"
                                 << "  unstructured meshes must have unique names.\n"
                                 << "  a Cartesian multiblock named ``" << name << "'' has already been registered.\n");
    }

    // Extract the list of vertices from the list of edges.
    std::set<int> vertices;
    for (const auto& edge_pair : edge_map)
    {
        const std::pair<int, int>& e = edge_pair.second;
        vertices.insert(e.first);
        vertices.insert(e.second);
    }

    // Check to see if we are updating a previously registered unstructured
    // mesh.
    for (int k = 0; k < d_nucd_meshes[obj_num]; ++k)
    {
        if (d_ucd_mesh_names[obj_num][k] == name)
        {
            d_ucd_mesh_vertices[obj_num][k] = vertices;
            d_ucd_mesh_edge_maps[obj_num][k] = edge_map;
            return;
        }
    }

    // Record the layout of the unstructured mesh.
    ++d_nucd_meshes[obj_num];
    d_ucd_mesh_names[obj_num].push_back(name);
    d_ucd_mesh_vertices[obj_num].push_back(vertices);
    d_ucd_mesh_edge_maps[obj_num].push_back(edge_map);
    return;
} // registerUnstructuredMesh

void
LSiloDataWriter::registerCrossUnstructuredMesh(const std::string& name,
                                               const std::vector<CrossEdge>& edges,
                                               const std::set<int>& obj_nums)
{
    int max_obj_num = *std::max_element(obj_nums.begin(), obj_nums.end());
    if (max_obj_num >= d_num_objs) resetNumObjs(max_obj_num + 1);

#if ~defined(NDEBUG)
    TBOX_ASSERT(max_obj_num < d_num_objs);
#endif

    // Extract the list of vertices from the list of edges
    std::map<std::pair<int, int>, std::set<int> > vertices;
    for (const auto& edge : edges)
    {
        const std::pair<std::pair<int, int>, int>& pt1 = edge.pt1;
        vertices[pt1.first].insert(pt1.second);
        const std::pair<std::pair<int, int>, int>& pt2 = edge.pt2;
        vertices[pt2.first].insert(pt2.second);
    }

    // Check to see if we are updating a previously registered unstructured mesh.
    for (int k = 0; k < d_nucd_cross_meshes; ++k)
    {
        if (d_ucd_cross_mesh_names[k] == name)
        {
            d_ucd_cross_edges[k] = edges;
            d_ucd_cross_mesh_vertices[k] = vertices;
            d_ucd_cross_mesh_obj_nums[k] = obj_nums;
            return;
        }
    }

    // Record the layout of the mesh.
    ++d_nucd_cross_meshes;
    d_ucd_cross_mesh_names.push_back(name);
    d_ucd_cross_edges.push_back(edges);
    d_ucd_cross_mesh_vertices.push_back(vertices);
    d_ucd_cross_mesh_obj_nums.push_back(obj_nums);
}

void
LSiloDataWriter::registerCoordsData(Pointer<LData> coords_data, const int obj_num)
{
    if (obj_num >= d_num_objs) resetNumObjs(obj_num + 1);
#if !defined(NDEBUG)
    TBOX_ASSERT(coords_data);
    TBOX_ASSERT(coords_data->getDepth() == NDIM);
    TBOX_ASSERT(obj_num < d_num_objs);
#endif
    d_coords_data[obj_num] = coords_data;
    return;
} // registerCoordsData

void
LSiloDataWriter::registerVariableData(const std::string& var_name, Pointer<LData> var_data, const int obj_num)
{
    const int start_depth = 0;
    const int var_depth = var_data->getDepth();
    registerVariableData(var_name, var_data, start_depth, var_depth, obj_num);
    return;
} // registerVariableData

void
LSiloDataWriter::registerVariableData(const std::string& var_name,
                                      Pointer<LData> var_data,
                                      const int start_depth,
                                      const int var_depth,
                                      const int obj_num)
{
    if (obj_num >= d_num_objs) resetNumObjs(obj_num + 1);

#if !defined(NDEBUG)
    TBOX_ASSERT(!var_name.empty());
    TBOX_ASSERT(var_data);
    TBOX_ASSERT(obj_num < d_num_objs);
#endif
    if (find(d_var_names[obj_num].begin(), d_var_names[obj_num].end(), var_name) != d_var_names[obj_num].end())
    {
        TBOX_ERROR(d_object_name << "::registerVariableData()\n"
                                 << "  variable with name " << var_name << " already registered for plotting\n"
                                 << "  on object number " << obj_num << std::endl);
    }
    ++d_nvars[obj_num];
    d_var_names[obj_num].push_back(var_name);
    d_var_start_depths[obj_num].push_back(start_depth);
    d_var_plot_depths[obj_num].push_back(var_depth);
    d_var_depths[obj_num].push_back(var_data->getDepth());
    d_var_data[obj_num].push_back(var_data);
    return;
} // registerVariableData

void
LSiloDataWriter::registerLagrangianAO(AO& ao, const int obj_num)
{
    if (obj_num == -1)
    {
        for (int i = 0; i < d_num_objs; ++i)
        {
            d_ao[i] = ao;
            d_build_vec_scatters[obj_num] = true;
        }
    }
    else
    {
        if (obj_num >= d_num_objs) resetNumObjs(obj_num + 1);

#if !defined(NDEBUG)
        TBOX_ASSERT(obj_num < d_num_objs);
#endif
        d_ao[obj_num] = ao;
        d_build_vec_scatters[obj_num] = true;
    }
    return;
} // registerLagrangianAO

void
LSiloDataWriter::registerLagrangianAO(std::vector<AO>& ao, std::vector<int>& obj_nums)
{
    for (auto obj_num : obj_nums) registerLagrangianAO(ao[obj_num], obj_num);
    return;
} // registerLagrangianAO

void
LSiloDataWriter::writePlotData(const int time_step_number, const double simulation_time)
{
#if defined(IBTK_HAVE_SILO)
#if !defined(NDEBUG)
    TBOX_ASSERT(time_step_number >= 0);
    TBOX_ASSERT(!d_dump_directory_name.empty());
#endif

    if (time_step_number <= d_time_step_number)
    {
        TBOX_ERROR(d_object_name << "::writePlotData()\n"
                                 << "  data writer with name " << d_object_name << "\n"
                                 << "  time step number: " << time_step_number
                                 << " is <= last time step number: " << d_time_step_number << std::endl);
    }
    d_time_step_number = time_step_number;

    if (d_dump_directory_name.empty())
    {
        TBOX_ERROR(d_object_name << "::writePlotData()\n"
                                 << "  data writer with name " << d_object_name << "\n"
                                 << "  dump directory name is empty" << std::endl);
    }

    int ierr;
    char temp_buf[SILO_NAME_BUFSIZE];
    std::string current_file_name;
    DBfile* dbfile;
    const int mpi_rank = IBTK_MPI::getRank();
    const int mpi_nodes = IBTK_MPI::getNodes();
    DBSetAllowEmptyObjects(1);

    // Construct the VecScatter objects required to write the plot data.
    for (int obj_num = 0; obj_num < d_num_objs; ++obj_num)
    {
        if (d_build_vec_scatters[obj_num])
        {
            buildVecScatters(d_ao[obj_num], obj_num);
        }
        d_build_vec_scatters[obj_num] = false;
    }

    // Create the working directory.
    std::snprintf(temp_buf, sizeof(temp_buf), "%06d", d_time_step_number);
    std::string current_dump_directory_name = SILO_DUMP_DIR_PREFIX + temp_buf;
    std::string dump_dirname = d_dump_directory_name + "/" + current_dump_directory_name;

    Utilities::recursiveMkdir(dump_dirname);

    // Create one local DBfile per MPI process.
    std::snprintf(temp_buf, sizeof(temp_buf), "%04d", mpi_rank);
    current_file_name = dump_dirname + "/" + SILO_PROCESSOR_FILE_PREFIX;
    current_file_name += temp_buf;
    current_file_name += SILO_PROCESSOR_FILE_POSTFIX;

    if (!(dbfile = DBCreate(current_file_name.c_str(), DB_CLOBBER, DB_LOCAL, nullptr, DB_PDB)))
    {
        TBOX_ERROR(d_object_name << "::writePlotData()\n"
                                 << "  Could not create DBfile named " << current_file_name << std::endl);
    }

    std::vector<std::vector<int> > meshtype(d_num_objs), vartype(d_num_objs);
    std::vector<std::vector<std::vector<int> > > multimeshtype(d_num_objs), multivartype(d_num_objs);

    // Set the local data.
    // We need to store the data across all objects drawn so we can draw cross structures.
    std::vector<Vec> local_X_vecs(d_num_objs);
    std::vector<double*> local_X_arrs(d_num_objs);
    std::vector<std::vector<Vec> > local_v_vecs_vec(d_num_objs);
    std::vector<std::vector<double*> > local_v_arrs_vec(d_num_objs);

    for (int obj_num = 0; obj_num < d_num_objs; ++obj_num)
    {
        if (d_coords_data[obj_num])
        {
            // Scatter the data from "global" to "local" form.
            ierr = VecDuplicate(d_dst_vec[obj_num][NDIM], &local_X_vecs[obj_num]);
            IBTK_CHKERRQ(ierr);

            Vec global_X_vec = d_coords_data[obj_num]->getVec();
            ierr = VecScatterBegin(
                d_vec_scatter[obj_num][NDIM], global_X_vec, local_X_vecs[obj_num], INSERT_VALUES, SCATTER_FORWARD);
            IBTK_CHKERRQ(ierr);
            ierr = VecScatterEnd(
                d_vec_scatter[obj_num][NDIM], global_X_vec, local_X_vecs[obj_num], INSERT_VALUES, SCATTER_FORWARD);
            IBTK_CHKERRQ(ierr);

            ierr = VecGetArray(local_X_vecs[obj_num], &local_X_arrs[obj_num]);
            IBTK_CHKERRQ(ierr);

            for (int v = 0; v < d_nvars[obj_num]; ++v)
            {
                const int var_depth = d_var_depths[obj_num][v];
                Vec local_v_vec;
                ierr = VecDuplicate(d_dst_vec[obj_num][var_depth], &local_v_vec);
                IBTK_CHKERRQ(ierr);

                Vec global_v_vec = d_var_data[obj_num][v]->getVec();
                ierr = VecScatterBegin(
                    d_vec_scatter[obj_num][var_depth], global_v_vec, local_v_vec, INSERT_VALUES, SCATTER_FORWARD);
                IBTK_CHKERRQ(ierr);
                ierr = VecScatterEnd(
                    d_vec_scatter[obj_num][var_depth], global_v_vec, local_v_vec, INSERT_VALUES, SCATTER_FORWARD);
                IBTK_CHKERRQ(ierr);

                double* local_v_arr;
                ierr = VecGetArray(local_v_vec, &local_v_arr);
                IBTK_CHKERRQ(ierr);

                local_v_vecs_vec[obj_num].push_back(local_v_vec);
                local_v_arrs_vec[obj_num].push_back(local_v_arr);
            }

            // Keep track of the current offset in the local Vec data.
            int offset = 0;

            // Add the local clouds to the local DBfile.
            for (int cloud = 0; cloud < d_nclouds[obj_num]; ++cloud)
            {
                const int nmarks = d_cloud_nmarks[obj_num][cloud];

                std::string dirname = "obj_num_" + std::to_string(obj_num) + "_cloud_" + std::to_string(cloud);

                if (DBMkDir(dbfile, dirname.c_str()) == -1)
                {
                    TBOX_ERROR(d_object_name << "::writePlotData()\n"
                                             << "  Could not create directory named " << dirname << std::endl);
                }

                const double* const X = local_X_arrs[obj_num] + NDIM * offset;
                std::vector<const double*> var_vals(d_nvars[obj_num]);
                for (int v = 0; v < d_nvars[obj_num]; ++v)
                {
                    var_vals[v] = local_v_arrs_vec[obj_num][v] + d_var_depths[obj_num][v] * offset;
                }

                build_local_marker_cloud(dbfile,
                                         dirname,
                                         nmarks,
                                         X,
                                         d_nvars[obj_num],
                                         d_var_names[obj_num],
                                         d_var_start_depths[obj_num],
                                         d_var_plot_depths[obj_num],
                                         d_var_depths[obj_num],
                                         var_vals,
                                         time_step_number,
                                         simulation_time);

                offset += nmarks;
            }

            // Add the local blocks to the local DBfile.
            for (int block = 0; block < d_nblocks[obj_num]; ++block)
            {
                const IntVector<NDIM>& nelem = d_block_nelems[obj_num][block];
                const IntVector<NDIM>& periodic = d_block_periodic[obj_num][block];
                const int ntot = nelem.getProduct();

                std::string dirname = "obj_num_" + std::to_string(obj_num) + "_block_" + std::to_string(block);

                if (DBMkDir(dbfile, dirname.c_str()) == -1)
                {
                    TBOX_ERROR(d_object_name << "::writePlotData()\n"
                                             << "  Could not create directory named " << dirname << std::endl);
                }

                const double* const X = local_X_arrs[obj_num] + NDIM * offset;
                std::vector<const double*> var_vals(d_nvars[obj_num]);
                for (int v = 0; v < d_nvars[obj_num]; ++v)
                {
                    var_vals[v] = local_v_arrs_vec[obj_num][v] + d_var_depths[obj_num][v] * offset;
                }

                build_local_curv_block(dbfile,
                                       dirname,
                                       nelem,
                                       periodic,
                                       X,
                                       d_nvars[obj_num],
                                       d_var_names[obj_num],
                                       d_var_start_depths[obj_num],
                                       d_var_plot_depths[obj_num],
                                       d_var_depths[obj_num],
                                       var_vals,
                                       time_step_number,
                                       simulation_time);
                meshtype[obj_num].push_back(DB_QUAD_CURV);
                vartype[obj_num].push_back(DB_QUADVAR);

                offset += ntot;
            }

            // Add the local multiblocks to the local DBfile.
            multimeshtype[obj_num].resize(d_nmbs[obj_num]);
            multivartype[obj_num].resize(d_nmbs[obj_num]);
            for (int mb = 0; mb < d_nmbs[obj_num]; ++mb)
            {
                for (int block = 0; block < d_mb_nblocks[obj_num][mb]; ++block)
                {
                    const IntVector<NDIM>& nelem = d_mb_nelems[obj_num][mb][block];
                    const IntVector<NDIM>& periodic = d_mb_periodic[obj_num][mb][block];
                    const int ntot = nelem.getProduct();

                    std::string dirname = "obj_num_" + std::to_string(obj_num) + "_mb_" + std::to_string(mb) +
                                          "_block_" + std::to_string(block);

                    if (DBMkDir(dbfile, dirname.c_str()) == -1)
                    {
                        TBOX_ERROR(d_object_name << "::writePlotData()\n"
                                                 << "  Could not create directory named " << dirname << std::endl);
                    }

                    const double* const X = local_X_arrs[obj_num] + NDIM * offset;
                    std::vector<const double*> var_vals(d_nvars[obj_num]);
                    for (int v = 0; v < d_nvars[obj_num]; ++v)
                    {
                        var_vals[v] = local_v_arrs_vec[obj_num][v] + d_var_depths[obj_num][v] * offset;
                    }

                    build_local_curv_block(dbfile,
                                           dirname,
                                           nelem,
                                           periodic,
                                           X,
                                           d_nvars[obj_num],
                                           d_var_names[obj_num],
                                           d_var_start_depths[obj_num],
                                           d_var_plot_depths[obj_num],
                                           d_var_depths[obj_num],
                                           var_vals,
                                           time_step_number,
                                           simulation_time);
                    multimeshtype[obj_num][mb].push_back(DB_QUAD_CURV);
                    multivartype[obj_num][mb].push_back(DB_QUADVAR);

                    offset += ntot;
                }
            }

            // Add the local UCD meshes to the local DBfile.
            for (int mesh = 0; mesh < d_nucd_meshes[obj_num]; ++mesh)
            {
                const std::set<int>& vertices = d_ucd_mesh_vertices[obj_num][mesh];
                const std::multimap<int, std::pair<int, int> >& edge_map = d_ucd_mesh_edge_maps[obj_num][mesh];
                const size_t ntot = vertices.size();

                std::string dirname = "obj_num_" + std::to_string(obj_num) + "_mesh_" + std::to_string(mesh);

                if (DBMkDir(dbfile, dirname.c_str()) == -1)
                {
                    TBOX_ERROR(d_object_name << "::writePlotData()\n"
                                             << "  Could not create directory named " << dirname << std::endl);
                }

                const double* const X = local_X_arrs[obj_num] + NDIM * offset;
                std::vector<const double*> var_vals(d_nvars[obj_num]);
                for (int v = 0; v < d_nvars[obj_num]; ++v)
                {
                    var_vals[v] = local_v_arrs_vec[obj_num][v] + d_var_depths[obj_num][v] * offset;
                }

                build_local_ucd_mesh(dbfile,
                                     dirname,
                                     vertices,
                                     edge_map,
                                     X,
                                     d_nvars[obj_num],
                                     d_var_names[obj_num],
                                     d_var_start_depths[obj_num],
                                     d_var_plot_depths[obj_num],
                                     d_var_depths[obj_num],
                                     var_vals,
                                     time_step_number,
                                     simulation_time);

                offset += ntot;
            }
        }
    }

    // Now build local ucd cross meshes.
    for (int mesh = 0; mesh < d_nucd_cross_meshes; ++mesh)
    {
        const std::set<int>& cross_obj_nums = d_ucd_cross_mesh_obj_nums[mesh];
        const std::map<std::pair<int, int>, std::set<int> >& cross_vertices = d_ucd_cross_mesh_vertices[mesh];
        const std::vector<CrossEdge>& cross_edges = d_ucd_cross_edges[mesh];

        std::string dirname = "cross_obj_num_" + std::to_string(mesh);
        if (DBMkDir(dbfile, dirname.c_str()) == -1)
        {
            TBOX_ERROR(d_object_name << "::writePlotData()\n"
                                     << "   Could not create directory named " << dirname << std::endl);
        }

        int ntot = std::accumulate(cross_vertices.begin(),
                                   cross_vertices.end(),
                                   0,
                                   [](int value, const std::pair<std::pair<int, int>, std::set<int> >& p) -> int {
                                       return value + p.second.size();
                                   });
        // Compute the coordinates.
        std::vector<float> block_X(NDIM * ntot);
        int block_offset = 0;
        std::map<std::pair<std::pair<int, int>, int>, int> local_vertex_map;
        for (const auto& cross_vertex : cross_vertices)
        {
            int obj_num = cross_vertex.first.first;
            int cloud_num = cross_vertex.first.second;
            const std::set<int>& idxs = cross_vertex.second;

            const double* const X = local_X_arrs[obj_num];
            int offset =
                std::accumulate(d_cloud_nmarks[obj_num].begin(), d_cloud_nmarks[obj_num].begin() + cloud_num, 0);
            // Compute offset.
            for (const auto& idx : idxs)
            {
                for (int d = 0; d < NDIM; ++d)
                    block_X[d * ntot + block_offset] = static_cast<float>(X[NDIM * (offset + idx) + d]);
                auto vertex = std::make_pair(cross_vertex.first, idx);
                local_vertex_map[vertex] = block_offset;
                ++block_offset;
            }
        }

        // Now compute a pruned local edge list.
        // Note the integers in this map correspond to locations in block_X
        std::multimap<int, int> local_edge_map;
        for (const auto& cross_edge : cross_edges)
        {
            // Determine first link
            int e1 = local_vertex_map[cross_edge.pt1];
            int e2 = local_vertex_map[cross_edge.pt2];
            local_edge_map.emplace(e1, e2);
        }

        // Set the working directory
        if (DBSetDir(dbfile, dirname.c_str()) == -1)
            TBOX_ERROR(d_object_name << "::writePlotData()\n"
                                     << "   Could not set directory " << dirname << std::endl);

        // Now setup the silo arguments
        int ndims = NDIM;
        int cycle_num = time_step_number;
        auto time = static_cast<float>(simulation_time);
        double dtime = simulation_time;
        static const int MAX_OPTS = 3;
        DBoptlist* optlist = DBMakeOptlist(MAX_OPTS);
        DBAddOption(optlist, DBOPT_CYCLE, &cycle_num);
        DBAddOption(optlist, DBOPT_TIME, &time);
        DBAddOption(optlist, DBOPT_DTIME, &dtime);

        const char* meshname = "mesh";
        const char* coordnames[3] = { "xcoords", "ycoords", "zcoords" };
        std::vector<float*> coords(NDIM);
        for (int d = 0; d < NDIM; ++d) coords[d] = ntot > 0 ? &block_X[d * ntot] : nullptr;
        const int nnodes = ntot;

        std::vector<int> nodelist;
        nodelist.reserve(2 * local_edge_map.size());
        for (const auto& edge_pair : local_edge_map)
        {
            nodelist.push_back(edge_pair.first);
            nodelist.push_back(edge_pair.second);
        }
        int lnodelist = static_cast<int>(nodelist.size());
        int nshapetypes = 1;
        int shapecnt[] = { static_cast<int>(local_edge_map.size()) };
        int shapesize[] = { 2 };
        int shapetype[] = { DB_ZONETYPE_BEAM };
        int nzones = static_cast<int>(local_edge_map.size());

        const int origin = 0;
        const int lo_offset = 0;
        const int hi_offset = 0;

        // Write file
        DBPutZonelist2(dbfile,
                       "zonelist",
                       nzones,
                       ndims,
                       lnodelist > 0 ? &nodelist[0] : nullptr,
                       lnodelist,
                       origin,
                       lo_offset,
                       hi_offset,
                       shapetype,
                       shapesize,
                       shapecnt,
                       nshapetypes,
                       optlist);
        DBPutUcdmesh(dbfile,
                     meshname,
                     ndims,
                     const_cast<char**>(coordnames),
                     &coords[0],
                     nnodes,
                     nzones,
                     "zonelist",
                     nullptr,
                     DB_FLOAT,
                     nullptr);

        DBFreeOptlist(optlist);

        // Reset working directory
        if (DBSetDir(dbfile, "..") == -1)
            TBOX_ERROR(d_object_name << "::writePlotData()\n"
                                     << "   Could not return to the base directory from subdirectory " << dirname
                                     << std::endl);
    }

    // Clean up allocated data.
    for (int obj_num = 0; obj_num < d_num_objs; ++obj_num)
    {
        if (d_coords_data[obj_num])
        {
            ierr = VecRestoreArray(local_X_vecs[obj_num], &local_X_arrs[obj_num]);
            IBTK_CHKERRQ(ierr);
            ierr = VecDestroy(&local_X_vecs[obj_num]);
            IBTK_CHKERRQ(ierr);
            for (int v = 0; v < d_nvars[obj_num]; ++v)
            {
                ierr = VecRestoreArray(local_v_vecs_vec[obj_num][v], &local_v_arrs_vec[obj_num][v]);
                IBTK_CHKERRQ(ierr);
                ierr = VecDestroy(&local_v_vecs_vec[obj_num][v]);
                IBTK_CHKERRQ(ierr);
            }
        }
    }

    DBClose(dbfile);

    // Send data to the root MPI process required to create the multimesh and
    // multivar objects.
    std::vector<std::vector<int> > nclouds_per_proc, nblocks_per_proc, nmbs_per_proc, nucd_meshes_per_proc;
    std::vector<std::vector<std::vector<int> > > meshtypes_per_proc, vartypes_per_proc, mb_nblocks_per_proc;
    std::vector<std::vector<std::vector<std::vector<int> > > > multimeshtypes_per_proc, multivartypes_per_proc;
    std::vector<std::vector<std::vector<std::string> > > cloud_names_per_proc, block_names_per_proc, mb_names_per_proc,
        ucd_mesh_names_per_proc;

    if (mpi_rank == SILO_MPI_ROOT)
    {
        nclouds_per_proc.resize(d_num_objs);
        nblocks_per_proc.resize(d_num_objs);
        nmbs_per_proc.resize(d_num_objs);
        nucd_meshes_per_proc.resize(d_num_objs);
        meshtypes_per_proc.resize(d_num_objs);
        vartypes_per_proc.resize(d_num_objs);
        mb_nblocks_per_proc.resize(d_num_objs);
        multimeshtypes_per_proc.resize(d_num_objs);
        multivartypes_per_proc.resize(d_num_objs);
        cloud_names_per_proc.resize(d_num_objs);
        block_names_per_proc.resize(d_num_objs);
        mb_names_per_proc.resize(d_num_objs);
        ucd_mesh_names_per_proc.resize(d_num_objs);
    }

    for (int obj_num = 0; obj_num < d_num_objs; ++obj_num)
    {
        if (mpi_rank == SILO_MPI_ROOT)
        {
            nclouds_per_proc[obj_num].resize(mpi_nodes);
            nblocks_per_proc[obj_num].resize(mpi_nodes);
            nmbs_per_proc[obj_num].resize(mpi_nodes);
            nucd_meshes_per_proc[obj_num].resize(mpi_nodes);
            meshtypes_per_proc[obj_num].resize(mpi_nodes);
            vartypes_per_proc[obj_num].resize(mpi_nodes);
            mb_nblocks_per_proc[obj_num].resize(mpi_nodes);
            multimeshtypes_per_proc[obj_num].resize(mpi_nodes);
            multivartypes_per_proc[obj_num].resize(mpi_nodes);
            cloud_names_per_proc[obj_num].resize(mpi_nodes);
            block_names_per_proc[obj_num].resize(mpi_nodes);
            mb_names_per_proc[obj_num].resize(mpi_nodes);
            ucd_mesh_names_per_proc[obj_num].resize(mpi_nodes);

            // Set the values for the root process.
            nclouds_per_proc[obj_num][mpi_rank] = d_nclouds[obj_num];
            nblocks_per_proc[obj_num][mpi_rank] = d_nblocks[obj_num];
            nmbs_per_proc[obj_num][mpi_rank] = d_nmbs[obj_num];
            nucd_meshes_per_proc[obj_num][mpi_rank] = d_nucd_meshes[obj_num];
            meshtypes_per_proc[obj_num][mpi_rank] = meshtype[obj_num];
            vartypes_per_proc[obj_num][mpi_rank] = vartype[obj_num];
            mb_nblocks_per_proc[obj_num][mpi_rank] = d_mb_nblocks[obj_num];
            multimeshtypes_per_proc[obj_num][mpi_rank] = multimeshtype[obj_num];
            multivartypes_per_proc[obj_num][mpi_rank] = multivartype[obj_num];
            cloud_names_per_proc[obj_num][mpi_rank] = d_cloud_names[obj_num];
            block_names_per_proc[obj_num][mpi_rank] = d_block_names[obj_num];
            mb_names_per_proc[obj_num][mpi_rank] = d_mb_names[obj_num];
            ucd_mesh_names_per_proc[obj_num][mpi_rank] = d_ucd_mesh_names[obj_num];
        }

        // Get the values for the non-root processes.
        int one = 1;
        for (int proc = 0; proc < mpi_nodes; ++proc)
        {
            // Skip the root process; we already have those values.
            if (proc == SILO_MPI_ROOT)
            {
                proc += 1;
                if (proc >= mpi_nodes) break;
            }

            if (mpi_rank == proc)
            {
                IBTK_MPI::send(&d_nclouds[obj_num], one, SILO_MPI_ROOT, false);
            }
            if (mpi_rank == SILO_MPI_ROOT)
            {
                IBTK_MPI::recv(&nclouds_per_proc[obj_num][proc], one, proc, false);
            }

            if (mpi_rank == proc && d_nclouds[obj_num] > 0)
            {
                int num_bytes;
                for (int cloud = 0; cloud < d_nclouds[obj_num]; ++cloud)
                {
                    num_bytes = static_cast<int>((d_cloud_names[obj_num][cloud].size() + 1) * sizeof(char));
                    IBTK_MPI::send(&num_bytes, one, SILO_MPI_ROOT, false);
                    IBTK_MPI::sendBytes(
                        static_cast<const void*>(d_cloud_names[obj_num][cloud].c_str()), num_bytes, SILO_MPI_ROOT);
                }
            }
            if (mpi_rank == SILO_MPI_ROOT && nclouds_per_proc[obj_num][proc] > 0)
            {
                cloud_names_per_proc[obj_num][proc].resize(nclouds_per_proc[obj_num][proc]);

                int num_bytes;
                for (int cloud = 0; cloud < nclouds_per_proc[obj_num][proc]; ++cloud)
                {
                    IBTK_MPI::recv(&num_bytes, one, proc, false);
                    std::vector<char> name(num_bytes / sizeof(char));
                    IBTK_MPI::recvBytes(static_cast<void*>(name.data()), num_bytes);
                    cloud_names_per_proc[obj_num][proc][cloud].assign(name.data());
                }
            }

            if (mpi_rank == proc)
            {
                IBTK_MPI::send(&d_nblocks[obj_num], one, SILO_MPI_ROOT, false);
            }
            if (mpi_rank == SILO_MPI_ROOT)
            {
                IBTK_MPI::recv(&nblocks_per_proc[obj_num][proc], one, proc, false);
            }

            if (mpi_rank == proc && d_nblocks[obj_num] > 0)
            {
                IBTK_MPI::send(&meshtype[obj_num][0], d_nblocks[obj_num], SILO_MPI_ROOT, false);
                IBTK_MPI::send(&vartype[obj_num][0], d_nblocks[obj_num], SILO_MPI_ROOT, false);

                int num_bytes;
                for (int block = 0; block < d_nblocks[obj_num]; ++block)
                {
                    num_bytes = static_cast<int>((d_block_names[obj_num][block].size() + 1) * sizeof(char));
                    IBTK_MPI::send(&num_bytes, one, SILO_MPI_ROOT, false);
                    IBTK_MPI::sendBytes(
                        static_cast<const void*>(d_block_names[obj_num][block].c_str()), num_bytes, SILO_MPI_ROOT);
                }
            }
            if (mpi_rank == SILO_MPI_ROOT && nblocks_per_proc[obj_num][proc] > 0)
            {
                meshtypes_per_proc[obj_num][proc].resize(nblocks_per_proc[obj_num][proc]);
                vartypes_per_proc[obj_num][proc].resize(nblocks_per_proc[obj_num][proc]);
                block_names_per_proc[obj_num][proc].resize(nblocks_per_proc[obj_num][proc]);

                IBTK_MPI::recv(&meshtypes_per_proc[obj_num][proc][0], nblocks_per_proc[obj_num][proc], proc, false);
                IBTK_MPI::recv(&vartypes_per_proc[obj_num][proc][0], nblocks_per_proc[obj_num][proc], proc, false);

                int num_bytes;
                for (int block = 0; block < nblocks_per_proc[obj_num][proc]; ++block)
                {
                    IBTK_MPI::recv(&num_bytes, one, proc, false);
                    std::vector<char> name(num_bytes / sizeof(char));
                    IBTK_MPI::recvBytes(static_cast<void*>(name.data()), num_bytes);
                    block_names_per_proc[obj_num][proc][block].assign(name.data());
                }
            }

            if (mpi_rank == proc)
            {
                IBTK_MPI::send(&d_nmbs[obj_num], one, SILO_MPI_ROOT, false);
            }
            if (mpi_rank == SILO_MPI_ROOT)
            {
                IBTK_MPI::recv(&nmbs_per_proc[obj_num][proc], one, proc, false);
            }

            if (mpi_rank == proc && d_nmbs[obj_num] > 0)
            {
                IBTK_MPI::send(&d_mb_nblocks[obj_num][0], d_nmbs[obj_num], SILO_MPI_ROOT, false);

                int num_bytes;
                for (int mb = 0; mb < d_nmbs[obj_num]; ++mb)
                {
                    IBTK_MPI::send(&multimeshtype[obj_num][mb][0], d_mb_nblocks[obj_num][mb], SILO_MPI_ROOT, false);
                    IBTK_MPI::send(&multivartype[obj_num][mb][0], d_mb_nblocks[obj_num][mb], SILO_MPI_ROOT, false);

                    num_bytes = static_cast<int>(d_mb_names[obj_num][mb].size()) + 1;
                    IBTK_MPI::send(&num_bytes, one, SILO_MPI_ROOT, false);

                    const void* mb_name_ptr = d_mb_names[obj_num][mb].c_str();
                    MPI_Send(const_cast<void*>(mb_name_ptr),
                             num_bytes,
                             MPI_CHAR,
                             SILO_MPI_ROOT,
                             SILO_MPI_TAG,
                             IBTK_MPI::getCommunicator());
                }
            }
            if (mpi_rank == SILO_MPI_ROOT && nmbs_per_proc[obj_num][proc] > 0)
            {
                mb_nblocks_per_proc[obj_num][proc].resize(nmbs_per_proc[obj_num][proc]);
                multimeshtypes_per_proc[obj_num][proc].resize(nmbs_per_proc[obj_num][proc]);
                multivartypes_per_proc[obj_num][proc].resize(nmbs_per_proc[obj_num][proc]);
                mb_names_per_proc[obj_num][proc].resize(nmbs_per_proc[obj_num][proc]);

                IBTK_MPI::recv(&mb_nblocks_per_proc[obj_num][proc][0], nmbs_per_proc[obj_num][proc], proc, false);

                int num_bytes;
                for (int mb = 0; mb < nmbs_per_proc[obj_num][proc]; ++mb)
                {
                    multimeshtypes_per_proc[obj_num][proc][mb].resize(mb_nblocks_per_proc[obj_num][proc][mb]);
                    multivartypes_per_proc[obj_num][proc][mb].resize(mb_nblocks_per_proc[obj_num][proc][mb]);

                    IBTK_MPI::recv(&multimeshtypes_per_proc[obj_num][proc][mb][0],
                                   mb_nblocks_per_proc[obj_num][proc][mb],
                                   proc,
                                   false);
                    IBTK_MPI::recv(&multivartypes_per_proc[obj_num][proc][mb][0],
                                   mb_nblocks_per_proc[obj_num][proc][mb],
                                   proc,
                                   false);

                    IBTK_MPI::recv(&num_bytes, one, proc, false);
                    std::vector<char> name(num_bytes / sizeof(char));

                    MPI_Status status;
                    MPI_Recv(
                        name.data(), num_bytes, MPI_CHAR, proc, SILO_MPI_TAG, IBTK_MPI::getCommunicator(), &status);

                    mb_names_per_proc[obj_num][proc][mb].assign(name.data());
                }
            }

            if (mpi_rank == proc)
            {
                IBTK_MPI::send(&d_nucd_meshes[obj_num], one, SILO_MPI_ROOT, false);
            }
            if (mpi_rank == SILO_MPI_ROOT)
            {
                IBTK_MPI::recv(&nucd_meshes_per_proc[obj_num][proc], one, proc, false);
            }

            if (mpi_rank == proc && d_nucd_meshes[obj_num] > 0)
            {
                int num_bytes;
                for (int mesh = 0; mesh < d_nucd_meshes[obj_num]; ++mesh)
                {
                    num_bytes = static_cast<int>((d_ucd_mesh_names[obj_num][mesh].size() + 1) * sizeof(char));
                    IBTK_MPI::send(&num_bytes, one, SILO_MPI_ROOT, false);
                    IBTK_MPI::sendBytes(
                        static_cast<const void*>(d_ucd_mesh_names[obj_num][mesh].c_str()), num_bytes, SILO_MPI_ROOT);
                }
            }
            if (mpi_rank == SILO_MPI_ROOT && nucd_meshes_per_proc[obj_num][proc] > 0)
            {
                ucd_mesh_names_per_proc[obj_num][proc].resize(nucd_meshes_per_proc[obj_num][proc]);
                int num_bytes;
                for (int mesh = 0; mesh < nucd_meshes_per_proc[obj_num][proc]; ++mesh)
                {
                    IBTK_MPI::recv(&num_bytes, one, proc, false);
                    std::vector<char> name(num_bytes / sizeof(char));
                    IBTK_MPI::recvBytes(static_cast<void*>(name.data()), num_bytes);
                    ucd_mesh_names_per_proc[obj_num][proc][mesh].assign(name.data());
                }
            }

            IBTK_MPI::barrier();
        }
    }

    std::vector<int> nucd_cross_meshes_per_proc;
    std::vector<std::vector<std::string> > ucd_cross_mesh_names_per_proc;
    if (mpi_rank == SILO_MPI_ROOT)
    {
        // Fill in data for the root process
        nucd_cross_meshes_per_proc.resize(mpi_nodes);
        ucd_cross_mesh_names_per_proc.resize(mpi_nodes);

        nucd_cross_meshes_per_proc[mpi_rank] = d_nucd_cross_meshes;
        ucd_cross_mesh_names_per_proc[mpi_rank] = d_ucd_cross_mesh_names;
    }

    int one = 1;
    for (int proc = 0; proc < mpi_nodes; ++proc)
    {
        // Skip root process.
        if (proc == SILO_MPI_ROOT)
        {
            proc += 1;
            if (proc >= mpi_nodes) break;
        }

        if (mpi_rank == proc) IBTK_MPI::send(&d_nucd_cross_meshes, one, SILO_MPI_ROOT, false);
        if (mpi_rank == SILO_MPI_ROOT) IBTK_MPI::recv(&nucd_cross_meshes_per_proc[proc], one, proc, false);

        if (mpi_rank == proc && d_nucd_cross_meshes > 0)
        {
            int num_bytes = 0;
            for (int mesh = 0; mesh < d_nucd_cross_meshes; ++mesh)
            {
                num_bytes = static_cast<int>((d_ucd_cross_mesh_names[mesh].size() + 1) * sizeof(char));
                IBTK_MPI::send(&num_bytes, one, SILO_MPI_ROOT, false);
                IBTK_MPI::sendBytes(
                    static_cast<const void*>(d_ucd_cross_mesh_names[mesh].c_str()), num_bytes, SILO_MPI_ROOT);
            }
        }
        if (mpi_rank == SILO_MPI_ROOT && nucd_cross_meshes_per_proc[proc] > 0)
        {
            ucd_cross_mesh_names_per_proc[proc].resize(nucd_cross_meshes_per_proc[proc]);
            int num_bytes;
            for (int mesh = 0; mesh < nucd_cross_meshes_per_proc[proc]; ++mesh)
            {
                IBTK_MPI::recv(&num_bytes, one, proc, false);
                std::vector<char> name(num_bytes / sizeof(char));
                IBTK_MPI::recvBytes(static_cast<void*>(name.data()), num_bytes);
                ucd_cross_mesh_names_per_proc[proc][mesh].assign(name.data());
            }
        }
    }
    IBTK_MPI::barrier();

    if (mpi_rank == SILO_MPI_ROOT)
    {
        // Create and initialize the multimesh Silo database on the root MPI
        // process.
        std::snprintf(temp_buf, sizeof(temp_buf), "%06d", d_time_step_number);
        std::string summary_file_name =
            dump_dirname + "/" + SILO_SUMMARY_FILE_PREFIX + temp_buf + SILO_SUMMARY_FILE_POSTFIX;
        if (!(dbfile = DBCreate(summary_file_name.c_str(), DB_CLOBBER, DB_LOCAL, nullptr, DB_PDB)))
        {
            TBOX_ERROR(d_object_name << "::writePlotData()\n"
                                     << "  Could not create DBfile named " << summary_file_name << std::endl);
        }

        int cycle = time_step_number;
        auto time = static_cast<float>(simulation_time);
        double dtime = simulation_time;

        static const int MAX_OPTS = 3;
        DBoptlist* optlist = DBMakeOptlist(MAX_OPTS);
        DBAddOption(optlist, DBOPT_CYCLE, &cycle);
        DBAddOption(optlist, DBOPT_TIME, &time);
        DBAddOption(optlist, DBOPT_DTIME, &dtime);

        for (int proc = 0; proc < mpi_nodes; ++proc)
        {
            for (int obj_num = 0; obj_num < d_num_objs; ++obj_num)
            {
                for (int cloud = 0; cloud < nclouds_per_proc[obj_num][proc]; ++cloud)
                {
                    std::snprintf(temp_buf, sizeof(temp_buf), "%04d", proc);
                    current_file_name = SILO_PROCESSOR_FILE_PREFIX;
                    current_file_name += temp_buf;
                    current_file_name += SILO_PROCESSOR_FILE_POSTFIX;

                    std::string meshname = current_file_name + ":obj_num_" + std::to_string(obj_num) + "_cloud_" +
                                           std::to_string(cloud) + "/mesh";
                    auto meshname_ptr = const_cast<char*>(meshname.c_str());
                    int meshtype = DB_POINTMESH;

                    std::string& cloud_name = cloud_names_per_proc[obj_num][proc][cloud];

                    DBPutMultimesh(dbfile, cloud_name.c_str(), 1, &meshname_ptr, &meshtype, optlist);

                    if (DBMkDir(dbfile, cloud_name.c_str()) == -1)
                    {
                        TBOX_ERROR(d_object_name << "::writePlotData()\n"
                                                 << "  Could not create directory named " << cloud_name << std::endl);
                    }
                }

                for (int block = 0; block < nblocks_per_proc[obj_num][proc]; ++block)
                {
                    std::snprintf(temp_buf, sizeof(temp_buf), "%04d", proc);
                    current_file_name = SILO_PROCESSOR_FILE_PREFIX;
                    current_file_name += temp_buf;
                    current_file_name += SILO_PROCESSOR_FILE_POSTFIX;

                    std::string meshname = current_file_name + ":obj_num_" + std::to_string(obj_num) + "_block_" +
                                           std::to_string(block) + "/mesh";
                    auto meshname_ptr = const_cast<char*>(meshname.c_str());
                    int meshtype = meshtypes_per_proc[obj_num][proc][block];

                    std::string& block_name = block_names_per_proc[obj_num][proc][block];

                    DBPutMultimesh(dbfile, block_name.c_str(), 1, &meshname_ptr, &meshtype, optlist);

                    if (DBMkDir(dbfile, block_name.c_str()) == -1)
                    {
                        TBOX_ERROR(d_object_name << "::writePlotData()\n"
                                                 << "  Could not create directory named " << block_name << std::endl);
                    }
                }

                for (int mb = 0; mb < nmbs_per_proc[obj_num][proc]; ++mb)
                {
                    std::snprintf(temp_buf, sizeof(temp_buf), "%04d", proc);
                    current_file_name = SILO_PROCESSOR_FILE_PREFIX;
                    current_file_name += temp_buf;
                    current_file_name += SILO_PROCESSOR_FILE_POSTFIX;

                    const int nblocks = mb_nblocks_per_proc[obj_num][proc][mb];
                    std::vector<std::string> meshnames;
                    for (int block = 0; block < nblocks; ++block)
                    {
                        meshnames.push_back(current_file_name + ":obj_num_" + std::to_string(obj_num) + "_mb_" +
                                            std::to_string(mb) + "_block_" + std::to_string(block) + "/mesh");
                    }
                    std::vector<const char*> meshnames_ptrs;
                    for (int block = 0; block < nblocks; ++block)
                    {
                        meshnames_ptrs.push_back(meshnames[block].c_str());
                    }

                    std::string& mb_name = mb_names_per_proc[obj_num][proc][mb];

                    DBPutMultimesh(dbfile,
                                   mb_name.c_str(),
                                   nblocks,
                                   meshnames_ptrs.data(),
                                   &multimeshtypes_per_proc[obj_num][proc][mb][0],
                                   optlist);

                    if (DBMkDir(dbfile, mb_name.c_str()) == -1)
                    {
                        TBOX_ERROR(d_object_name << "::writePlotData()\n"
                                                 << "  Could not create directory named " << mb_name << std::endl);
                    }
                }

                for (int mesh = 0; mesh < nucd_meshes_per_proc[obj_num][proc]; ++mesh)
                {
                    std::snprintf(temp_buf, sizeof(temp_buf), "%04d", proc);
                    current_file_name = SILO_PROCESSOR_FILE_PREFIX;
                    current_file_name += temp_buf;
                    current_file_name += SILO_PROCESSOR_FILE_POSTFIX;

                    std::string meshname = current_file_name + ":obj_num_" + std::to_string(obj_num) + "_mesh_" +
                                           std::to_string(mesh) + "/mesh";
                    auto meshname_ptr = const_cast<char*>(meshname.c_str());
                    int meshtype = DB_UCDMESH;

                    std::string& mesh_name = ucd_mesh_names_per_proc[obj_num][proc][mesh];

                    DBPutMultimesh(dbfile, mesh_name.c_str(), 1, &meshname_ptr, &meshtype, optlist);

                    if (DBMkDir(dbfile, mesh_name.c_str()) == -1)
                    {
                        TBOX_ERROR(d_object_name << "::writePlotData()\n"
                                                 << "  Could not create directory named " << mesh_name << std::endl);
                    }
                }

                for (int v = 0; v < d_nvars[obj_num]; ++v)
                {
                    for (int cloud = 0; cloud < nclouds_per_proc[obj_num][proc]; ++cloud)
                    {
                        std::snprintf(temp_buf, sizeof(temp_buf), "%04d", proc);
                        current_file_name = SILO_PROCESSOR_FILE_PREFIX;
                        current_file_name += temp_buf;
                        current_file_name += SILO_PROCESSOR_FILE_POSTFIX;

                        std::string varname = current_file_name + ":obj_num_" + std::to_string(obj_num) + "_cloud_" +
                                              std::to_string(cloud) + "/" + d_var_names[obj_num][v];
                        auto varname_ptr = const_cast<char*>(varname.c_str());
                        int vartype = DB_POINTVAR;

                        std::string& cloud_name = cloud_names_per_proc[obj_num][proc][cloud];

                        std::string var_name = cloud_name + "/" + d_var_names[obj_num][v];

                        DBPutMultivar(dbfile, var_name.c_str(), 1, &varname_ptr, &vartype, optlist);
                    }

                    for (int block = 0; block < nblocks_per_proc[obj_num][proc]; ++block)
                    {
                        std::snprintf(temp_buf, sizeof(temp_buf), "%04d", proc);
                        current_file_name = SILO_PROCESSOR_FILE_PREFIX;
                        current_file_name += temp_buf;
                        current_file_name += SILO_PROCESSOR_FILE_POSTFIX;

                        std::string varname = current_file_name + ":obj_num_" + std::to_string(obj_num) + "_block_" +
                                              std::to_string(block) + "/" + d_var_names[obj_num][v];
                        auto varname_ptr = const_cast<char*>(varname.c_str());
                        int vartype = vartypes_per_proc[obj_num][proc][block];

                        std::string& block_name = block_names_per_proc[obj_num][proc][block];

                        std::string var_name = block_name + "/" + d_var_names[obj_num][v];

                        DBPutMultivar(dbfile, var_name.c_str(), 1, &varname_ptr, &vartype, optlist);
                    }

                    for (int mb = 0; mb < nmbs_per_proc[obj_num][proc]; ++mb)
                    {
                        std::snprintf(temp_buf, sizeof(temp_buf), "%04d", proc);
                        current_file_name = SILO_PROCESSOR_FILE_PREFIX;
                        current_file_name += temp_buf;
                        current_file_name += SILO_PROCESSOR_FILE_POSTFIX;

                        const int nblocks = mb_nblocks_per_proc[obj_num][proc][mb];

                        std::vector<std::string> varnames;
                        for (int block = 0; block < nblocks; ++block)
                        {
                            varnames.push_back(current_file_name + ":obj_num_" + std::to_string(obj_num) + "_mb_" +
                                               std::to_string(mb) + "_block_" + std::to_string(block) +
                                               d_var_names[obj_num][v]);
                        }
                        std::vector<const char*> varnames_ptrs;
                        for (int block = 0; block < nblocks; ++block)
                        {
                            varnames_ptrs.push_back(varnames[block].c_str());
                        }

                        std::string& mb_name = mb_names_per_proc[obj_num][proc][mb];

                        std::string var_name = mb_name + "/" + d_var_names[obj_num][v];

                        DBPutMultivar(dbfile,
                                      var_name.c_str(),
                                      nblocks,
                                      varnames_ptrs.data(),
                                      &multivartypes_per_proc[obj_num][proc][mb][0],
                                      optlist);
                    }

                    for (int mesh = 0; mesh < nucd_meshes_per_proc[obj_num][proc]; ++mesh)
                    {
                        std::snprintf(temp_buf, sizeof(temp_buf), "%04d", proc);
                        current_file_name = SILO_PROCESSOR_FILE_PREFIX;
                        current_file_name += temp_buf;
                        current_file_name += SILO_PROCESSOR_FILE_POSTFIX;

                        std::string varname = current_file_name + ":obj_num_" + std::to_string(obj_num) + "_mesh_" +
                                              std::to_string(mesh) + "/" + d_var_names[obj_num][v];
                        auto varname_ptr = const_cast<char*>(varname.c_str());
                        int vartype = DB_UCDVAR;

                        std::string& mesh_name = ucd_mesh_names_per_proc[obj_num][proc][mesh];

                        std::string var_name = mesh_name + "/" + d_var_names[obj_num][v];

                        DBPutMultivar(dbfile, var_name.c_str(), 1, &varname_ptr, &vartype, optlist);
                    }
                }
            }

            for (int mesh = 0; mesh < nucd_cross_meshes_per_proc[proc]; ++mesh)
            {
                std::snprintf(temp_buf, sizeof(temp_buf), "%04d", proc);
                current_file_name = SILO_PROCESSOR_FILE_PREFIX;
                current_file_name += temp_buf;
                current_file_name += SILO_PROCESSOR_FILE_POSTFIX;

                std::string meshname = current_file_name + ":cross_obj_num_" + std::to_string(mesh) + "/mesh";
                auto meshname_ptr = const_cast<char*>(meshname.c_str());
                int meshtype = DB_UCDMESH;
                std::string& mesh_name = ucd_cross_mesh_names_per_proc[proc][mesh];
                DBPutMultimesh(dbfile, mesh_name.c_str(), 1, &meshname_ptr, &meshtype, optlist);

                if (DBMkDir(dbfile, mesh_name.c_str()) == -1)
                    TBOX_ERROR(d_object_name << "::writePlotData()\n"
                                             << "   Could not create directory named " << mesh_name << std::endl);
            }
        }

        DBFreeOptlist(optlist);
        DBClose(dbfile);

        // Create or update the dumps file on the root MPI process.
        static bool summary_file_opened = false;
        std::string path = d_dump_directory_name + "/" + VISIT_DUMPS_FILENAME;
        std::snprintf(temp_buf, sizeof(temp_buf), "%06d", d_time_step_number);
        std::string file =
            current_dump_directory_name + "/" + SILO_SUMMARY_FILE_PREFIX + temp_buf + SILO_SUMMARY_FILE_POSTFIX;
        if (!summary_file_opened)
        {
            summary_file_opened = true;
            std::ofstream sfile(path.c_str(), std::ios::out);
            sfile << file << std::endl;
            sfile.close();
        }
        else
        {
            std::ofstream sfile(path.c_str(), std::ios::app);
            sfile << file << std::endl;
            sfile.close();
        }
    }
    IBTK_MPI::barrier();
#else
    NULL_USE(SILO_MPI_ROOT);
    NULL_USE(SILO_MPI_TAG);
    NULL_USE(SILO_NAME_BUFSIZE);
    NULL_USE(d_time_step_number);
    NULL_USE(time_step_number);
    NULL_USE(simulation_time);
    TBOX_WARNING("LSiloDataWriter::writePlotData(): SILO is not installed; cannot write data." << std::endl);
#endif // if defined(IBTK_HAVE_SILO)
    return;
} // writePlotData

void
LSiloDataWriter::putToDatabase(Pointer<Database> db)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(db);
#endif
    db->putInteger("LAG_SILO_DATA_WRITER_VERSION", LAG_SILO_DATA_WRITER_VERSION);

    db->putInteger("d_num_objs", d_num_objs);

    for (int obj_num = 0; obj_num < d_num_objs; ++obj_num)
    {
        const std::string obj_num_string = "_" + std::to_string(obj_num);

        db->putInteger("d_nclouds" + obj_num_string, d_nclouds[obj_num]);
        if (d_nclouds[obj_num] > 0)
        {
            db->putStringArray("d_cloud_names" + obj_num_string,
                               &d_cloud_names[obj_num][0],
                               static_cast<int>(d_cloud_names[obj_num].size()));
            db->putIntegerArray("d_cloud_nmarks" + obj_num_string,
                                &d_cloud_nmarks[obj_num][0],
                                static_cast<int>(d_cloud_nmarks[obj_num].size()));
            db->putIntegerArray("d_cloud_first_lag_idx" + obj_num_string,
                                &d_cloud_first_lag_idx[obj_num][0],
                                static_cast<int>(d_cloud_first_lag_idx[obj_num].size()));
        }

        db->putInteger("d_nblocks" + obj_num_string, d_nblocks[obj_num]);
        if (d_nblocks[obj_num] > 0)
        {
            db->putStringArray("d_block_names" + obj_num_string,
                               &d_block_names[obj_num][0],
                               static_cast<int>(d_block_names[obj_num].size()));

            std::vector<int> flattened_block_nelems;
            flattened_block_nelems.reserve(NDIM * d_block_nelems.size());
            for (const auto& block : d_block_nelems[obj_num])
            {
                flattened_block_nelems.insert(flattened_block_nelems.end(), &block[0], &block[0] + NDIM);
            }
            db->putIntegerArray("flattened_block_nelems" + obj_num_string,
                                &flattened_block_nelems[0],
                                static_cast<int>(flattened_block_nelems.size()));

            std::vector<int> flattened_block_periodic;
            flattened_block_periodic.reserve(NDIM * d_block_periodic.size());
            for (const auto& block : d_block_periodic[obj_num])
            {
                flattened_block_periodic.insert(flattened_block_periodic.end(), &block[0], &block[0] + NDIM);
            }
            db->putIntegerArray("flattened_block_periodic" + obj_num_string,
                                &flattened_block_periodic[0],
                                static_cast<int>(flattened_block_periodic.size()));

            db->putIntegerArray("d_block_first_lag_idx" + obj_num_string,
                                &d_block_first_lag_idx[obj_num][0],
                                static_cast<int>(d_block_first_lag_idx[obj_num].size()));
        }

        db->putInteger("d_nmbs" + obj_num_string, d_nmbs[obj_num]);
        if (d_nmbs[obj_num] > 0)
        {
            db->putStringArray(
                "d_mb_names" + obj_num_string, &d_mb_names[obj_num][0], static_cast<int>(d_mb_names[obj_num].size()));

            for (int mb = 0; mb < d_nmbs[obj_num]; ++mb)
            {
                const std::string mb_string = "_" + std::to_string(mb);

                db->putInteger("d_mb_nblocks" + obj_num_string + mb_string, d_mb_nblocks[obj_num][mb]);
                if (d_mb_nblocks[obj_num][mb] > 0)
                {
                    std::vector<int> flattened_mb_nelems;
                    flattened_mb_nelems.reserve(NDIM * d_mb_nelems.size());
                    for (const auto& block : d_mb_nelems[obj_num][mb])
                    {
                        flattened_mb_nelems.insert(flattened_mb_nelems.end(), &block[0], &block[0] + NDIM);
                    }
                    db->putIntegerArray("flattened_mb_nelems" + obj_num_string + mb_string,
                                        &flattened_mb_nelems[0],
                                        static_cast<int>(flattened_mb_nelems.size()));

                    std::vector<int> flattened_mb_periodic;
                    flattened_mb_periodic.reserve(NDIM * d_mb_periodic.size());
                    for (const auto& vec : d_mb_periodic[obj_num][mb])
                    {
                        flattened_mb_periodic.insert(flattened_mb_periodic.end(), &vec[0], &vec[0] + NDIM);
                    }
                    db->putIntegerArray("flattened_mb_periodic" + obj_num_string + mb_string,
                                        &flattened_mb_periodic[0],
                                        static_cast<int>(flattened_mb_periodic.size()));

                    db->putIntegerArray("d_mb_first_lag_idx" + obj_num_string + mb_string,
                                        &d_mb_first_lag_idx[obj_num][mb][0],
                                        static_cast<int>(d_mb_first_lag_idx[obj_num][mb].size()));
                }
            }
        }

        db->putInteger("d_nucd_meshes" + obj_num_string, d_nucd_meshes[obj_num]);
        if (d_nucd_meshes[obj_num] > 0)
        {
            db->putStringArray("d_ucd_mesh_names" + obj_num_string,
                               &d_ucd_mesh_names[obj_num][0],
                               static_cast<int>(d_ucd_mesh_names[obj_num].size()));

            for (int mesh = 0; mesh < d_nucd_meshes[obj_num]; ++mesh)
            {
                const std::string mesh_string = "_" + std::to_string(mesh);

                std::vector<int> ucd_mesh_vertices_vector;
                ucd_mesh_vertices_vector.reserve(d_ucd_mesh_vertices[obj_num][mesh].size());
                for (const auto& vertex : d_ucd_mesh_vertices[obj_num][mesh])
                {
                    ucd_mesh_vertices_vector.push_back(vertex);
                }
                db->putInteger("ucd_mesh_vertices_vector.size()" + obj_num_string + mesh_string,
                               static_cast<int>(ucd_mesh_vertices_vector.size()));
                db->putIntegerArray("ucd_mesh_vertices_vector" + obj_num_string + mesh_string,
                                    &ucd_mesh_vertices_vector[0],
                                    static_cast<int>(ucd_mesh_vertices_vector.size()));

                std::vector<int> ucd_mesh_edge_maps_vector;
                ucd_mesh_edge_maps_vector.reserve(3 * d_ucd_mesh_edge_maps[obj_num][mesh].size());
                for (const auto& edge_pair : d_ucd_mesh_edge_maps[obj_num][mesh])
                {
                    const int i = edge_pair.first;
                    std::pair<int, int> e = edge_pair.second;
                    ucd_mesh_edge_maps_vector.push_back(i);
                    ucd_mesh_edge_maps_vector.push_back(e.first);
                    ucd_mesh_edge_maps_vector.push_back(e.second);
                }
                db->putInteger("ucd_mesh_edge_maps_vector.size()" + obj_num_string + mesh_string,
                               static_cast<int>(ucd_mesh_edge_maps_vector.size()));
                db->putIntegerArray("ucd_mesh_edge_maps_vector" + obj_num_string + mesh_string,
                                    &ucd_mesh_edge_maps_vector[0],
                                    static_cast<int>(ucd_mesh_edge_maps_vector.size()));
            }
        }
    }
    return;
} // putToDatabase

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

void
LSiloDataWriter::buildVecScatters(AO& ao, const int obj_num)
{
    if (!d_coords_data[obj_num]) return;

    int ierr;

    // Setup the IS data used to generate the VecScatters that redistribute the
    // distributed data into local marker clouds, local logically Cartesian
    // blocks, and local UCD meshes.
    std::vector<int> ref_is_idxs;
    for (int cloud = 0; cloud < d_nclouds[obj_num]; ++cloud)
    {
        const int nmarks = d_cloud_nmarks[obj_num][cloud];
        const int first_lag_idx = d_cloud_first_lag_idx[obj_num][cloud];
        ref_is_idxs.reserve(ref_is_idxs.size() + nmarks);

        for (int idx = first_lag_idx; idx < first_lag_idx + nmarks; ++idx)
        {
            ref_is_idxs.push_back(idx);
        }
    }

    for (int block = 0; block < d_nblocks[obj_num]; ++block)
    {
        const IntVector<NDIM>& nelem = d_block_nelems[obj_num][block];
        const int ntot = nelem.getProduct();
        const int first_lag_idx = d_block_first_lag_idx[obj_num][block];
        ref_is_idxs.reserve(ref_is_idxs.size() + ntot);

        for (int idx = first_lag_idx; idx < first_lag_idx + ntot; ++idx)
        {
            ref_is_idxs.push_back(idx);
        }
    }

    for (int mb = 0; mb < d_nmbs[obj_num]; ++mb)
    {
        for (int block = 0; block < d_mb_nblocks[obj_num][mb]; ++block)
        {
            const IntVector<NDIM>& nelem = d_mb_nelems[obj_num][mb][block];
            const int ntot = nelem.getProduct();
            const int first_lag_idx = d_mb_first_lag_idx[obj_num][mb][block];
            ref_is_idxs.reserve(ref_is_idxs.size() + ntot);

            for (int idx = first_lag_idx; idx < first_lag_idx + ntot; ++idx)
            {
                ref_is_idxs.push_back(idx);
            }
        }
    }

    for (int mesh = 0; mesh < d_nucd_meshes[obj_num]; ++mesh)
    {
        ref_is_idxs.insert(
            ref_is_idxs.end(), d_ucd_mesh_vertices[obj_num][mesh].begin(), d_ucd_mesh_vertices[obj_num][mesh].end());
    }

    // Map Lagrangian indices to PETSc indices.
    std::vector<int> ao_dummy(1, -1);
    ierr = AOApplicationToPetsc(
        ao,
        (!ref_is_idxs.empty() ? static_cast<int>(ref_is_idxs.size()) : static_cast<int>(ao_dummy.size())),
        (!ref_is_idxs.empty() ? &ref_is_idxs[0] : &ao_dummy[0]));
    IBTK_CHKERRQ(ierr);

    // Setup IS indices for all necessary data depths.
    std::map<int, std::vector<int> > src_is_idxs;

    src_is_idxs[NDIM] = ref_is_idxs;
    d_src_vec[obj_num][NDIM] = d_coords_data[obj_num]->getVec();

    for (int v = 0; v < d_nvars[obj_num]; ++v)
    {
        const int var_depth = d_var_depths[obj_num][v];
        if (src_is_idxs.find(var_depth) == src_is_idxs.end())
        {
            src_is_idxs[var_depth] = ref_is_idxs;
            d_src_vec[obj_num][var_depth] = d_var_data[obj_num][v]->getVec();
        }
    }

    // Create the VecScatters to scatter data from the global PETSc Vec to
    // contiguous local subgrids.  VecScatter objects are individually created
    // for data depths as necessary.
    for (const auto& src_is_idx : src_is_idxs)
    {
        const int depth = src_is_idx.first;
        const std::vector<int>& idxs = src_is_idx.second;
        const int idxs_sz = static_cast<int>(idxs.size());

        IS src_is;
        ierr = ISCreateBlock(
            PETSC_COMM_WORLD, depth, idxs_sz, (idxs.empty() ? nullptr : &idxs[0]), PETSC_COPY_VALUES, &src_is);
        IBTK_CHKERRQ(ierr);

        Vec& src_vec = d_src_vec[obj_num][depth];
        Vec& dst_vec = d_dst_vec[obj_num][depth];
        if (dst_vec)
        {
            ierr = VecDestroy(&dst_vec);
            IBTK_CHKERRQ(ierr);
        }
        ierr = VecCreateMPI(PETSC_COMM_WORLD, depth * idxs_sz, PETSC_DETERMINE, &dst_vec);
        IBTK_CHKERRQ(ierr);

        VecScatter& vec_scatter = d_vec_scatter[obj_num][depth];
        if (vec_scatter)
        {
            ierr = VecScatterDestroy(&vec_scatter);
            IBTK_CHKERRQ(ierr);
        }
        ierr = VecScatterCreate(src_vec, src_is, dst_vec, nullptr, &vec_scatter);
        IBTK_CHKERRQ(ierr);

        ierr = ISDestroy(&src_is);
        IBTK_CHKERRQ(ierr);
    }
    return;
} // buildVecScatters

void
LSiloDataWriter::getFromRestart()
{
    Pointer<Database> restart_db = RestartManager::getManager()->getRootDatabase();
    Pointer<Database> db;
    if (restart_db->isDatabase(d_object_name))
    {
        db = restart_db->getDatabase(d_object_name);
    }
    else
    {
        TBOX_ERROR("Restart database corresponding to " << d_object_name << " not found in restart file.");
    }

    int ver = db->getInteger("LAG_SILO_DATA_WRITER_VERSION");
    if (ver != LAG_SILO_DATA_WRITER_VERSION)
    {
        TBOX_ERROR(d_object_name << ":  "
                                 << "Restart file version different than class version.");
    }

    const int num_objs = db->getInteger("d_num_objs");
    resetNumObjs(num_objs);

    for (int obj_num = 0; obj_num < d_num_objs; ++obj_num)
    {
        const std::string obj_num_string = "_" + std::to_string(obj_num);

        d_nclouds[obj_num] = db->getInteger("d_nclouds" + obj_num_string);
        if (d_nclouds[obj_num] > 0)
        {
            d_cloud_names[obj_num].resize(d_nclouds[obj_num]);
            db->getStringArray("d_cloud_names" + obj_num_string,
                               &d_cloud_names[obj_num][0],
                               static_cast<int>(d_cloud_names[obj_num].size()));

            d_cloud_nmarks[obj_num].resize(d_nclouds[obj_num]);
            db->getIntegerArray("d_cloud_nmarks" + obj_num_string,
                                &d_cloud_nmarks[obj_num][0],
                                static_cast<int>(d_cloud_nmarks[obj_num].size()));

            d_cloud_first_lag_idx[obj_num].resize(d_nclouds[obj_num]);
            db->getIntegerArray("d_cloud_first_lag_idx" + obj_num_string,
                                &d_cloud_first_lag_idx[obj_num][0],
                                static_cast<int>(d_cloud_first_lag_idx[obj_num].size()));
        }

        d_nblocks[obj_num] = db->getInteger("d_nblocks" + obj_num_string);
        if (d_nblocks[obj_num] > 0)
        {
            d_block_names[obj_num].resize(d_nblocks[obj_num]);
            db->getStringArray("d_block_names" + obj_num_string,
                               &d_block_names[obj_num][0],
                               static_cast<int>(d_block_names[obj_num].size()));

            d_block_nelems[obj_num].resize(d_nblocks[obj_num]);
            std::vector<int> flattened_block_nelems;
            flattened_block_nelems.resize(NDIM * d_block_nelems[obj_num].size());
            db->getIntegerArray("flattened_block_nelems" + obj_num_string,
                                &flattened_block_nelems[0],
                                static_cast<int>(flattened_block_nelems.size()));
            for (unsigned int l = 0; l < d_block_nelems[obj_num].size(); ++l)
            {
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    d_block_nelems[obj_num][l](d) = flattened_block_nelems[NDIM * l + d];
                }
            }

            d_block_periodic[obj_num].resize(d_nblocks[obj_num]);
            std::vector<int> flattened_block_periodic;
            flattened_block_periodic.resize(NDIM * d_block_periodic[obj_num].size());
            db->getIntegerArray("flattened_block_periodic" + obj_num_string,
                                &flattened_block_periodic[0],
                                static_cast<int>(flattened_block_periodic.size()));
            for (unsigned int l = 0; l < d_block_periodic[obj_num].size(); ++l)
            {
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    d_block_periodic[obj_num][l](d) = flattened_block_periodic[NDIM * l + d];
                }
            }

            d_block_first_lag_idx[obj_num].resize(d_nblocks[obj_num]);
            db->getIntegerArray("d_block_first_lag_idx" + obj_num_string,
                                &d_block_first_lag_idx[obj_num][0],
                                static_cast<int>(d_block_first_lag_idx[obj_num].size()));
        }

        d_nmbs[obj_num] = db->getInteger("d_nmbs" + obj_num_string);
        if (d_nmbs[obj_num] > 0)
        {
            d_mb_names[obj_num].resize(d_nmbs[obj_num]);
            db->getStringArray(
                "d_mb_names" + obj_num_string, &d_mb_names[obj_num][0], static_cast<int>(d_mb_names[obj_num].size()));

            d_mb_nblocks.resize(d_nmbs[obj_num]);
            d_mb_nelems.resize(d_nmbs[obj_num]);
            d_mb_periodic.resize(d_nmbs[obj_num]);
            d_mb_first_lag_idx.resize(d_nmbs[obj_num]);
            for (int mb = 0; mb < d_nmbs[obj_num]; ++mb)
            {
                const std::string mb_string = "_" + std::to_string(mb);

                d_mb_nblocks[obj_num][mb] = db->getInteger("d_mb_nblocks" + obj_num_string + mb_string);
                if (d_mb_nblocks[obj_num][mb] > 0)
                {
                    d_mb_nelems[obj_num][mb].resize(d_mb_nblocks[obj_num][mb]);
                    std::vector<int> flattened_mb_nelems;
                    flattened_mb_nelems.resize(NDIM * d_mb_nelems[obj_num][mb].size());
                    db->getIntegerArray("flattened_mb_nelems" + obj_num_string + mb_string,
                                        &flattened_mb_nelems[0],
                                        static_cast<int>(flattened_mb_nelems.size()));
                    for (unsigned int l = 0; l < d_mb_nelems[obj_num][mb].size(); ++l)
                    {
                        for (unsigned int d = 0; d < NDIM; ++d)
                        {
                            d_mb_nelems[obj_num][mb][l](d) = flattened_mb_nelems[NDIM * l + d];
                        }
                    }

                    d_mb_periodic[obj_num][mb].resize(d_mb_nblocks[obj_num][mb]);
                    std::vector<int> flattened_mb_periodic;
                    flattened_mb_periodic.resize(NDIM * d_mb_periodic[obj_num][mb].size());
                    db->getIntegerArray("flattened_mb_periodic" + obj_num_string + mb_string,
                                        &flattened_mb_periodic[0],
                                        static_cast<int>(flattened_mb_periodic.size()));
                    for (unsigned int l = 0; l < d_mb_periodic[obj_num][mb].size(); ++l)
                    {
                        for (unsigned int d = 0; d < NDIM; ++d)
                        {
                            d_mb_periodic[obj_num][mb][l](d) = flattened_mb_periodic[NDIM * l + d];
                        }
                    }

                    d_mb_first_lag_idx[obj_num][mb].resize(d_mb_nblocks[obj_num][mb]);
                    db->getIntegerArray("d_mb_first_lag_idx" + obj_num_string + mb_string,
                                        &d_mb_first_lag_idx[obj_num][mb][0],
                                        static_cast<int>(d_mb_first_lag_idx[obj_num][mb].size()));
                }
            }
        }

        d_nucd_meshes[obj_num] = db->getInteger("d_nucd_meshes" + obj_num_string);
        if (d_nucd_meshes[obj_num] > 0)
        {
            d_ucd_mesh_names[obj_num].resize(d_nucd_meshes[obj_num]);
            db->getStringArray("d_ucd_mesh_names" + obj_num_string,
                               &d_ucd_mesh_names[obj_num][0],
                               static_cast<int>(d_ucd_mesh_names[obj_num].size()));

            d_ucd_mesh_vertices[obj_num].resize(d_nucd_meshes[obj_num]);
            d_ucd_mesh_edge_maps[obj_num].resize(d_nucd_meshes[obj_num]);
            for (int mesh = 0; mesh < d_nucd_meshes[obj_num]; ++mesh)
            {
                const std::string mesh_string = "_" + std::to_string(mesh);

                const int ucd_mesh_vertices_vector_size =
                    db->getInteger("ucd_mesh_vertices_vector.size()" + obj_num_string + mesh_string);
                std::vector<int> ucd_mesh_vertices_vector(ucd_mesh_vertices_vector_size);
                db->getIntegerArray("ucd_mesh_vertices_vector" + obj_num_string + mesh_string,
                                    &ucd_mesh_vertices_vector[0],
                                    static_cast<int>(ucd_mesh_vertices_vector.size()));
                d_ucd_mesh_vertices[obj_num][mesh].insert(ucd_mesh_vertices_vector.begin(),
                                                          ucd_mesh_vertices_vector.end());

                const int ucd_mesh_edge_maps_vector_size =
                    db->getInteger("ucd_mesh_edge_maps_vector.size()" + obj_num_string + mesh_string);
                std::vector<int> ucd_mesh_edge_maps_vector(ucd_mesh_edge_maps_vector_size);
                db->getIntegerArray("ucd_mesh_edge_maps_vector" + obj_num_string + mesh_string,
                                    &ucd_mesh_edge_maps_vector[0],
                                    static_cast<int>(ucd_mesh_edge_maps_vector.size()));
                for (int l = 0; l < ucd_mesh_edge_maps_vector_size / 3; ++l)
                {
                    const int idx1 = ucd_mesh_edge_maps_vector[3 * l];
                    const std::pair<int, int> e(ucd_mesh_edge_maps_vector[3 * l + 1],
                                                ucd_mesh_edge_maps_vector[3 * l + 2]);
#if !defined(NDEBUG)
                    TBOX_ASSERT(idx1 == e.first);
#endif
                    d_ucd_mesh_edge_maps[obj_num][mesh].insert(std::make_pair(idx1, e));
                }
            }
        }
    }
    return;
} // getFromRestart
} // namespace IBTK
//////////////////////////////////////////////////////////////////////////////
