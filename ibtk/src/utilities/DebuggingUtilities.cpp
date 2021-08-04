// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2019 by the IBAMR developers
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

#include <IBTK_config.h>

#include "ibtk/DebuggingUtilities.h"
#include "ibtk/IBTK_MPI.h"
#include "ibtk/LData.h"
#include "ibtk/namespaces.h" // IWYU pragma: keep

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/CellGeometry.h"
#include "SAMRAI/pdat/CellIndex.h"
#include "SAMRAI/pdat/FaceData.h"
#include "SAMRAI/pdat/FaceGeometry.h"
#include "SAMRAI/pdat/FaceIndex.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/pdat/NodeData.h"
#include "SAMRAI/pdat/NodeGeometry.h"
#include "SAMRAI/pdat/NodeIndex.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/pdat/SideGeometry.h"
#include "SAMRAI/pdat/SideIndex.h"
#include "SAMRAI/tbox/PIO.h"


IBTK_DISABLE_EXTRA_WARNINGS
#include <boost/multi_array.hpp>
IBTK_ENABLE_EXTRA_WARNINGS

#include <cmath>
#include <ostream>
#include <string>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

bool
DebuggingUtilities::checkCellDataForNaNs(const int patch_data_idx,
                                         const std::shared_ptr<PatchHierarchy > hierarchy,
                                         const bool interior_only,
                                         const int coarsest_ln_in,
                                         const int finest_ln_in)
{
    int num_nans = 0;
    const int coarsest_ln = coarsest_ln_in < 0 ? 0 : coarsest_ln_in;
    const int finest_ln = finest_ln_in < 0 ? hierarchy->getFinestLevelNumber() : finest_ln_in;
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        std::shared_ptr<PatchLevel > level = hierarchy->getPatchLevel(ln);
        for (PatchLevel::Iterator p = level->begin(); p != level->end(); p++)
        {
            const int patch_num = p();
            std::shared_ptr<Patch > patch = level->getPatch(patch_num);
            std::shared_ptr<CellData<double> > patch_data = std::static_pointer_cast<CellData<double> >(patch->getPatchData(patch_data_idx));
            const Box& data_box = interior_only ? patch_data->getBox() : patch_data->getGhostBox();
            for (Box::Iterator it(data_box); it; it++)
            {
                const hier::Index& i = it();
                for (int d = 0; d < patch_data->getDepth(); ++d)
                {
                    if ((*patch_data)(i, d) != (*patch_data)(i, d) || std::isnan((*patch_data)(i, d)))
                    {
                        ++num_nans;
                        plog << "found NaN!\n"
                             << "level number = " << ln << "\n"
                             << "index = " << i << "\n"
                             << "depth = " << d << "\n"
                             << "data value = " << (*patch_data)(i, d) << std::endl;
                    }
                    if (std::abs((*patch_data)(i, d)) > 1.0e12)
                    {
                        plog << "found large value!\n"
                             << "level number = " << ln << "\n"
                             << "index = " << i << "\n"
                             << "depth = " << d << "\n"
                             << "data value = " << (*patch_data)(i, d) << std::endl;
                    }
                }
            }
        }
    }
    return IBTK_MPI::minReduction(num_nans) > 0;
} // checkCellDataForNaNs

bool
DebuggingUtilities::checkFaceDataForNaNs(const int patch_data_idx,
                                         const std::shared_ptr<PatchHierarchy > hierarchy,
                                         const bool interior_only,
                                         const int coarsest_ln_in,
                                         const int finest_ln_in)
{
    int num_nans = 0;
    const int coarsest_ln = coarsest_ln_in < 0 ? 0 : coarsest_ln_in;
    const int finest_ln = finest_ln_in < 0 ? hierarchy->getFinestLevelNumber() : finest_ln_in;
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        std::shared_ptr<PatchLevel > level = hierarchy->getPatchLevel(ln);
        for (PatchLevel::Iterator p = level->begin(); p != level->end(); p++)
        {
            const int patch_num = p();
            std::shared_ptr<Patch > patch = level->getPatch(patch_num);
            std::shared_ptr<FaceData<double> > patch_data = std::static_pointer_cast<FaceData<double> >(patch->getPatchData(patch_data_idx));
            const Box& data_box = interior_only ? patch_data->getBox() : patch_data->getGhostBox();
            for (unsigned int axis = 0; axis < NDIM; ++axis)
            {
                for (Box::Iterator it(FaceGeometry::toFaceBox(data_box, axis)); it; it++)
                {
                    const hier::Index& i = it();
                    const FaceIndex i_f(i, axis, FaceIndex::Lower);
                    for (int d = 0; d < patch_data->getDepth(); ++d)
                    {
                        if ((*patch_data)(i_f, d) != (*patch_data)(i_f, d) || std::isnan((*patch_data)(i_f, d)))
                        {
                            ++num_nans;
                            plog << "found NaN!\n"
                                 << "level number = " << ln << "\n"
                                 << "index = " << i << "\n"
                                 << "depth = " << d << "\n"
                                 << "data value = " << (*patch_data)(i_f, d) << std::endl;
                        }
                        if (std::abs((*patch_data)(i_f, d)) > 1.0e12)
                        {
                            plog << "found large value!\n"
                                 << "level number = " << ln << "\n"
                                 << "index = " << i << "\n"
                                 << "depth = " << d << "\n"
                                 << "data value = " << (*patch_data)(i_f, d) << std::endl;
                        }
                    }
                }
            }
        }
    }
    return IBTK_MPI::minReduction(num_nans) > 0;
} // checkFaceDataForNaNs

bool
DebuggingUtilities::checkNodeDataForNaNs(const int patch_data_idx,
                                         const std::shared_ptr<PatchHierarchy > hierarchy,
                                         const bool interior_only,
                                         const int coarsest_ln_in,
                                         const int finest_ln_in)
{
    int num_nans = 0;
    const int coarsest_ln = coarsest_ln_in < 0 ? 0 : coarsest_ln_in;
    const int finest_ln = finest_ln_in < 0 ? hierarchy->getFinestLevelNumber() : finest_ln_in;
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        std::shared_ptr<PatchLevel > level = hierarchy->getPatchLevel(ln);
        for (PatchLevel::Iterator p = level->begin(); p != level->end(); p++)
        {
            const int patch_num = p();
            std::shared_ptr<Patch > patch = level->getPatch(patch_num);
            std::shared_ptr<NodeData<double> > patch_data = std::static_pointer_cast<NodeData<double> >(patch->getPatchData(patch_data_idx));
            const Box& data_box = interior_only ? patch_data->getBox() : patch_data->getGhostBox();
            for (Box::Iterator it(NodeGeometry::toNodeBox(data_box)); it; it++)
            {
                const hier::Index& i = it();
                const NodeIndex i_n(i, 0);
                for (int d = 0; d < patch_data->getDepth(); ++d)
                {
                    if ((*patch_data)(i_n, d) != (*patch_data)(i_n, d) || std::isnan((*patch_data)(i_n, d)))
                    {
                        ++num_nans;
                        plog << "found NaN!\n"
                             << "level number = " << ln << "\n"
                             << "index = " << i << "\n"
                             << "depth = " << d << "\n"
                             << "data value = " << (*patch_data)(i_n, d) << std::endl;
                    }
                    if (std::abs((*patch_data)(i_n, d)) > 1.0e12)
                    {
                        plog << "found large value!\n"
                             << "level number = " << ln << "\n"
                             << "index = " << i << "\n"
                             << "depth = " << d << "\n"
                             << "data value = " << (*patch_data)(i_n, d) << std::endl;
                    }
                }
            }
        }
    }
    return IBTK_MPI::minReduction(num_nans) > 0;
} // checkNodeDataForNaNs

bool
DebuggingUtilities::checkSideDataForNaNs(const int patch_data_idx,
                                         const std::shared_ptr<PatchHierarchy > hierarchy,
                                         const bool interior_only,
                                         const int coarsest_ln_in,
                                         const int finest_ln_in)
{
    int num_nans = 0;
    const int coarsest_ln = coarsest_ln_in < 0 ? 0 : coarsest_ln_in;
    const int finest_ln = finest_ln_in < 0 ? hierarchy->getFinestLevelNumber() : finest_ln_in;
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        std::shared_ptr<PatchLevel > level = hierarchy->getPatchLevel(ln);
        for (PatchLevel::Iterator p = level->begin(); p != level->end(); p++)
        {
            const int patch_num = p();
            std::shared_ptr<Patch > patch = level->getPatch(patch_num);
            std::shared_ptr<SideData<double> > patch_data = std::static_pointer_cast<SideData<double> >(patch->getPatchData(patch_data_idx));
            const Box& data_box = interior_only ? patch_data->getBox() : patch_data->getGhostBox();
            for (unsigned int axis = 0; axis < NDIM; ++axis)
            {
                for (Box::Iterator it(SideGeometry::toSideBox(data_box, axis)); it; it++)
                {
                    const hier::Index& i = it();
                    const SideIndex i_s(i, axis, SideIndex::Lower);
                    for (int d = 0; d < patch_data->getDepth(); ++d)
                    {
                        if ((*patch_data)(i_s, d) != (*patch_data)(i_s, d) || std::isnan((*patch_data)(i_s, d)))
                        {
                            ++num_nans;
                            plog << "found NaN!\n"
                                 << "level number = " << ln << "\n"
                                 << "index = " << i << "\n"
                                 << "depth = " << d << "\n"
                                 << "data value = " << (*patch_data)(i_s, d) << std::endl;
                        }
                        if (std::abs((*patch_data)(i_s, d)) > 1.0e12)
                        {
                            plog << "found large value!\n"
                                 << "level number = " << ln << "\n"
                                 << "index = " << i << "\n"
                                 << "depth = " << d << "\n"
                                 << "data value = " << (*patch_data)(i_s, d) << std::endl;
                        }
                    }
                }
            }
        }
    }
    return IBTK_MPI::minReduction(num_nans) > 0;
} // checkSideDataForNaNs

void
DebuggingUtilities::saveCellData(const int patch_data_idx,
                                 const std::shared_ptr<PatchHierarchy > hierarchy,
                                 const std::string& filename,
                                 const std::string& dirname)
{
    std::string truncated_dirname = dirname;
    while (truncated_dirname[truncated_dirname.size() - 1] == '/')
    {
        truncated_dirname = std::string(truncated_dirname, truncated_dirname.size() - 1);
    }
    Utilities::recursiveMkdir(truncated_dirname);

    const int rank = IBTK_MPI::getRank();
    const int nodes = IBTK_MPI::getNodes();
    for (int n = 0; n < nodes; ++n)
    {
        if (n == rank)
        {
            for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln)
            {
                std::shared_ptr<PatchLevel > level = hierarchy->getPatchLevel(ln);
                for (PatchLevel::Iterator p = level->begin(); p != level->end(); p++)
                {
                    const int patch_num = p();
                    std::shared_ptr<Patch > patch = level->getPatch(patch_num);
                    const Box& patch_box = patch->getBox();
                    std::shared_ptr<CellData<double> > data = std::static_pointer_cast<CellData<double> >(patch->getPatchData(patch_data_idx));

                    const std::string patch_filename = truncated_dirname + '/' + filename + '_' +
                                                       Utilities::levelToString(ln) + '_' +
                                                       Utilities::patchToString(patch_num);
                    std::ofstream of(patch_filename.c_str(), std::ios::out | std::ios::trunc | std::ios::binary);
                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        of.write(reinterpret_cast<const char*>(&patch_box.lower()(d)), sizeof(int));
                        of.write(reinterpret_cast<const char*>(&patch_box.upper()(d)), sizeof(int));
                    }
                    const int depth = data->getDepth();
                    of.write(reinterpret_cast<const char*>(&depth), sizeof(int));
                    for (int d = 0; d < depth; ++d)
                    {
                        for (Box::Iterator it(CellGeometry::toCellBox(patch_box)); it; it++)
                        {
                            const CellIndex i(it());
                            of.write(reinterpret_cast<const char*>(&(*data)(i, d)), sizeof(double));
                        }
                    }
                    of.close();
                }
            }
        }
        IBTK_MPI::barrier();
    }
    return;
} // saveCellData

void
DebuggingUtilities::saveFaceData(const int patch_data_idx,
                                 const std::shared_ptr<PatchHierarchy > hierarchy,
                                 const std::string& filename,
                                 const std::string& dirname)
{
    std::string truncated_dirname = dirname;
    while (truncated_dirname[truncated_dirname.size() - 1] == '/')
    {
        truncated_dirname = std::string(truncated_dirname, truncated_dirname.size() - 1);
    }
    Utilities::recursiveMkdir(truncated_dirname);

    const int rank = IBTK_MPI::getRank();
    const int nodes = IBTK_MPI::getNodes();
    for (int n = 0; n < nodes; ++n)
    {
        if (n == rank)
        {
            for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln)
            {
                std::shared_ptr<PatchLevel > level = hierarchy->getPatchLevel(ln);
                for (PatchLevel::Iterator p = level->begin(); p != level->end(); p++)
                {
                    const int patch_num = p();
                    std::shared_ptr<Patch > patch = level->getPatch(patch_num);
                    const Box& patch_box = patch->getBox();
                    std::shared_ptr<FaceData<double> > data = std::static_pointer_cast<FaceData<double> >(patch->getPatchData(patch_data_idx));

                    const std::string patch_filename = truncated_dirname + '/' + filename + '_' +
                                                       Utilities::levelToString(ln) + '_' +
                                                       Utilities::patchToString(patch_num);
                    std::ofstream of(patch_filename.c_str(), std::ios::out | std::ios::trunc | std::ios::binary);
                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        of.write(reinterpret_cast<const char*>(&patch_box.lower()(d)), sizeof(int));
                        of.write(reinterpret_cast<const char*>(&patch_box.upper()(d)), sizeof(int));
                    }
                    const int depth = data->getDepth();
                    of.write(reinterpret_cast<const char*>(&depth), sizeof(int));
                    for (unsigned int face = 0; face < NDIM; ++face)
                    {
                        for (int d = 0; d < depth; ++d)
                        {
                            for (Box::Iterator it(FaceGeometry::toFaceBox(patch_box, face)); it; it++)
                            {
                                const FaceIndex i(it(), face, FaceIndex::Lower);
                                of.write(reinterpret_cast<const char*>(&(*data)(i, d)), sizeof(double));
                            }
                        }
                    }
                    of.close();
                }
            }
        }
        IBTK_MPI::barrier();
    }
    return;
} // saveFaceData

void
DebuggingUtilities::saveNodeData(const int patch_data_idx,
                                 const std::shared_ptr<PatchHierarchy > hierarchy,
                                 const std::string& filename,
                                 const std::string& dirname)
{
    std::string truncated_dirname = dirname;
    while (truncated_dirname[truncated_dirname.size() - 1] == '/')
    {
        truncated_dirname = std::string(truncated_dirname, truncated_dirname.size() - 1);
    }
    Utilities::recursiveMkdir(truncated_dirname);

    const int rank = IBTK_MPI::getRank();
    const int nodes = IBTK_MPI::getNodes();
    for (int n = 0; n < nodes; ++n)
    {
        if (n == rank)
        {
            for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln)
            {
                std::shared_ptr<PatchLevel > level = hierarchy->getPatchLevel(ln);
                for (PatchLevel::Iterator p = level->begin(); p != level->end(); p++)
                {
                    const int patch_num = p();
                    std::shared_ptr<Patch > patch = level->getPatch(patch_num);
                    const Box& patch_box = patch->getBox();
                    std::shared_ptr<NodeData<double> > data = std::static_pointer_cast<NodeData<double> >(patch->getPatchData(patch_data_idx));

                    const std::string patch_filename = truncated_dirname + '/' + filename + '_' +
                                                       Utilities::levelToString(ln) + '_' +
                                                       Utilities::patchToString(patch_num);
                    std::ofstream of(patch_filename.c_str(), std::ios::out | std::ios::trunc | std::ios::binary);
                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        of.write(reinterpret_cast<const char*>(&patch_box.lower()(d)), sizeof(int));
                        of.write(reinterpret_cast<const char*>(&patch_box.upper()(d)), sizeof(int));
                    }
                    const int depth = data->getDepth();
                    of.write(reinterpret_cast<const char*>(&depth), sizeof(int));
                    for (int d = 0; d < depth; ++d)
                    {
                        for (Box::Iterator it(NodeGeometry::toNodeBox(patch_box)); it; it++)
                        {
                            const NodeIndex i(it(), IntVector(Dimension(NDIM), 0));
                            of.write(reinterpret_cast<const char*>(&(*data)(i, d)), sizeof(double));
                        }
                    }
                    of.close();
                }
            }
        }
        IBTK_MPI::barrier();
    }
    return;
} // saveNodeData

void
DebuggingUtilities::saveSideData(const int patch_data_idx,
                                 const std::shared_ptr<PatchHierarchy > hierarchy,
                                 const std::string& filename,
                                 const std::string& dirname)
{
    std::string truncated_dirname = dirname;
    while (truncated_dirname[truncated_dirname.size() - 1] == '/')
    {
        truncated_dirname = std::string(truncated_dirname, truncated_dirname.size() - 1);
    }
    Utilities::recursiveMkdir(truncated_dirname);

    const int rank = IBTK_MPI::getRank();
    const int nodes = IBTK_MPI::getNodes();
    for (int n = 0; n < nodes; ++n)
    {
        if (n == rank)
        {
            for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln)
            {
                std::shared_ptr<PatchLevel > level = hierarchy->getPatchLevel(ln);
                for (PatchLevel::Iterator p = level->begin(); p != level->end(); p++)
                {
                    const int patch_num = p();
                    std::shared_ptr<Patch > patch = level->getPatch(patch_num);
                    const Box& patch_box = patch->getBox();
                    std::shared_ptr<SideData<double> > data = std::static_pointer_cast<SideData<double> >(patch->getPatchData(patch_data_idx));

                    const std::string patch_filename = truncated_dirname + '/' + filename + '_' +
                                                       Utilities::levelToString(ln) + '_' +
                                                       Utilities::patchToString(patch_num);
                    std::ofstream of(patch_filename.c_str(), std::ios::out | std::ios::trunc | std::ios::binary);
                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        of.write(reinterpret_cast<const char*>(&patch_box.lower()(d)), sizeof(int));
                        of.write(reinterpret_cast<const char*>(&patch_box.upper()(d)), sizeof(int));
                    }
                    const int depth = data->getDepth();
                    of.write(reinterpret_cast<const char*>(&depth), sizeof(int));
                    for (unsigned int side = 0; side < NDIM; ++side)
                    {
                        for (int d = 0; d < depth; ++d)
                        {
                            for (Box::Iterator it(SideGeometry::toSideBox(patch_box, side)); it; it++)
                            {
                                const SideIndex i(it(), side, SideIndex::Lower);
                                of.write(reinterpret_cast<const char*>(&(*data)(i, d)), sizeof(double));
                            }
                        }
                    }
                    of.close();
                }
            }
        }
        IBTK_MPI::barrier();
    }
    return;
} // saveSideData

void
DebuggingUtilities::saveLagrangianData(const std::shared_ptr<LData> lag_data,
                                       const bool save_ghost_nodes,
                                       const std::string& filename,
                                       const std::string& dirname)
{
    std::string truncated_dirname = dirname;
    while (truncated_dirname[truncated_dirname.size() - 1] == '/')
    {
        truncated_dirname = std::string(truncated_dirname, truncated_dirname.size() - 1);
    }
    Utilities::recursiveMkdir(truncated_dirname);

    const boost::multi_array_ref<double, 2>& array_data = *lag_data->getGhostedLocalFormVecArray();
    const int rank = IBTK_MPI::getRank();
    const int nodes = IBTK_MPI::getNodes();
    for (int n = 0; n < nodes; ++n)
    {
        if (n == rank)
        {
            const std::string data_filename =
                truncated_dirname + '/' + filename + '_' + Utilities::processorToString(n);
            std::ofstream of(data_filename.c_str(), std::ios::out | std::ios::trunc | std::ios::binary);
            const int depth = lag_data->getDepth();
            of.write(reinterpret_cast<const char*>(&depth), sizeof(int));
            const int num_local_nodes = lag_data->getLocalNodeCount();
            of.write(reinterpret_cast<const char*>(&num_local_nodes), sizeof(int));
            for (int i = 0; i < num_local_nodes; ++i)
            {
                for (int d = 0; d < depth; ++d)
                {
                    of.write(reinterpret_cast<const char*>(&(array_data[i][d])), sizeof(double));
                }
            }
            if (save_ghost_nodes)
            {
                const int num_ghost_nodes = lag_data->getGhostNodeCount();
                of.write(reinterpret_cast<const char*>(&num_ghost_nodes), sizeof(int));
                for (int i = 0; i < num_ghost_nodes; ++i)
                {
                    for (int d = 0; d < depth; ++d)
                    {
                        of.write(reinterpret_cast<const char*>(&(array_data[i + num_local_nodes][d])), sizeof(double));
                    }
                }
            }
            of.close();
        }
        IBTK_MPI::barrier();
    }
    lag_data->restoreArrays();
    return;
} // saveLagrangianData

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
