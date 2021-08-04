// ---------------------------------------------------------------------
//
// Copyright (c) 2019 - 2019 by the IBAMR developers
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

#include "ibtk/SAMRAIDataCache.h"
#include "ibtk/ibtk_utilities.h"
#include "ibtk/namespaces.h" // IWYU pragma: keep

#include "SAMRAI/pdat/CellDataFactory.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/pdat/EdgeDataFactory.h"
#include "SAMRAI/pdat/EdgeVariable.h"
#include "SAMRAI/pdat/FaceDataFactory.h"
#include "SAMRAI/pdat/FaceVariable.h"
#include "MultiblockDataTranslator.h"
#include "SAMRAI/pdat/NodeDataFactory.h"
#include "SAMRAI/pdat/NodeVariable.h"
#include "SAMRAI/pdat/OuteredgeDataFactory.h"
#include "SAMRAI/pdat/OuteredgeVariable.h"
#include "SAMRAI/pdat/OuterfaceDataFactory.h"
#include "SAMRAI/pdat/OuterfaceVariable.h"
#include "SAMRAI/pdat/OuternodeDataFactory.h"
#include "SAMRAI/pdat/OuternodeVariable.h"
#include "SAMRAI/pdat/OutersideDataFactory.h"
#include "SAMRAI/pdat/OutersideVariable.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/pdat/SideDataFactory.h"
#include "SAMRAI/pdat/SideVariable.h"
#include "SAMRAI/hier/Variable.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/tbox/Utilities.h"

#include <utility>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
template <typename U, typename T>
inline bool
can_convert_to(T* t)
{
    return dynamic_cast<U*>(t) != nullptr;
}

template <typename U, typename T>
inline bool
can_convert_to(const std::shared_ptr<T>& t)
{
    return dynamic_cast<U*>(t.getPointer()) != nullptr;
}

template <typename FactoryType>
std::pair</*depth*/ int, /*ghost_width*/ int>
get_characteristics(std::shared_ptr<PatchDescriptor > patch_descriptor, const int idx)
{
    std::shared_ptr<FactoryType> pdat_fac = patch_descriptor->getPatchDataFactory(idx);
    const int depth = pdat_fac->getDefaultDepth();
    const int ghost_width = pdat_fac->getGhostCellWidth().max();
#if !defined NDEBUG
    TBOX_ASSERT(ghost_width == pdat_fac->getGhostCellWidth().min());
#endif
    return std::make_pair(depth, ghost_width);
}

template <typename T>
inline bool
get_data_characteristics(const int idx, int& depth, int& ghost_width)
{
    auto var_db = VariableDatabase::getDatabase();
    std::shared_ptr<Variable > var;
    var_db->mapIndexToVariable(idx, var);
    auto patch_descriptor = var_db->getPatchDescriptor();

    // Determine the data depth and ghost width from the patch data factory.
    //
    // TODO: Make this more generic and extensible.
    bool convertable = false;
    std::pair<int, int> characteristics;
    if (can_convert_to<CellVariable<T> >(var))
    {
        convertable = true;
        characteristics = get_characteristics<CellDataFactory<T> >(patch_descriptor, idx);
    }
    else if (can_convert_to<EdgeVariable<T> >(var))
    {
        convertable = true;
        characteristics = get_characteristics<EdgeDataFactory<T> >(patch_descriptor, idx);
    }
    else if (can_convert_to<FaceVariable<T> >(var))
    {
        convertable = true;
        characteristics = get_characteristics<FaceDataFactory<T> >(patch_descriptor, idx);
    }
    else if (can_convert_to<NodeVariable<T> >(var))
    {
        convertable = true;
        characteristics = get_characteristics<NodeDataFactory<T> >(patch_descriptor, idx);
    }
    else if (can_convert_to<OuteredgeVariable<T> >(var))
    {
        convertable = true;
        characteristics = get_characteristics<OuteredgeDataFactory<T> >(patch_descriptor, idx);
    }
    else if (can_convert_to<OuterfaceVariable<T> >(var))
    {
        convertable = true;
        characteristics = get_characteristics<OuterfaceDataFactory<T> >(patch_descriptor, idx);
    }
    else if (can_convert_to<OuternodeVariable<T> >(var))
    {
        convertable = true;
        characteristics = get_characteristics<OuternodeDataFactory<T> >(patch_descriptor, idx);
    }
    else if (can_convert_to<OutersideVariable<T> >(var))
    {
        convertable = true;
        characteristics = get_characteristics<OutersideDataFactory<T> >(patch_descriptor, idx);
    }
    else if (can_convert_to<SideVariable<T> >(var))
    {
        convertable = true;
        characteristics = get_characteristics<SideDataFactory<T> >(patch_descriptor, idx);
    }
    depth = characteristics.first;
    ghost_width = characteristics.second;
    return convertable;
}
} // namespace

/////////////////////////////// PUBLIC ///////////////////////////////////////

SAMRAIDataCache::~SAMRAIDataCache()
{
    setPatchHierarchy(nullptr);
    auto var_db = VariableDatabase::getDatabase();
    for (auto cloned_idx : d_all_cloned_patch_data_idxs)
    {
        var_db->removePatchDataIndex(cloned_idx);
    }
}

void
SAMRAIDataCache::setPatchHierarchy(std::shared_ptr<PatchHierarchy > hierarchy)
{
    if (hierarchy != d_hierarchy && d_hierarchy && (d_coarsest_ln != IBTK::invalid_level_number) &&
        (d_finest_ln != IBTK::invalid_level_number))
    {
        // Clean up allocated patch data on the old hierarchy.
        for (auto cloned_idx : d_all_cloned_patch_data_idxs)
        {
            for (auto ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
            {
                std::shared_ptr<PatchLevel > level = d_hierarchy->getPatchLevel(ln);
                if (level->checkAllocated(cloned_idx)) level->deallocatePatchData(cloned_idx);
            }
        }
    }
    d_hierarchy = hierarchy;
    if (!hierarchy) resetLevels(IBTK::invalid_level_number, IBTK::invalid_level_number);
}

void
SAMRAIDataCache::resetLevels(const int coarsest_ln, const int finest_ln)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(
        (d_hierarchy &&
         ((0 <= coarsest_ln) && (coarsest_ln <= finest_ln) && (finest_ln <= d_hierarchy->getFinestLevelNumber()))) ||
        (!d_hierarchy && (coarsest_ln == IBTK::invalid_level_number) && (finest_ln == IBTK::invalid_level_number)));
#endif
    if (d_hierarchy)
    {
        // Clean up allocated patch data on the old range of levels.
        for (auto cloned_idx : d_all_cloned_patch_data_idxs)
        {
            for (auto ln = d_coarsest_ln; ln < coarsest_ln; ++ln)
            {
                std::shared_ptr<PatchLevel > level = d_hierarchy->getPatchLevel(ln);
                if (level->checkAllocated(cloned_idx)) level->deallocatePatchData(cloned_idx);
            }
            for (auto ln = finest_ln + 1; ln <= d_finest_ln; ++ln)
            {
                std::shared_ptr<PatchLevel > level = d_hierarchy->getPatchLevel(ln);
                if (level->checkAllocated(cloned_idx)) level->deallocatePatchData(cloned_idx);
            }
        }
    }
    d_coarsest_ln = coarsest_ln;
    d_finest_ln = finest_ln;
}

/////////////////////////////// PRIVATE //////////////////////////////////////

int
SAMRAIDataCache::lookupCachedPatchDataIndex(const int idx)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_hierarchy && (0 <= d_coarsest_ln) && (d_coarsest_ln <= d_finest_ln) &&
                (d_finest_ln <= d_hierarchy->getFinestLevelNumber()));
#endif
    // Find the patch data index corresponding to the data descriptor.
    int cloned_idx;
    const SAMRAIDataCache::key_type data_descriptor = construct_data_descriptor(idx);
    auto it = d_available_data_idx_map.find(data_descriptor);
    if (it == d_available_data_idx_map.end())
    {
        auto var_db = VariableDatabase::getDatabase();
        std::shared_ptr<Variable > var;
        var_db->mapIndexToVariable(idx, var);
        cloned_idx = var_db->registerClonedPatchDataIndex(var, idx);
        d_all_cloned_patch_data_idxs.insert(cloned_idx);
        it = d_available_data_idx_map.emplace(data_descriptor, idx);
    }
    else
    {
        cloned_idx = it->second;
    }

    // Remove the index from the collection of available data indices and store it in the collection of checked-out
    // indices.
    d_available_data_idx_map.erase(it);
    d_unavailable_data_idx_map.emplace(data_descriptor, cloned_idx);

    // Allocate data if needed.
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        std::shared_ptr<PatchLevel > level = d_hierarchy->getPatchLevel(ln);
        if (!level->checkAllocated(cloned_idx)) level->allocatePatchData(cloned_idx);
    }
    return cloned_idx;
}

void
SAMRAIDataCache::restoreCachedPatchDataIndex(const int cached_idx)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_hierarchy && (0 <= d_coarsest_ln) && (d_coarsest_ln <= d_finest_ln) &&
                (d_finest_ln <= d_hierarchy->getFinestLevelNumber()));
#endif
    // Find the index in the collection of checked-out indices.
    const SAMRAIDataCache::key_type data_descriptor = construct_data_descriptor(cached_idx);
    const auto data_range = d_unavailable_data_idx_map.equal_range(data_descriptor);
    auto it = data_range.first;
    for (; it != data_range.second; ++it)
    {
        if (it->second == cached_idx) break;
    }
    if (it == data_range.second)
    {
        TBOX_ERROR("could not find cached index: " << cached_idx << "\n");
    }
    d_available_data_idx_map.emplace(data_descriptor, cached_idx);
    d_unavailable_data_idx_map.erase(it);
    return;
}

SAMRAIDataCache::key_type
SAMRAIDataCache::construct_data_descriptor(const int idx)
{
    // TODO: Make this more generic and extensible.
    std::shared_ptr<Variable > var;
    auto var_db = VariableDatabase::getDatabase();
    var_db->mapIndexToVariable(idx, var);
    Variable& var_ref = *var;
    const std::type_index var_type_id = typeid(var_ref);
    int depth = 0, ghost_width = 0;
    bool convertable = get_data_characteristics<double>(idx, depth, ghost_width) ||
                       get_data_characteristics<float>(idx, depth, ghost_width) ||
                       get_data_characteristics<int>(idx, depth, ghost_width);
    if (!convertable)
    {
        TBOX_ERROR("unsupported data alignment for SAMRAIDataCache: " << var_type_id.name() << "\n");
    }
    return SAMRAIDataCache::key_type{ var_type_id, depth, ghost_width };
}

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
