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

#include "ibtk/EdgeDataSynchronization.h"
#include "ibtk/EdgeSynchCopyFillPattern.h"
#include "ibtk/namespaces.h" // IWYU pragma: keep

#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/xfer/CoarsenAlgorithm.h"
#include "SAMRAI/hier/CoarsenOperator.h"
#include "SAMRAI/xfer/CoarsenSchedule.h"
#include "SAMRAI/pdat/EdgeVariable.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/xfer/RefineAlgorithm.h"
#include "SAMRAI/hier/RefineOperator.h"
#include "SAMRAI/xfer/RefineSchedule.h"
#include "SAMRAI/hier/Variable.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/xfer/VariableFillPattern.h"

#include "SAMRAI/tbox/Utilities.h"

#include <algorithm>
#include <memory>
#include <ostream>
#include <string>
#include <vector>

namespace SAMRAI
{
namespace xfer
{

class CoarsenPatchStrategy;
} // namespace xfer
} // namespace SAMRAI

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

EdgeDataSynchronization::~EdgeDataSynchronization()
{
    if (d_is_initialized) deallocateOperatorState();
    return;
} // ~EdgeDataSynchronization

void
EdgeDataSynchronization::initializeOperatorState(const SynchronizationTransactionComponent& transaction_comp,
                                                 std::shared_ptr<PatchHierarchy > hierarchy)
{
    initializeOperatorState(std::vector<SynchronizationTransactionComponent>(1, transaction_comp), hierarchy);
    return;
} // initializeOperatorState

void
EdgeDataSynchronization::initializeOperatorState(
    const std::vector<SynchronizationTransactionComponent>& transaction_comps,
    std::shared_ptr<PatchHierarchy > hierarchy)
{
    // Deallocate the operator state if the operator is already initialized.
    if (d_is_initialized) deallocateOperatorState();

    // Reset the transaction components.
    d_transaction_comps = transaction_comps;

    // Cache hierarchy data.
    d_hierarchy = hierarchy;
    d_grid_geom = d_hierarchy->getGridGeometry();
    d_coarsest_ln = 0;
    d_finest_ln = d_hierarchy->getFinestLevelNumber();

    // Setup cached coarsen algorithms and schedules.
    VariableDatabase* var_db = VariableDatabase::getDatabase();
    bool registered_coarsen_op = false;
    d_coarsen_alg = std::make_shared<CoarsenAlgorithm>(Dimension(NDIM));
    for (const auto& transaction_comp : d_transaction_comps)
    {
        const std::string& coarsen_op_name = transaction_comp.d_coarsen_op_name;
        if (coarsen_op_name != "NONE")
        {
            const int data_idx = transaction_comp.d_data_idx;
            std::shared_ptr<Variable > var;
            var_db->mapIndexToVariable(data_idx, var);
#if !defined(NDEBUG)
            TBOX_ASSERT(var);
#endif
            std::shared_ptr<CoarsenOperator > coarsen_op = d_grid_geom->lookupCoarsenOperator(var, coarsen_op_name);
#if !defined(NDEBUG)
            TBOX_ASSERT(coarsen_op);
#endif
            d_coarsen_alg->registerCoarsen(data_idx, // destination
                                           data_idx, // source
                                           coarsen_op);
            registered_coarsen_op = true;
        }
    }

    CoarsenPatchStrategy* coarsen_strategy = nullptr;
    d_coarsen_scheds.resize(d_finest_ln + 1);
    if (registered_coarsen_op)
    {
        for (int ln = d_coarsest_ln + 1; ln <= d_finest_ln; ++ln)
        {
            std::shared_ptr<PatchLevel > level = d_hierarchy->getPatchLevel(ln);
            std::shared_ptr<PatchLevel > coarser_level = d_hierarchy->getPatchLevel(ln - 1);
            d_coarsen_scheds[ln] = d_coarsen_alg->createSchedule(coarser_level, level, coarsen_strategy);
        }
    }

    // Setup cached refine algorithms and schedules.
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        d_refine_alg[axis] = std::make_shared<RefineAlgorithm>();
        for (const auto& transaction_comp : d_transaction_comps)
        {
            const int data_idx = transaction_comp.d_data_idx;
            std::shared_ptr<Variable > var;
            var_db->mapIndexToVariable(data_idx, var);
            std::shared_ptr<EdgeVariable<double> > ec_var = var;
            if (!ec_var)
            {
                TBOX_ERROR("EdgeDataSynchronization::initializeOperatorState():\n"
                           << "  only double-precision edge-centered data is supported." << std::endl);
            }
            std::shared_ptr<RefineOperator > refine_op = nullptr;
            std::shared_ptr<VariableFillPattern > fill_pattern = std::make_shared<EdgeSynchCopyFillPattern>(axis);
            d_refine_alg[axis]->registerRefine(data_idx, // destination
                                               data_idx, // source
                                               data_idx, // temporary work space
                                               refine_op,
                                               fill_pattern);
        }

        d_refine_scheds[axis].resize(d_finest_ln + 1);
        for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
        {
            std::shared_ptr<PatchLevel > level = d_hierarchy->getPatchLevel(ln);
            d_refine_scheds[axis][ln] = d_refine_alg[axis]->createSchedule(level);
        }
    }

    // Indicate the operator is initialized.
    d_is_initialized = true;
    return;
} // initializeOperatorState

void
EdgeDataSynchronization::resetTransactionComponent(const SynchronizationTransactionComponent& transaction_comp)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_is_initialized);
#endif
    if (d_transaction_comps.size() != 1)
    {
        TBOX_ERROR("EdgeDataSynchronization::resetTransactionComponent():"
                   << "  invalid reset operation.  attempting to change the number of registered "
                      "synchronization transaction components.\n");
    }
    resetTransactionComponents(std::vector<SynchronizationTransactionComponent>(1, transaction_comp));
    return;
} // resetTransactionComponent

void
EdgeDataSynchronization::resetTransactionComponents(
    const std::vector<SynchronizationTransactionComponent>& transaction_comps)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_is_initialized);
#endif
    if (d_transaction_comps.size() != transaction_comps.size())
    {
        TBOX_ERROR("EdgeDataSynchronization::resetTransactionComponents():"
                   << "  invalid reset operation.  attempting to change the number of registered "
                      "synchronization transaction components.\n");
    }

    // Reset the transaction components.
    d_transaction_comps = transaction_comps;

    // Reset cached coarsen algorithms and schedules.
    VariableDatabase* var_db = VariableDatabase::getDatabase();
    bool registered_coarsen_op = false;
    d_coarsen_alg = std::make_shared<CoarsenAlgorithm>(Dimension(NDIM));
    for (const auto& transaction_comp : d_transaction_comps)
    {
        const std::string& coarsen_op_name = transaction_comp.d_coarsen_op_name;
        if (coarsen_op_name != "NONE")
        {
            const int data_idx = transaction_comp.d_data_idx;
            std::shared_ptr<Variable > var;
            var_db->mapIndexToVariable(data_idx, var);
#if !defined(NDEBUG)
            TBOX_ASSERT(var);
#endif
            std::shared_ptr<CoarsenOperator > coarsen_op = d_grid_geom->lookupCoarsenOperator(var, coarsen_op_name);
#if !defined(NDEBUG)
            TBOX_ASSERT(coarsen_op);
#endif
            d_coarsen_alg->registerCoarsen(data_idx, // destination
                                           data_idx, // source
                                           coarsen_op);
            registered_coarsen_op = true;
        }
    }

    if (registered_coarsen_op)
    {
        for (int ln = d_coarsest_ln + 1; ln <= d_finest_ln; ++ln)
        {
            d_coarsen_alg->resetSchedule(d_coarsen_scheds[ln]);
        }
    }

    // Reset cached refine algorithms and schedules.
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        d_refine_alg[axis] = std::make_shared<RefineAlgorithm>();
        for (const auto& transaction_comp : d_transaction_comps)
        {
            const int data_idx = transaction_comp.d_data_idx;
            std::shared_ptr<Variable > var;
            var_db->mapIndexToVariable(data_idx, var);
            std::shared_ptr<EdgeVariable<double> > ec_var = var;
            if (!ec_var)
            {
                TBOX_ERROR("EdgeDataSynchronization::resetTransactionComponents():\n"
                           << "  only double-precision edge-centered data is supported." << std::endl);
            }
            std::shared_ptr<RefineOperator > refine_op = nullptr;
            std::shared_ptr<VariableFillPattern > fill_pattern = std::make_shared<EdgeSynchCopyFillPattern>(axis);
            d_refine_alg[axis]->registerRefine(data_idx, // destination
                                               data_idx, // source
                                               data_idx, // temporary work space
                                               refine_op,
                                               fill_pattern);
        }

        for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
        {
            d_refine_alg[axis]->resetSchedule(d_refine_scheds[axis][ln]);
        }
    }
    return;
} // resetTransactionComponents

void
EdgeDataSynchronization::deallocateOperatorState()
{
    if (!d_is_initialized) return;

    // Clear cached communication schedules.
    d_coarsen_alg.setNull();
    d_coarsen_scheds.clear();

    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        d_refine_alg[axis].setNull();
        d_refine_scheds[axis].clear();
    }

    // Indicate that the operator is NOT initialized.
    d_is_initialized = false;
    return;
} // deallocateOperatorState

void
EdgeDataSynchronization::synchronizeData(const double fill_time)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_is_initialized);
#endif
    for (int ln = d_finest_ln; ln >= d_coarsest_ln; --ln)
    {
        std::shared_ptr<PatchLevel > level = d_hierarchy->getPatchLevel(ln);

        // Synchronize data on the current level.
        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            d_refine_scheds[axis][ln]->fillData(fill_time);
        }

        // When appropriate, coarsen data from the current level to the next
        // coarser level.
        if (ln > d_coarsest_ln && d_coarsen_scheds[ln]) d_coarsen_scheds[ln]->coarsenData();
    }
    return;
} // synchronizeData

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
