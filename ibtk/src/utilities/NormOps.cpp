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

#include "ibtk/IBTK_MPI.h"
#include "ibtk/NormOps.h"
#include "ibtk/namespaces.h" // IWYU pragma: keep

#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/math/PatchCellDataNormOpsReal.h"
#include "SAMRAI/hier/PatchData.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/math/PatchSideDataNormOpsReal.h"
#include "SAMRAI/solv/SAMRAIVectorReal.h"
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/pdat/SideVariable.h"


#include <algorithm>
#include <cmath>
#include <functional>
#include <numeric>
#include <string>
#include <vector>

namespace SAMRAI
{
namespace hier
{

class Variable;

class Box;
} // namespace hier
} // namespace SAMRAI

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// WARNING: This function will sort the input vector in ascending order.
inline double
accurate_sum(std::vector<double>& vec)
{
    if (vec.size() == 1) return vec[0];
    std::sort(vec.begin(), vec.end(), std::less<double>());
    return std::accumulate(vec.begin(), vec.end(), 0.0);
} // accurate_sum

// WARNING: This function will sort the input vector in ascending order.
inline double
accurate_sum_of_squares(std::vector<double>& vec)
{
    if (vec.size() == 1) return vec[0] * vec[0];
    std::sort(vec.begin(), vec.end(), std::less<double>());
    return std::inner_product(vec.begin(), vec.end(), vec.begin(), 0.0);
} // accurate_sum_of_squares
} // namespace

/////////////////////////////// PUBLIC ///////////////////////////////////////

double
NormOps::L1Norm(const SAMRAIVectorReal<double>* const samrai_vector, const bool local_only)
{
    const double L1_norm_local = L1Norm_local(samrai_vector);
    if (local_only) return L1_norm_local;

    const int nprocs = IBTK_MPI::getNodes();
    std::vector<double> L1_norm_proc(nprocs, 0.0);
    IBTK_MPI::allGather(L1_norm_local, &L1_norm_proc[0]);
    const double ret_val = accurate_sum(L1_norm_proc);
    return ret_val;
} // L1Norm

double
NormOps::L2Norm(const SAMRAIVectorReal<double>* const samrai_vector, const bool local_only)
{
    const double L2_norm_local = L2Norm_local(samrai_vector);
    if (local_only) return L2_norm_local;

    const int nprocs = IBTK_MPI::getNodes();
    std::vector<double> L2_norm_proc(nprocs, 0.0);
    IBTK_MPI::allGather(L2_norm_local, &L2_norm_proc[0]);
    const double ret_val = std::sqrt(accurate_sum_of_squares(L2_norm_proc));
    return ret_val;
} // L2Norm

double
NormOps::maxNorm(const SAMRAIVectorReal<double>* const samrai_vector, const bool local_only)
{
    return samrai_vector->maxNorm(local_only);
} // maxNorm

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

double
NormOps::L1Norm_local(const SAMRAIVectorReal<double>* const samrai_vector)
{
    std::vector<double> L1_norm_local_patch;
    std::shared_ptr<PatchHierarchy > hierarchy = samrai_vector->getPatchHierarchy();
    const int coarsest_ln = samrai_vector->getCoarsestLevelNumber();
    const int finest_ln = samrai_vector->getFinestLevelNumber();
    const int ncomp = samrai_vector->getNumberOfComponents();
    for (int comp = 0; comp < ncomp; ++comp)
    {
        const std::shared_ptr<Variable >& comp_var = samrai_vector->getComponentVariable(comp);
        const int comp_idx = samrai_vector->getComponentDescriptorIndex(comp);
        const int cvol_idx = samrai_vector->getControlVolumeIndex(comp);
        const bool has_cvol = cvol_idx >= 0;

        std::shared_ptr<CellVariable<double> > comp_cc_var = comp_var;
        if (comp_cc_var)
        {
            PatchCellDataNormOpsReal<double> patch_ops;
            for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
            {
                std::shared_ptr<PatchLevel > level = hierarchy->getPatchLevel(ln);
                for (PatchLevel::Iterator p = level->begin(); p != level->end(); p++)
                {
                    std::shared_ptr<Patch > patch = *p;
                    const Box& patch_box = patch->getBox();
                    std::shared_ptr<CellData<double> > comp_data = std::static_pointer_cast<CellData<double> >(patch->getPatchData(comp_idx));
                    std::shared_ptr<CellData<double> > cvol_data =
                        (has_cvol ? std::static_pointer_cast<PatchData >(patch->getPatchData(cvol_idx)) : std::shared_ptr<PatchData >(nullptr));
                    L1_norm_local_patch.push_back(patch_ops.L1Norm(comp_data, patch_box, cvol_data));
                }
            }
        }

        std::shared_ptr<SideVariable<double> > comp_sc_var = comp_var;
        if (comp_sc_var)
        {
            PatchSideDataNormOpsReal<double> patch_ops;
            for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
            {
                std::shared_ptr<PatchLevel > level = hierarchy->getPatchLevel(ln);
                for (PatchLevel::Iterator p = level->begin(); p != level->end(); p++)
                {
                    std::shared_ptr<Patch > patch = *p;
                    const Box& patch_box = patch->getBox();
                    std::shared_ptr<SideData<double> > comp_data = std::static_pointer_cast<SideData<double> >(patch->getPatchData(comp_idx));
                    std::shared_ptr<SideData<double> > cvol_data =
                        (has_cvol ? std::static_pointer_cast<PatchData >(patch->getPatchData(cvol_idx)) : std::shared_ptr<PatchData >(nullptr));
                    L1_norm_local_patch.push_back(patch_ops.L1Norm(comp_data, patch_box, cvol_data));
                }
            }
        }
    }
    return accurate_sum(L1_norm_local_patch);
} // L1Norm_local

double
NormOps::L2Norm_local(const SAMRAIVectorReal<double>* const samrai_vector)
{
    std::vector<double> L2_norm_local_patch;
    std::shared_ptr<PatchHierarchy > hierarchy = samrai_vector->getPatchHierarchy();
    const int coarsest_ln = samrai_vector->getCoarsestLevelNumber();
    const int finest_ln = samrai_vector->getFinestLevelNumber();
    const int ncomp = samrai_vector->getNumberOfComponents();
    for (int comp = 0; comp < ncomp; ++comp)
    {
        const std::shared_ptr<Variable >& comp_var = samrai_vector->getComponentVariable(comp);
        const int comp_idx = samrai_vector->getComponentDescriptorIndex(comp);
        const int cvol_idx = samrai_vector->getControlVolumeIndex(comp);
        const bool has_cvol = cvol_idx >= 0;

        std::shared_ptr<CellVariable<double> > comp_cc_var = comp_var;
        if (comp_cc_var)
        {
            PatchCellDataNormOpsReal<double> patch_ops;
            for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
            {
                std::shared_ptr<PatchLevel > level = hierarchy->getPatchLevel(ln);
                for (PatchLevel::Iterator p = level->begin(); p != level->end(); p++)
                {
                    std::shared_ptr<Patch > patch = *p;
                    const Box& patch_box = patch->getBox();
                    std::shared_ptr<CellData<double> > comp_data = std::static_pointer_cast<CellData<double> >(patch->getPatchData(comp_idx));
                    std::shared_ptr<CellData<double> > cvol_data =
                        (has_cvol ? std::static_pointer_cast<PatchData >(patch->getPatchData(cvol_idx)) : std::shared_ptr<PatchData >(nullptr));
                    L2_norm_local_patch.push_back(patch_ops.L2Norm(comp_data, patch_box, cvol_data));
                }
            }
        }

        std::shared_ptr<SideVariable<double> > comp_sc_var = comp_var;
        if (comp_sc_var)
        {
            PatchSideDataNormOpsReal<double> patch_ops;
            for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
            {
                std::shared_ptr<PatchLevel > level = hierarchy->getPatchLevel(ln);
                for (PatchLevel::Iterator p = level->begin(); p != level->end(); p++)
                {
                    std::shared_ptr<Patch > patch = *p;
                    const Box& patch_box = patch->getBox();
                    std::shared_ptr<SideData<double> > comp_data = std::static_pointer_cast<SideData<double> >(patch->getPatchData(comp_idx));
                    std::shared_ptr<SideData<double> > cvol_data =
                        (has_cvol ? std::static_pointer_cast<PatchData >(patch->getPatchData(cvol_idx)) : std::shared_ptr<PatchData >(nullptr));
                    L2_norm_local_patch.push_back(patch_ops.L2Norm(comp_data, patch_box, cvol_data));
                }
            }
        }
    }
    return std::sqrt(accurate_sum_of_squares(L2_norm_local_patch));
} // L2Norm_local

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
