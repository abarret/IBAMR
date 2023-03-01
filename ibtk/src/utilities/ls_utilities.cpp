#include <ibtk/IBTK_MPI.h>
#include <ibtk/ibtk_utilities.h>
#include <ibtk/ls_utilities.h>

#include <CartesianGridGeometry.h>
#include <CellData.h>
#include <CellVariable.h>
#include <CoarsenAlgorithm.h>
#include <CoarsenSchedule.h>
#include <RefineAlgorithm.h>
#include <RefineSchedule.h>
#include <SideData.h>
#include <SideVariable.h>

#include <queue>

#include <ibtk/app_namespaces.h>

namespace IBTK
{

void
flood_fill_on_level_side(const int ls_idx, Pointer<PatchLevel<NDIM> > level, const double base_val, const double time)
{
    RefineAlgorithm<NDIM> ghost_fill_alg;
    ghost_fill_alg.registerRefine(ls_idx, ls_idx, ls_idx, nullptr);
    Pointer<RefineSchedule<NDIM> > ghost_fill_sched = ghost_fill_alg.createSchedule(level);

    // Perform flood filling
    std::vector<int> patch_filled_vec(level->getNumberOfPatches());
    unsigned int patch_num = 0;
    for (PatchLevel<NDIM>::Iterator p(level); p; p++, ++patch_num)
    {
        Pointer<Patch<NDIM> > patch = level->getPatch(p());
        const Box<NDIM>& box = patch->getBox();
        Pointer<SideData<NDIM, double> > ls_data = patch->getPatchData(ls_idx);
        for (int axis = 0; axis < NDIM; ++axis)
        {
            IntVector<NDIM> side = 0;
            side(axis) = 1;
            SideData<NDIM, int> idx_touched(box, 1, ls_data->getGhostCellWidth(), side);
            idx_touched.fillAll(0);
            std::queue<SideIndex<NDIM> > idx_queue;
            bool found_pt;
            for (SideIterator<NDIM> si(box, axis); si; si++)
            {
                const SideIndex<NDIM>& idx = si();
                if ((*ls_data)(idx) <= 0.0)
                {
                    idx_queue.push(idx);
                    found_pt = true;
                }
            }
            patch_filled_vec[patch_num] = found_pt ? 1 : 0;

            while (idx_queue.size() > 0)
            {
                const SideIndex<NDIM>& idx = idx_queue.front();
                // If this point is interior and we haven't visited it yet
                if (idx_touched(idx) == 0 && ((*ls_data)(idx) == base_val || (*ls_data)(idx) <= 0.0))
                {
                    // Then insert the point into the touched list
                    idx_touched(idx) = 1;
                    // If it's unitialized, it's interior
                    if ((*ls_data)(idx) == base_val) (*ls_data)(idx) = -base_val;
                    // Now add neighboring points if they also have not been touched
                    SideIndex<NDIM> idx_s = idx + IntVector<NDIM>(0, -1);
                    SideIndex<NDIM> idx_n = idx + IntVector<NDIM>(0, 1);
                    SideIndex<NDIM> idx_e = idx + IntVector<NDIM>(1, 0);
                    SideIndex<NDIM> idx_w = idx + IntVector<NDIM>(-1, 0);
                    if (box.contains(idx_s) && ((*ls_data)(idx_s) == base_val || (*ls_data)(idx_s) < 0.0) &&
                        idx_touched(idx_s) == 0)
                        idx_queue.push(idx_s);
                    if (box.contains(idx_n) && ((*ls_data)(idx_n) == base_val || (*ls_data)(idx_n) < 0.0) &&
                        idx_touched(idx_n) == 0)
                        idx_queue.push(idx_n);
                    if (box.contains(idx_e) && ((*ls_data)(idx_e) == base_val || (*ls_data)(idx_e) < 0.0) &&
                        idx_touched(idx_e) == 0)
                        idx_queue.push(idx_e);
                    if (box.contains(idx_w) && ((*ls_data)(idx_w) == base_val || (*ls_data)(idx_w) < 0.0) &&
                        idx_touched(idx_w) == 0)
                        idx_queue.push(idx_w);
                }
                // Finished with this index
                idx_queue.pop();
            }

            int num_negative_found = 1;
            while (num_negative_found > 0)
            {
                num_negative_found = 0;
                // We fill ghost cells
                ghost_fill_sched->fillData(time);
                patch_num = 0;
                for (PatchLevel<NDIM>::Iterator p(level); p; p++)
                {
                    bool found_negative = false;
                    // If a patch has been filled, skip it.
                    if (patch_filled_vec[patch_num] == 1) continue;
                    Pointer<Patch<NDIM> > patch = level->getPatch(p());
                    // Now loop through ghost cells.
                    Pointer<CellData<NDIM, double> > ls_data = patch->getPatchData(ls_idx);
                    for (SideIterator<NDIM> si(ls_data->getGhostBox(), axis); si; si++)
                    {
                        const SideIndex<NDIM>& idx = si();
                        if (patch->getBox().contains(idx)) continue;
                        if ((*ls_data)(idx) == -base_val)
                        {
                            found_negative = true;
                            break;
                        }
                    }
                    if (found_negative)
                    {
                        ls_data->fillAll(-base_val, ls_data->getGhostBox());
                        patch_filled_vec[patch_num] = 1;
                        ++num_negative_found;
                    }
                }
                num_negative_found = IBTK_MPI::sumReduction(num_negative_found);
            }
        }
    }
}

void
flood_fill_on_level_cell(const int ls_idx, Pointer<PatchLevel<NDIM> > level, const double base_val, const double time)
{
    RefineAlgorithm<NDIM> ghost_fill_alg;
    ghost_fill_alg.registerRefine(ls_idx, ls_idx, ls_idx, nullptr);
    Pointer<RefineSchedule<NDIM> > ghost_fill_sched = ghost_fill_alg.createSchedule(level);

    // Perform flood filling
    std::vector<int> patch_filled_vec(level->getNumberOfPatches());
    unsigned int patch_num = 0;
    for (PatchLevel<NDIM>::Iterator p(level); p; p++, ++patch_num)
    {
        Pointer<Patch<NDIM> > patch = level->getPatch(p());
        const Box<NDIM>& box = patch->getBox();
        Pointer<CellData<NDIM, double> > ls_data = patch->getPatchData(ls_idx);
        CellData<NDIM, int> idx_touched(box, 1, ls_data->getGhostCellWidth());
        idx_touched.fillAll(0);
        std::queue<CellIndex<NDIM> > idx_queue;
        bool found_pt = false;
        for (CellIterator<NDIM> ci(box); ci; ci++)
        {
            const CellIndex<NDIM>& idx = ci();
            if ((*ls_data)(idx) <= 0.0)
            {
                idx_queue.push(idx);
                found_pt = true;
            }
        }
        patch_filled_vec[patch_num] = found_pt ? 1 : 0;

        // We have the boundaries found, now we loop through the queue.
        while (idx_queue.size() > 0)
        {
            const CellIndex<NDIM>& idx = idx_queue.front();
            // If this point is interior and we haven't visited it yet
            if (idx_touched(idx) == 0 && ((*ls_data)(idx) == base_val || (*ls_data)(idx) <= 0.0))
            {
                // Then insert the point into the touched list
                idx_touched(idx) = 1;
                // If it's unitialized, it's interior
                if ((*ls_data)(idx) == base_val) (*ls_data)(idx) = -base_val;
                // Now add neighboring points if they also have not been touched
                CellIndex<NDIM> idx_s = idx + IntVector<NDIM>(0, -1);
                CellIndex<NDIM> idx_n = idx + IntVector<NDIM>(0, 1);
                CellIndex<NDIM> idx_e = idx + IntVector<NDIM>(1, 0);
                CellIndex<NDIM> idx_w = idx + IntVector<NDIM>(-1, 0);
                if (box.contains(idx_s) && ((*ls_data)(idx_s) == base_val || (*ls_data)(idx_s) < 0.0) &&
                    idx_touched(idx_s) == 0)
                    idx_queue.push(idx_s);
                if (box.contains(idx_n) && ((*ls_data)(idx_n) == base_val || (*ls_data)(idx_n) < 0.0) &&
                    idx_touched(idx_n) == 0)
                    idx_queue.push(idx_n);
                if (box.contains(idx_e) && ((*ls_data)(idx_e) == base_val || (*ls_data)(idx_e) < 0.0) &&
                    idx_touched(idx_e) == 0)
                    idx_queue.push(idx_e);
                if (box.contains(idx_w) && ((*ls_data)(idx_w) == base_val || (*ls_data)(idx_w) < 0.0) &&
                    idx_touched(idx_w) == 0)
                    idx_queue.push(idx_w);
            }
            // Finished with this index
            idx_queue.pop();
        }

        // At this point, if there are any boxes that haven't been filled, they are either entirely inside or entirely
        // outside.
        int num_negative_found = 1;
        while (num_negative_found > 0)
        {
            num_negative_found = 0;
            // We fill ghost cells
            ghost_fill_sched->fillData(time);
            patch_num = 0;
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                bool found_negative = false;
                // If a patch has been filled, skip it.
                if (patch_filled_vec[patch_num] == 1) continue;
                Pointer<Patch<NDIM> > patch = level->getPatch(p());
                // Now loop through ghost cells.
                Pointer<CellData<NDIM, double> > ls_data = patch->getPatchData(ls_idx);
                for (CellIterator<NDIM> ci(ls_data->getGhostBox()); ci; ci++)
                {
                    const CellIndex<NDIM>& idx = ci();
                    if (patch->getBox().contains(idx)) continue;
                    if ((*ls_data)(idx) == -base_val)
                    {
                        found_negative = true;
                        break;
                    }
                }
                if (found_negative)
                {
                    ls_data->fillAll(-base_val, ls_data->getGhostBox());
                    patch_filled_vec[patch_num] = 1;
                    ++num_negative_found;
                }
            }
            num_negative_found = IBTK_MPI::sumReduction(num_negative_found);
        }
    }
}
} // namespace IBTK
