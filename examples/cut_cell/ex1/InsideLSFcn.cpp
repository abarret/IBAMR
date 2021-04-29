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

#include "ibamr/config.h"

#include "ibamr/cut_cells/ls_functions.h"

#include "InsideLSFcn.h"

#include <SAMRAI_config.h>

#include <array>

#include "ibamr/app_namespaces.h"

/////////////////////////////// PUBLIC ///////////////////////////////////////

InsideLSFcn::InsideLSFcn(const string& object_name, Pointer<Database> input_db) : CartGridFunction(object_name)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(!d_object_name.empty());
#endif
    d_theta = input_db->getDouble("theta");
    d_x_trans = input_db->getDouble("x_trans");
    d_y_trans = input_db->getDouble("y_trans");
    return;
} // InsideLSFcn

void
InsideLSFcn::setDataOnPatch(const int data_idx,
                            Pointer<hier::Variable<NDIM> > var,
                            Pointer<Patch<NDIM> > patch,
                            const double data_time,
                            const bool initial_time,
                            Pointer<PatchLevel<NDIM> > /*level*/)
{
    Pointer<NodeData<NDIM, double> > n_data = patch->getPatchData(data_idx);
    Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
    const double* const dx = pgeom->getDx();
    const double* const x_low = pgeom->getXLower();
    const hier::Index<NDIM>& idx_low = patch->getBox().lower();
    Matrix2d rot;
    rot(0, 0) = cos(d_theta);
    rot(0, 1) = sin(d_theta);
    rot(1, 0) = -sin(d_theta);
    rot(1, 1) = cos(d_theta);
    for (NodeIterator<NDIM> ni(patch->getBox()); ni; ni++)
    {
        const NodeIndex<NDIM>& idx = ni();

        VectorNd x_pt;
        for (int d = 0; d < NDIM; ++d) x_pt(d) = x_low[d] + dx[d] * static_cast<double>(idx(d) - idx_low(d));
        // Rotate point back
        x_pt = rot * x_pt;
        x_pt(0) -= d_x_trans;
        x_pt(1) -= d_y_trans;
        // Now calculate distances
        double y_low = -0.5, y_up = 0.5;
        double x_low = -0.5, x_up = 0.5;
        double dist_y = std::min(x_pt(1) - y_low, y_up - x_pt(1));
        double dist_x = std::min(x_pt(0) - x_low, x_up - x_pt(0));
        (*n_data)(idx) = -std::min(dist_y, dist_x);
    }
} // setDataOnPatch
//////////////////////////////////////////////////////////////////////////////
