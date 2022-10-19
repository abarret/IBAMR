// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2020 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#include <tbox/MemoryDatabase.h>

#include <libmesh/auto_ptr.h>
#include <libmesh/point.h>

#include <SAMRAI_config.h>

#include <array>

#include <ibamr/app_namespaces.h>

// Local includes
#include "ExtraStressExact.h"

/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

ExtraStressExact::ExtraStressExact(const string object_name,
                                   Pointer<Database> input_db,
                                   const int depth,
                                   const bool use_inside_solution)
    : CartGridFunction(object_name), d_depth(depth), d_use_inside_solution(use_inside_solution)
{
    d_dpdx = input_db->getDouble("dpdx");
    d_mu = input_db->getDouble("mu");
    d_mup = input_db->getDouble("mup");
    d_la = input_db->getDouble("lambda");
    d_theta = input_db->getDouble("theta");
    d_ylow = input_db->getDouble("ylow");
    d_yup = input_db->getDouble("yup");
    d_Q(0, 0) = cos(d_theta);
    d_Q(0, 1) = -sin(d_theta);
    d_Q(1, 0) = sin(d_theta);
    d_Q(1, 1) = cos(d_theta);
    return;
} // ExtraStressExact

void
ExtraStressExact::setDataOnPatch(const int data_idx,
                                 Pointer<Variable<NDIM> > /*var*/,
                                 Pointer<Patch<NDIM> > patch,
                                 const double data_time,
                                 const bool initial_time,
                                 Pointer<PatchLevel<NDIM> > /*level*/)
{
    Pointer<CellData<NDIM, double> > Q_data = patch->getPatchData(data_idx);
#if !defined(NDEBUG)
    TBOX_ASSERT(Q_data);
#endif
    if (initial_time)
    {
        MatrixNd sig = MatrixNd::Identity();
        sig = d_Q * sig * d_Q.transpose();
        Q_data->fill(sig(0, 0), 0);
        Q_data->fill(sig(1, 1), 1);
        Q_data->fill(sig(0, 1), 2);
        return;
    }
    const Box<NDIM>& patch_box = patch->getBox();
    const hier::Index<NDIM>& idx_low = patch_box.lower();
    Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();

    const double* const xlow = pgeom->getXLower();
    const double* const dx = pgeom->getDx();

    for (CellIterator<NDIM> ci(patch_box); ci; ci++)
    {
        const CellIndex<NDIM>& idx = ci();
        VectorNd x;
        for (int d = 0; d < NDIM; ++d) x[d] = xlow[d] + dx[d] * (static_cast<double>(idx(d) - idx_low(d)) + 0.5);
        MatrixNd sig = exactValue(x, data_time);
        // Now fill in data
        (*Q_data)(idx, 0) = sig(0, 0);
        (*Q_data)(idx, 1) = sig(1, 1);
        (*Q_data)(idx, 2) = sig(0, 1);
    }
    return;
} // setDataOnPatch

void
ExtraStressExact::setBcCoefs(Pointer<ArrayData<NDIM, double> >& acoef_data,
                             Pointer<ArrayData<NDIM, double> >& bcoef_data,
                             Pointer<ArrayData<NDIM, double> >& gcoef_data,
                             const Pointer<Variable<NDIM> >& /*var*/,
                             const Patch<NDIM>& patch,
                             const BoundaryBox<NDIM>& bdry_box,
                             const double fill_time) const
{
    TBOX_ASSERT(acoef_data);
    const int location_index = bdry_box.getLocationIndex();
    const int axis = location_index / 2;
    const int side = location_index % 2;
    const Box<NDIM>& box = acoef_data->getBox();
    const Box<NDIM>& patch_box = patch.getBox();
    const hier::Index<NDIM>& idx_low = patch_box.lower();
    Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch.getPatchGeometry();
    const double* const dx = pgeom->getDx();
    const double* const xlow = pgeom->getXLower();

    for (Box<NDIM>::Iterator bc(box); bc; bc++)
    {
        const hier::Index<NDIM>& i = bc();
        if (!acoef_data.isNull()) (*acoef_data)(i, 0) = 1.0;
        if (!bcoef_data.isNull()) (*bcoef_data)(i, 0) = 0.0;
        if (axis == 0)
        {
            // Determine physical location
            VectorNd x;
            for (int d = 0; d < NDIM; ++d)
                x[d] = xlow[d] + dx[d] * (static_cast<double>(i(d) - idx_low(d)) + (d == axis ? 0.0 : 0.5));
            if (!gcoef_data.isNull()) (*gcoef_data)(i, 0) = exactValue(x, fill_time, d_depth);
        }
        else if (axis == 1)
        {
            if (!gcoef_data.isNull()) (*gcoef_data)(i, 0) = d_depth == 2 ? 0.0 : 1.0;
        }
        else
        {
            TBOX_ERROR("Should not get here\n");
        }
    }
}

IntVector<NDIM>
ExtraStressExact::numberOfExtensionsFillable() const
{
    return 128;
}

std::unique_ptr<FunctionBase<double> >
ExtraStressExact::clone() const
{
    Pointer<Database> input_db = new MemoryDatabase("DB");
    input_db->putDouble("dpdx", d_dpdx);
    input_db->putDouble("mu", d_mu);
    input_db->putDouble("mup", d_mup);
    input_db->putDouble("lambda", d_la);
    input_db->putDouble("theta", d_theta);
    input_db->putDouble("ylow", d_ylow);
    input_db->putDouble("yup", d_yup);
    auto clone = libmesh_make_unique<ExtraStressExact>(d_object_name, input_db, d_depth, d_use_inside_solution);
    return clone;
}

double
ExtraStressExact::operator()(const libMesh::Point& p, const double time)
{
    if (d_use_inside_solution)
    {
        MatrixNd sig = MatrixNd::Identity();
        sig(0, 0) = 2.0 * std::pow(d_la * d_dpdx * d_ylow / (d_mu + d_mup), 2.0) + 1.0;
        sig(1, 1) = 1.0;
        sig(0, 1) = sig(1, 0) = d_la * d_dpdx * d_ylow / (d_mu + d_mup);
        sig = d_Q * sig * d_Q.transpose();
        const std::pair<int, int> sig_idxs = voigt_to_tensor_idx(d_depth);
        return sig(sig_idxs.first, sig_idxs.second);
    }
    else
    {
        return d_depth == 2 ? 0.0 : 1.0;
    }
}

void
ExtraStressExact::operator()(const libMesh::Point& p, const double time, DenseVector<double>& output)
{
    MatrixNd sig = MatrixNd::Identity();
    if (d_use_inside_solution)
    {
        sig(0, 0) = 2.0 * std::pow(d_la * d_dpdx * d_ylow / (d_mu + d_mup), 2.0) + 1.0;
        sig(1, 1) = 1.0;
        sig(0, 1) = sig(1, 0) = d_la * d_dpdx * d_ylow / (d_mu + d_mup);
        sig = d_Q * sig * d_Q.transpose();
    }
    output.resize(3);
    output(0) = sig(0, 0);
    output(1) = sig(1, 1);
    output(2) = sig(0, 1);
}

double
ExtraStressExact::exactValue(VectorNd pt, double t, int idx) const
{
    MatrixNd sig = exactValue(pt, t);
    std::pair<int, int> sig_idxs = voigt_to_tensor_idx(idx);
    return sig(sig_idxs.first, sig_idxs.second);
}
IBTK::MatrixNd
ExtraStressExact::exactValue(VectorNd pt, double t) const
{
    pt = d_Q.transpose() * pt;
    MatrixNd sig = MatrixNd::Identity();
    if (pt[1] < d_yup && pt[1] > d_ylow)
    {
        sig(0, 0) = 2.0 * std::pow(d_la * d_dpdx * pt[1] / (d_mu + d_mup), 2.0) + 1.0;
        sig(1, 1) = 1.0;
        sig(0, 1) = sig(1, 0) = d_la * d_dpdx * pt[1] / (d_mu + d_mup);
    }
    // Now transform back to original configuration
    sig = d_Q * sig * d_Q.transpose();
    return sig;
}
//////////////////////////////////////////////////////////////////////////////
