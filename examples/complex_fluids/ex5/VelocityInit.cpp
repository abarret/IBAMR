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

#include <SAMRAI_config.h>

#include <array>

#include <ibamr/app_namespaces.h>

// local includes
#include "VelocityInit.h"

/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

VelocityInit::VelocityInit(string object_name, Pointer<Database> input_db, const int axis)
    : CartGridFunction(std::move(object_name)), d_axis(axis)
{
    d_dpdx = input_db->getDouble("dpdx");
    d_mu = input_db->getDouble("mu");
    d_mup = input_db->getDouble("mup");
    d_theta = input_db->getDouble("theta");
    d_ylow = input_db->getDouble("ylow");
    d_yup = input_db->getDouble("yup");
    d_Q(0, 0) = cos(d_theta);
    d_Q(0, 1) = -sin(d_theta);
    d_Q(1, 0) = sin(d_theta);
    d_Q(1, 1) = cos(d_theta);
    return;
} // VelocityInit

void
VelocityInit::setDataOnPatch(const int data_idx,
                             Pointer<Variable<NDIM> > /*var*/,
                             Pointer<Patch<NDIM> > patch,
                             const double data_time,
                             const bool initial_time,
                             Pointer<PatchLevel<NDIM> > /*level*/)
{
    Pointer<SideData<NDIM, double> > u_sc_data = patch->getPatchData(data_idx);
    Pointer<CellData<NDIM, double> > u_cc_data = patch->getPatchData(data_idx);

    if (initial_time)
    {
        if (u_sc_data) u_sc_data->fillAll(0.0);
        if (u_cc_data) u_cc_data->fillAll(0.0);
        return;
    }

    const Box<NDIM>& patch_box = patch->getBox();
    const hier::Index<NDIM>& idx_low = patch_box.lower();
    Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
    const double* const xlow = pgeom->getXLower();
    const double* const dx = pgeom->getDx();

    if (u_sc_data)
    {
        for (int axis = 0; axis < NDIM; ++axis)
        {
            for (SideIterator<NDIM> si(patch_box, axis); si; si++)
            {
                const SideIndex<NDIM>& idx = si();
                VectorNd x;
                for (int d = 0; d < NDIM; ++d)
                    x[d] = xlow[d] + dx[d] * (static_cast<double>(idx(d) - idx_low(d)) + (axis == d ? 0 : 0.5));
                // Determine reference location
                VectorNd x_ref = d_Q.transpose() * x;
                VectorNd u = VectorNd::Zero();
                if (x_ref[1] < d_yup && x_ref[1] > d_ylow)
                {
                    u[0] = d_dpdx / (8.0 * (d_mu + d_mup)) * (4.0 * x_ref[1] * x_ref[1] - 1.0);
                }
                // Now transform back to original configuration
                u = d_Q * u;
                // Now fill in data
                (*u_sc_data)(idx) = u(axis);
            }
        }
    }
    if (u_cc_data)
    {
        for (CellIterator<NDIM> ci(patch_box); ci; ci++)
        {
            const CellIndex<NDIM>& idx = ci();
            VectorNd x;
            for (int d = 0; d < NDIM; ++d) x[d] = xlow[d] + dx[d] * (static_cast<double>(idx(d) - idx_low(d)) + 0.5);
            // Determine reference location
            VectorNd x_ref = d_Q.transpose() * x;
            VectorNd u = VectorNd::Zero();
            if (x_ref[1] < d_yup && x_ref[1] > d_ylow)
            {
                u[0] = d_dpdx / (8.0 * (d_mu + d_mup)) * (4.0 * x_ref[1] * x_ref[1] - 1.0);
            }
            // Now transform back to original configuration
            u = d_Q * u;
            // Now fill in data
            for (int d = 0; d < NDIM; ++d) (*u_cc_data)(idx, d) = u[d];
        }
    }

    return;
} // setDataOnPatch

void
VelocityInit::setBcCoefs(Pointer<ArrayData<NDIM, double> >& acoef_data,
                         Pointer<ArrayData<NDIM, double> >& bcoef_data,
                         Pointer<ArrayData<NDIM, double> >& gcoef_data,
                         const Pointer<Variable<NDIM> >& /*var*/,
                         const Patch<NDIM>& patch,
                         const BoundaryBox<NDIM>& bdry_box,
                         const double /*fill_time*/) const
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
            VectorNd x_ref = d_Q.transpose() * x;
            VectorNd u = VectorNd::Zero();
            if (x_ref[1] < d_yup && x_ref[1] > d_ylow)
            {
                u[0] = d_dpdx / (8.0 * (d_mu + d_mup)) * (4.0 * x_ref[1] * x_ref[1] - 1.0);
            }
            // Rotate back
            u = d_Q * u;
            if (!gcoef_data.isNull()) (*gcoef_data)(i, 0) = u[d_axis];
            //            pout << "d_axis = " << d_axis << " with index " << i << "\n";
            //            pout << "vel " << u[d_axis] << "\n";
        }
        else if (axis == 1)
        {
            if (!gcoef_data.isNull()) (*gcoef_data)(i, 0) = 0.0;
        }
        else
        {
            TBOX_ERROR("Should not get here\n");
        }
    }
}

IntVector<NDIM>
VelocityInit::numberOfExtensionsFillable() const
{
    return 128;
}

//////////////////////////////////////////////////////////////////////////////
