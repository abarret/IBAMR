// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2018 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#ifndef included_VelocityInit
#define included_VelocityInit

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/CartGridFunction.h>
#include <ibtk/ibtk_utilities.h>

#include <libmesh/function_base.h>

#include <CartesianGridGeometry.h>
#include <RobinBcCoefStrategy.h>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

/*!
 * \brief Method to initialize the value of the advected scalar Q.
 */
class VelocityInit : public IBTK::CartGridFunction,
                     public SAMRAI::solv::RobinBcCoefStrategy<NDIM>,
                     public libMesh::FunctionBase<double>
{
public:
    /*!
     * \brief Constructor.
     */
    VelocityInit(std::string object_name, SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db, int axis);

    /*!
     * \brief Destructor.
     */
    ~VelocityInit() = default;

    /*!
     * Indicates whether the concrete CartGridFunction object is time dependent.
     */
    bool isTimeDependent() const
    {
        return true;
    }

    /*!
     * Set the data on the patch interior to the exact answer.
     */
    void setDataOnPatch(int data_idx,
                        SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > var,
                        SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                        double data_time,
                        bool initial_time = false,
                        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level =
                            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> >(NULL));

    void setBcCoefs(SAMRAI::tbox::Pointer<SAMRAI::pdat::ArrayData<NDIM, double> >& acoef_data,
                    SAMRAI::tbox::Pointer<SAMRAI::pdat::ArrayData<NDIM, double> >& bcoef_data,
                    SAMRAI::tbox::Pointer<SAMRAI::pdat::ArrayData<NDIM, double> >& gcoef_data,
                    const SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> >& variable,
                    const SAMRAI::hier::Patch<NDIM>& patch,
                    const SAMRAI::hier::BoundaryBox<NDIM>& bdry_box,
                    double fill_time = 0.0) const override;

    SAMRAI::hier::IntVector<NDIM> numberOfExtensionsFillable() const;

    /// FunctionBase functions
    std::unique_ptr<FunctionBase<double> > clone() const override;
    double operator()(const libMesh::Point& p, double time = 0.0) override;
    void operator()(const libMesh::Point& p, double time, DenseVector<double>& output) override;

protected:
private:
    double exactValue(IBTK::VectorNd pt, double t, int idx) const;
    IBTK::VectorNd exactValue(IBTK::VectorNd pt, double t) const;

    double d_theta = std::numeric_limits<double>::quiet_NaN();
    double d_ylow = std::numeric_limits<double>::quiet_NaN();
    double d_yup = std::numeric_limits<double>::quiet_NaN();
    IBTK::MatrixNd d_Q;

    double d_dpdx = std::numeric_limits<double>::quiet_NaN();
    double d_mu = std::numeric_limits<double>::quiet_NaN();
    double d_mup = std::numeric_limits<double>::quiet_NaN();
    int d_axis = -1;
};

/////////////////////////////// INLINE ///////////////////////////////////////

//#include "VelocityInit.I"

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_VelocityInit
