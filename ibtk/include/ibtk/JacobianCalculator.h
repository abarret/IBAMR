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

#ifndef included_IBTK_JacobianCalculator
#define included_IBTK_JacobianCalculator

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <IBTK_config.h>

#include <ibtk/FECache.h>
#include <ibtk/ibtk_macros.h>
#include <ibtk/ibtk_utilities.h>

#include "SAMRAI/tbox/Utilities.h"

#include <libmesh/dense_matrix.h>
#include <libmesh/elem.h>
#include <libmesh/enum_elem_type.h>
#include <libmesh/enum_order.h>
#include <libmesh/enum_quadrature_type.h>
#include <libmesh/fe.h>
#include <libmesh/point.h>

IBTK_DISABLE_EXTRA_WARNINGS
#include <boost/multi_array.hpp>
IBTK_ENABLE_EXTRA_WARNINGS

#include <array>
#include <iosfwd>
#include <tuple>
#include <vector>

namespace IBTK
{
class Tet10Mapping;
class Hex27Mapping;
} // namespace IBTK

namespace IBTK
{
/**
 * Internal class that computes mapped quadrature point locations for
 * Lagrange-type interpolatory elements.
 *
 * @tparam dim Logical dimension of the mesh.
 *
 * @tparam spacedim Spatial dimension of the mesh (i.e., nodes have spacedim
 * meaningful coordinates).
 *
 * @tparam n_nodes Number of nodes on an element. Defaults to -1, meaning a
 * run-time calculation of the number of nodes. This template parameter is
 * useful for first-order elements since the number is small and providing it
 * improves performance.
 */
template <int dim, int spacedim = dim, int n_nodes = -1>
class PointMap
{
public:
    PointMap(const libMesh::ElemType elem_type, const std::vector<libMesh::Point>& q_points);

    /**
     * Calculate mapped quadrature points.
     */
    void getMappedQuadraturePoints(const libMesh::Point* begin,
                                   const libMesh::Point* end,
                                   std::vector<libMesh::Point>& physical_q_points);

protected:
    /**
     * Quadrature points on the reference element.
     */
    std::vector<libMesh::Point> d_reference_q_points;

    /**
     * Table containing the values of 1D shape functions (which, with a tensor
     * product, define the mapping) at reference quadrature points.
     */
    libMesh::DenseMatrix<double> d_phi;
};

class JacobianCalculator
{
public:
    /**
     * Key type. Completely describes (excepting p-refinement) a libMesh
     * quadrature rule.
     */
    using key_type = std::tuple<libMesh::ElemType, libMesh::QuadratureType, libMesh::Order>;

    /**
     * Constructor.
     */
    JacobianCalculator(const key_type quad_key);

    /**
     * Calculate the JxW values on the given element and return a reference to
     * the result.
     */
    virtual const std::vector<double>& get_JxW(const libMesh::Elem* elem) = 0;

    virtual ~JacobianCalculator() = default;

protected:
    const key_type d_quad_key;

    std::vector<libMesh::Point> d_quad_points;
    std::vector<double> d_quad_weights;
};

/*!
 * Abstract class defining the interface to a mapping.
 */
template <int dim, int spacedim = dim>
class Mapping
{
public:
    /*!
     * Recalculate relevant quantities for the provided element.
     */
    virtual void reinit(const libMesh::Elem* elem) = 0;

    /*!
     * Get the current jacobian times quadrature weight (JxW) values.
     */
    virtual const std::vector<double>& getJxW() const = 0;

    /*!
     * Get the positions of the quadrature points on the current element.
     */
    virtual const std::vector<libMesh::Point>& getQuadraturePoints() const = 0;

    /*!
     * Get the contravariants.
     */
    virtual const EigenAlignedVector<Eigen::Matrix<double, spacedim, dim> >& getContravariants() const = 0;

    /*!
     * Get the covariants.
     */
    virtual const EigenAlignedVector<Eigen::Matrix<double, spacedim, dim> >& getCovariants() const = 0;

    /*!
     * Standard 'quadrature key' alias - all the information to completely
     * define a libMesh quadrature rule.
     */
    using key_type = std::tuple<libMesh::ElemType, libMesh::QuadratureType, libMesh::Order>;

    /*!
     * Return a pointer to the correct mapping for a given quadrature key and
     * update flags object.
     */
    static std::unique_ptr<Mapping<dim, spacedim> > build(const key_type key, const FEUpdateFlags update_flags);

    virtual ~Mapping<dim, spacedim>() = default;

protected:
    /*!
     * Compute the contravariants and covariants. In general each mapping will
     * have to overload this function.
     */
    virtual void fillTransforms(const libMesh::Elem* elem) = 0;

    /*!
     * Compute determinants of contravariants (the Jacobians).
     */
    virtual void fillJacobians() = 0;

    /*!
     * Compute JxW values.
     */
    virtual void fillJxW() = 0;

    /*!
     * Compute the positions of quadrature points on the current element.
     */
    virtual void fillQuadraturePoints(const libMesh::Elem* elem) = 0;
};

/*!
 * Base class for all nodal mappings (i.e., mappings corresponding to
 * Lagrange-type finite element spaces).
 *
 * @tparam n_nodes Number of nodes of the element: defaults to runtime
 * calculation (-1).
 */
template <int dim, int spacedim = dim, int n_nodes = -1>
class NodalMapping : public Mapping<dim, spacedim>, public JacobianCalculator
{
public:
    /*!
     * Key type. Completely describes (excepting p-refinement) a libMesh
     * quadrature rule.
     */
    using key_type = std::tuple<libMesh::ElemType, libMesh::QuadratureType, libMesh::Order>;

    /*!
     * Constructor.
     */
    NodalMapping(const key_type quad_key, const FEUpdateFlags update_flags);

    /*!
     * Recalculate relevant quantities for the provided element.
     */
    virtual void reinit(const libMesh::Elem* elem) override;

    /*!
     * Calculate the JxW values on the given element and return a reference to
     * the result.
     *
     * @deprecated This function only exists for compatibility with the older
     * JacobianCalculator class.
     */
    const std::vector<double>& get_JxW(const libMesh::Elem* elem) override
    {
        reinit(elem);
        return d_JxW;
    }

    virtual const std::vector<double>& getJxW() const override
    {
#ifndef NDEBUG
        TBOX_ASSERT(d_update_flags & FEUpdateFlags::update_JxW);
#endif
        return d_JxW;
    }

    virtual const std::vector<libMesh::Point>& getQuadraturePoints() const override
    {
#ifndef NDEBUG
        TBOX_ASSERT(d_update_flags & FEUpdateFlags::update_quadrature_points);
#endif
        return d_quadrature_points;
    }

    virtual const EigenAlignedVector<Eigen::Matrix<double, spacedim, dim> >& getContravariants() const override
    {
#ifndef NDEBUG
        TBOX_ASSERT(d_update_flags & FEUpdateFlags::update_contravariants);
#endif
        return d_contravariants;
    }

    virtual const EigenAlignedVector<Eigen::Matrix<double, spacedim, dim> >& getCovariants() const override
    {
#ifndef NDEBUG
        TBOX_ASSERT(d_update_flags & FEUpdateFlags::update_covariants);
#endif
        return d_covariants;
    }

protected:
    /*!
     * Computed update flags for the mapping.
     */
    FEUpdateFlags d_update_flags;

    /*!
     * Array of contravariants.
     */
    EigenAlignedVector<Eigen::Matrix<double, spacedim, dim> > d_contravariants;

    /*!
     * Array of covariants (i.e., the transpose of the inverse of the
     * covariants when dim == spacedim)
     */
    EigenAlignedVector<Eigen::Matrix<double, spacedim, dim> > d_covariants;

    /*!
     * Array of Jacobians.
     */
    std::vector<double> d_Jacobians;

    /*!
     * Array of JxW values.
     */
    std::vector<double> d_JxW;

    /*!
     * Array of mapped quadrature points.
     */
    std::vector<libMesh::Point> d_quadrature_points;

    /*!
     * Object that computes quadrature point locations. This is sufficiently
     * different from the rest of the mapping code that it is implemented in
     * another class.
     */
    PointMap<dim, spacedim, n_nodes> d_point_map;

    /*!
     * Boolean indicating that the mapping is affine - if it is we can skip
     * some computations. The default implementation returns false. Inheriting
     * classes should overload this they represent affine mappings.
     */
    virtual bool isAffine() const;

    /*!
     * Compute determinants of contravariants (the Jacobians). The default
     * implementation given here is usually the correct one.
     */
    virtual void fillJacobians() override;

    /*!
     * Compute JxW values. The default implementation given here is usually the
     * correct one.
     */
    virtual void fillJxW() override;

    /*!
     * Compute the positions of quadrature points on the current element.
     */
    virtual void fillQuadraturePoints(const libMesh::Elem* elem) override;
};

/*!
 * A generic implementation for Lagrange-type elements: works for all elements
 * in that family but is less efficient than the specialized classes for
 * lower-order or tensor-product elements. Supports nonzero codimension.
 */
template <int dim, int spacedim = dim, int n_nodes = -1>
class LagrangeMapping : public NodalMapping<dim, spacedim, n_nodes>
{
public:
    /**
     * Key type. Completely describes (excepting p-refinement) a libMesh
     * quadrature rule.
     */
    using key_type = std::tuple<libMesh::ElemType, libMesh::QuadratureType, libMesh::Order>;

    /**
     * Constructor.
     */
    LagrangeMapping(const key_type quad_key, const FEUpdateFlags update_flags);

protected:
    virtual void fillTransforms(const libMesh::Elem* elem) override;

    /**
     * Number of nodes for the considered element type.
     */
    const int d_n_nodes;

    /**
     * Values of shape function gradients on the reference element at
     * quadrature points.
     */
    boost::multi_array<std::array<double, dim>, 2> d_dphi;

    friend class Hex27Mapping;
};

/*!
 * Specialization for TRI3 elements with codimension zero.
 */
class Tri3Mapping : public NodalMapping<2, 2, 3>
{
public:
    /**
     * Explicitly use the base class' constructor (this class does not require
     * any additional setup).
     */
    using NodalMapping<2, 2, 3>::NodalMapping;

protected:
    virtual void fillTransforms(const libMesh::Elem* elem) override;

    virtual bool isAffine() const override;
};

/*!
 * Specialization for QUAD4 elements with codimension zero.
 */
class Quad4Mapping : public NodalMapping<2, 2, 4>
{
public:
    /**
     * Explicitly use the base class' constructor (this class does not require
     * any additional setup).
     */
    using NodalMapping<2, 2, 4>::NodalMapping;

protected:
    virtual void fillTransforms(const libMesh::Elem* elem) override;
};

/*!
 * Specialization for QUAD9 elements with codimension zero.
 */
class Quad9Mapping : public NodalMapping<2, 2, 9>
{
public:
    /**
     * Key type. Completely describes (excepting p-refinement) a libMesh
     * quadrature rule.
     */
    using key_type = std::tuple<libMesh::ElemType, libMesh::QuadratureType, libMesh::Order>;

    /**
     * Constructor.
     */
    Quad9Mapping(const key_type quad_key, const FEUpdateFlags update_flags);

protected:
    virtual void fillTransforms(const libMesh::Elem* elem) override;

    /**
     * Number of 1D quadrature points in the rule used to generate the 2D
     * tensor product rule.
     */
    std::size_t d_n_oned_q_points;

    /**
     * Table containing the values of 1D shape functions (which, with a tensor
     * product, define the mapping) at reference quadrature points.
     */
    libMesh::DenseMatrix<double> d_phi;

    /**
     * Table containing the derivatives of 1D shape functions (which, with a
     * tensor product, define the mapping) at reference quadrature points.
     */
    libMesh::DenseMatrix<double> d_dphi;
};

/*!
 * Specialization for TET4 elements.
 */
class Tet4Mapping : public NodalMapping<3, 3, 4>
{
public:
    /**
     * Explicitly use the base class' constructor (this class does not require
     * any additional setup).
     */
    using NodalMapping<3, 3, 4>::NodalMapping;

protected:
    virtual void fillTransforms(const libMesh::Elem* elem) override;

    virtual bool isAffine() const override;

    friend class Tet10Mapping;
};

/*!
 * Specialization for TET10 elements. Since, for most applications and in the
 * reference configuration, most TET10 elements are actually affine this class
 * tries use the TET4 mapping whenever possible.
 */
class Tet10Mapping : public LagrangeMapping<3, 3, 10>
{
public:
    /*!
     * Key type. Completely describes (excepting p-refinement) a libMesh
     * quadrature rule.
     */
    using key_type = std::tuple<libMesh::ElemType, libMesh::QuadratureType, libMesh::Order>;

    /*!
     * Constructor.
     */
    Tet10Mapping(const key_type quad_key, const FEUpdateFlags update_flags);

    virtual void reinit(const libMesh::Elem* elem) override;

protected:
    /*!
     * TET4 mapping that is used whenever the given elem is affine.
     */
    Tet4Mapping tet4_mapping;

    /*!
     * Utility function that determines if the element is affine (i.e., all
     * nodes at edge midpoints are averages of corners)
     */
    static bool elem_is_affine(const libMesh::Elem* elem);
};

/*!
 * Specialization for HEX27 elements. Since, for most applications and in the
 * reference configuration, most HEX27 elements are actually trilinear, this
 * class tries use the lower degree mapping whenever possible.
 */
class Hex27Mapping : public LagrangeMapping<3, 3, 27>
{
public:
    /*!
     * Key type. Completely describes (excepting p-refinement) a libMesh
     * quadrature rule.
     */
    using key_type = std::tuple<libMesh::ElemType, libMesh::QuadratureType, libMesh::Order>;

    /*!
     * Constructor.
     */
    Hex27Mapping(const key_type quad_key, const FEUpdateFlags update_flags);

    virtual void reinit(const libMesh::Elem* elem) override;

protected:
    /*!
     * HEX8 mapping that is used whenever the given elem is trilinear.
     */
    LagrangeMapping<3, 3, 8> hex8_mapping;

    /*!
     * Utility function that determines if the element is trilinear (i.e., all
     * nodes at edge midpoints are averages of corners)
     */
    static bool elem_is_trilinear(const libMesh::Elem* elem);
};

// Specialization of build for 2D
template <>
std::unique_ptr<Mapping<2, 2> > Mapping<2, 2>::build(const key_type key, const FEUpdateFlags update_flags);

// Specialization of build for 3D
template <>
std::unique_ptr<Mapping<3, 3> > Mapping<3, 3>::build(const key_type key, const FEUpdateFlags update_flags);

} // namespace IBTK

#endif //#ifndef included_IBTK_JacobianCalculator
