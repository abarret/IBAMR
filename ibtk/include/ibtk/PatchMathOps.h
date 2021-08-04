// ---------------------------------------------------------------------
//
// Copyright (c) 2011 - 2019 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#ifndef included_IBTK_PatchMathOps
#define included_IBTK_PatchMathOps

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "ibtk/CartCellRobinPhysBdryOp.h"
#include "ibtk/CartSideRobinPhysBdryOp.h"




namespace SAMRAI
{
namespace hier
{

class Patch;
} // namespace hier
namespace pdat
{
template <class TYPE>
class CellData;
template <class TYPE>
class FaceData;
template <class TYPE>
class NodeData;
template <class TYPE>
class EdgeData;
template <class TYPE>
class SideData;
} // namespace pdat
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class PatchMathOps provides functionality to perform mathematical
 * operations on \em individual SAMRAI::hier::Patch objects.
 *
 * \note Coarse-fine interface discretizations are handled in an implicit manner
 * via ghost cells.
 */
class PatchMathOps : public 
{
public:
    /*!
     * \brief Default constructor.
     */
    PatchMathOps() = default;

    /*!
     * \brief Destructor.
     */
    virtual ~PatchMathOps() = default;

    /*!
     * \name Mathematical operations.
     */
    //\{

    /*!
     * \brief Computes dst = curl src.
     *
     * Uses centered differences.
     */
    void curl(std::shared_ptr<SAMRAI::pdat::CellData<double> > dst,
              std::shared_ptr<SAMRAI::pdat::CellData<double> > src,
              std::shared_ptr<SAMRAI::hier::Patch > patch) const;

    /*!
     * \brief Computes dst = curl src.
     *
     * Uses centered differences.
     */
    void curl(std::shared_ptr<SAMRAI::pdat::CellData<double> > dst,
              std::shared_ptr<SAMRAI::pdat::FaceData<double> > src,
              std::shared_ptr<SAMRAI::hier::Patch > patch) const;

    /*!
     * \brief Computes dst = curl src.
     *
     * Uses centered differences.
     */
    void curl(std::shared_ptr<SAMRAI::pdat::FaceData<double> > dst,
              std::shared_ptr<SAMRAI::pdat::FaceData<double> > src,
              std::shared_ptr<SAMRAI::hier::Patch > patch) const;

    /*!
     * \brief Computes dst = curl src.
     *
     * Uses centered differences.
     */
    void curl(std::shared_ptr<SAMRAI::pdat::CellData<double> > dst,
              std::shared_ptr<SAMRAI::pdat::SideData<double> > src,
              std::shared_ptr<SAMRAI::hier::Patch > patch) const;

    /*!
     * \brief Computes dst = curl src.
     *
     * Uses centered differences.
     */
    void curl(std::shared_ptr<SAMRAI::pdat::SideData<double> > dst,
              std::shared_ptr<SAMRAI::pdat::SideData<double> > src,
              std::shared_ptr<SAMRAI::hier::Patch > patch) const;

    /*!
     * \brief Computes dst = curl src.
     *
     * Uses centered differences.
     */
    void curl(std::shared_ptr<SAMRAI::pdat::NodeData<double> > dst,
              std::shared_ptr<SAMRAI::pdat::SideData<double> > src,
              std::shared_ptr<SAMRAI::hier::Patch > patch) const;

    /*!
     * \brief Computes dst = curl src.
     *
     * Uses centered differences.
     */
    void curl(std::shared_ptr<SAMRAI::pdat::EdgeData<double> > dst,
              std::shared_ptr<SAMRAI::pdat::SideData<double> > src,
              std::shared_ptr<SAMRAI::hier::Patch > patch) const;

    /*!
     * \brief Computes dst = rot src.
     *
     * Uses centered differences.
     */
    void rot(std::shared_ptr<SAMRAI::pdat::SideData<double> > dst,
             std::shared_ptr<SAMRAI::pdat::NodeData<double> > src,
             std::shared_ptr<SAMRAI::hier::Patch > patch,
             CartSideRobinPhysBdryOp* bc_op = nullptr,
             double fill_time = 0.0) const;

    /*!
     * \brief Computes dst = rot src.
     *
     * Uses centered differences.
     */
    void rot(std::shared_ptr<SAMRAI::pdat::SideData<double> > dst,
             std::shared_ptr<SAMRAI::pdat::CellData<double> > src,
             std::shared_ptr<SAMRAI::hier::Patch > patch,
             CartSideRobinPhysBdryOp* bc_op = nullptr,
             double fill_time = 0.0) const;

    /*!
     * \brief Computes dst = rot src.
     *
     * Uses centered differences.
     */
    void rot(std::shared_ptr<SAMRAI::pdat::SideData<double> > dst,
             std::shared_ptr<SAMRAI::pdat::EdgeData<double> > src,
             std::shared_ptr<SAMRAI::hier::Patch > patch,
             CartSideRobinPhysBdryOp* bc_op = nullptr,
             double fill_time = 0.0) const;

    /*!
     * \brief Computes dst = rot src.
     *
     * Uses centered differences.
     */
    void rot(std::shared_ptr<SAMRAI::pdat::SideData<double> > dst,
             std::shared_ptr<SAMRAI::pdat::SideData<double> > src,
             std::shared_ptr<SAMRAI::hier::Patch > patch,
             CartSideRobinPhysBdryOp* bc_op = nullptr,
             double fill_time = 0.0) const;

    /*!
     * \brief Computes dst_l = alpha div src1 + beta src2_m.
     *
     * Uses centered differences.
     */
    void div(std::shared_ptr<SAMRAI::pdat::CellData<double> > dst,
             double alpha,
             std::shared_ptr<SAMRAI::pdat::CellData<double> > src1,
             double beta,
             std::shared_ptr<SAMRAI::pdat::CellData<double> > src2,
             std::shared_ptr<SAMRAI::hier::Patch > patch,
             int l = 0,
             int m = 0) const;

    /*!
     * \brief Computes dst_l = alpha div src1 + beta src2_m.
     *
     * Uses centered differences.
     */
    void div(std::shared_ptr<SAMRAI::pdat::CellData<double> > dst,
             double alpha,
             std::shared_ptr<SAMRAI::pdat::FaceData<double> > src1,
             double beta,
             std::shared_ptr<SAMRAI::pdat::CellData<double> > src2,
             std::shared_ptr<SAMRAI::hier::Patch > patch,
             int l = 0,
             int m = 0) const;

    /*!
     * \brief Computes dst_l = alpha div src1 + beta src2_m.
     *
     * Uses centered differences.
     */
    void div(std::shared_ptr<SAMRAI::pdat::CellData<double> > dst,
             double alpha,
             std::shared_ptr<SAMRAI::pdat::SideData<double> > src1,
             double beta,
             std::shared_ptr<SAMRAI::pdat::CellData<double> > src2,
             std::shared_ptr<SAMRAI::hier::Patch > patch,
             int l = 0,
             int m = 0) const;

    /*!
     * \brief Computes dst = alpha grad src1_l + beta src2.
     *
     * Uses centered differences.
     */
    void grad(std::shared_ptr<SAMRAI::pdat::CellData<double> > dst,
              double alpha,
              std::shared_ptr<SAMRAI::pdat::CellData<double> > src1,
              double beta,
              std::shared_ptr<SAMRAI::pdat::CellData<double> > src2,
              std::shared_ptr<SAMRAI::hier::Patch > patch,
              int l = 0) const;

    /*!
     * \brief Computes dst = alpha grad src1_l + beta src2.
     *
     * Uses centered differences.
     */
    void grad(std::shared_ptr<SAMRAI::pdat::FaceData<double> > dst,
              double alpha,
              std::shared_ptr<SAMRAI::pdat::CellData<double> > src1,
              double beta,
              std::shared_ptr<SAMRAI::pdat::FaceData<double> > src2,
              std::shared_ptr<SAMRAI::hier::Patch > patch,
              int l = 0) const;

    /*!
     * \brief Computes dst = alpha grad src1_l + beta src2.
     *
     * Uses centered differences.
     */
    void grad(std::shared_ptr<SAMRAI::pdat::SideData<double> > dst,
              double alpha,
              std::shared_ptr<SAMRAI::pdat::CellData<double> > src1,
              double beta,
              std::shared_ptr<SAMRAI::pdat::SideData<double> > src2,
              std::shared_ptr<SAMRAI::hier::Patch > patch,
              int l = 0) const;

    /*!
     * \brief Computes dst = alpha grad src1_l + beta src2.
     *
     * Uses centered differences.
     */
    void grad(std::shared_ptr<SAMRAI::pdat::FaceData<double> > dst,
              std::shared_ptr<SAMRAI::pdat::FaceData<double> > alpha,
              std::shared_ptr<SAMRAI::pdat::CellData<double> > src1,
              double beta,
              std::shared_ptr<SAMRAI::pdat::FaceData<double> > src2,
              std::shared_ptr<SAMRAI::hier::Patch > patch,
              int l = 0) const;

    /*!
     * \brief Computes dst = alpha grad src1_l + beta src2.
     *
     * Uses centered differences.
     */
    void grad(std::shared_ptr<SAMRAI::pdat::SideData<double> > dst,
              std::shared_ptr<SAMRAI::pdat::SideData<double> > alpha,
              std::shared_ptr<SAMRAI::pdat::CellData<double> > src1,
              double beta,
              std::shared_ptr<SAMRAI::pdat::SideData<double> > src2,
              std::shared_ptr<SAMRAI::hier::Patch > patch,
              int l = 0) const;

    /*!
     * \brief Computes the cell-centered vector field dst from the face-centered
     * vector field src by spatial averaging.
     */
    void interp(std::shared_ptr<SAMRAI::pdat::CellData<double> > dst,
                std::shared_ptr<SAMRAI::pdat::FaceData<double> > src,
                std::shared_ptr<SAMRAI::hier::Patch > patch) const;

    /*!
     * \brief Computes the cell-centered vector field dst from the side-centered
     * vector field src by spatial averaging.
     */
    void interp(std::shared_ptr<SAMRAI::pdat::CellData<double> > dst,
                std::shared_ptr<SAMRAI::pdat::SideData<double> > src,
                std::shared_ptr<SAMRAI::hier::Patch > patch) const;

    /*!
     * \brief Computes the face-centered vector field dst from the cell-centered
     * vector field src by spatial averaging.
     */
    void interp(std::shared_ptr<SAMRAI::pdat::FaceData<double> > dst,
                std::shared_ptr<SAMRAI::pdat::CellData<double> > src,
                std::shared_ptr<SAMRAI::hier::Patch > patch) const;

    /*!
     * \brief Computes the side-centered vector field dst from the cell-centered
     * vector field src by spatial averaging.
     */
    void interp(std::shared_ptr<SAMRAI::pdat::SideData<double> > dst,
                std::shared_ptr<SAMRAI::pdat::CellData<double> > src,
                std::shared_ptr<SAMRAI::hier::Patch > patch) const;

    /*!
     * \brief Computes the cell-centered vector field dst from the node-centered
     * vector field src by spatial averaging.
     */
    void interp(std::shared_ptr<SAMRAI::pdat::CellData<double> > dst,
                std::shared_ptr<SAMRAI::pdat::NodeData<double> > src,
                std::shared_ptr<SAMRAI::hier::Patch > patch) const;

    /*!
     * \brief Computes the cell-centered vector field dst from the edge-centered
     * vector field src by spatial averaging.
     */
    void interp(std::shared_ptr<SAMRAI::pdat::CellData<double> > dst,
                std::shared_ptr<SAMRAI::pdat::EdgeData<double> > src,
                std::shared_ptr<SAMRAI::hier::Patch > patch) const;

    /*!
     * \brief Computes the node-centered vector field dst from the cell-centered
     * vector field src by spatial averaging.
     */
    void interp(std::shared_ptr<SAMRAI::pdat::NodeData<double> > dst,
                std::shared_ptr<SAMRAI::pdat::CellData<double> > src,
                std::shared_ptr<SAMRAI::hier::Patch > patch,
                bool dst_ghost_interp) const;

    /*!
     * \brief Computes the edge-centered vector field dst from the cell-centered
     * vector field src by spatial averaging.
     */
    void interp(std::shared_ptr<SAMRAI::pdat::EdgeData<double> > dst,
                std::shared_ptr<SAMRAI::pdat::CellData<double> > src,
                std::shared_ptr<SAMRAI::hier::Patch > patch,
                bool dst_ghost_interp) const;

    /*!
     * \brief Computes the side-centered vector field dst from the cell-centered
     * vector field src by spatial harmonic averaging.
     */
    void harmonic_interp(std::shared_ptr<SAMRAI::pdat::SideData<double> > dst,
                         std::shared_ptr<SAMRAI::pdat::CellData<double> > src,
                         std::shared_ptr<SAMRAI::hier::Patch > patch) const;

    /*!
     * \brief Computes the node-centered vector field dst from the cell-centered
     * vector field src by spatial harmonic averaging.
     */
    void harmonic_interp(std::shared_ptr<SAMRAI::pdat::NodeData<double> > dst,
                         std::shared_ptr<SAMRAI::pdat::CellData<double> > src,
                         std::shared_ptr<SAMRAI::hier::Patch > patch,
                         bool dst_ghost_interp) const;

    /*!
     * \brief Computes the edge-centered vector field dst from the cell-centered
     * vector field src by spatial harmonic averaging.
     */
    void harmonic_interp(std::shared_ptr<SAMRAI::pdat::EdgeData<double> > dst,
                         std::shared_ptr<SAMRAI::pdat::CellData<double> > src,
                         std::shared_ptr<SAMRAI::hier::Patch > patch,
                         bool dst_ghost_interp) const;

    /*!
     * \brief Computes dst_l = alpha L src1_m + beta src1_m + gamma src2_n.
     *
     * Uses the standard 5 point stencil in 2D (7 point stencil in 3D).
     */
    void laplace(std::shared_ptr<SAMRAI::pdat::CellData<double> > dst,
                 double alpha,
                 double beta,
                 std::shared_ptr<SAMRAI::pdat::CellData<double> > src1,
                 double gamma,
                 std::shared_ptr<SAMRAI::pdat::CellData<double> > src2,
                 std::shared_ptr<SAMRAI::hier::Patch > patch,
                 int l = 0,
                 int m = 0,
                 int n = 0) const;

    /*!
     * \brief Computes dst_l = alpha L src1_m + beta src1_m + gamma src2_n.
     *
     * Uses the standard 5 point stencil in 2D (7 point stencil in 3D).
     */
    void laplace(std::shared_ptr<SAMRAI::pdat::SideData<double> > dst,
                 double alpha,
                 double beta,
                 std::shared_ptr<SAMRAI::pdat::SideData<double> > src1,
                 double gamma,
                 std::shared_ptr<SAMRAI::pdat::SideData<double> > src2,
                 std::shared_ptr<SAMRAI::hier::Patch > patch,
                 int l = 0,
                 int m = 0,
                 int n = 0) const;

    /*!
     * \brief Computes dst_l = div alpha grad src1_m + beta src1_m + gamma
     * src2_n.
     *
     * Uses a 9 point stencil in 2D (19 point stencil in 3D).
     */
    void laplace(std::shared_ptr<SAMRAI::pdat::CellData<double> > dst,
                 std::shared_ptr<SAMRAI::pdat::FaceData<double> > alpha,
                 double beta,
                 std::shared_ptr<SAMRAI::pdat::CellData<double> > src1,
                 double gamma,
                 std::shared_ptr<SAMRAI::pdat::CellData<double> > src2,
                 std::shared_ptr<SAMRAI::hier::Patch > patch,
                 int l = 0,
                 int m = 0,
                 int n = 0) const;

    /*!
     * \brief Computes dst_l = div alpha grad src1_m + beta src1_m + gamma
     * src2_n.
     *
     * Uses a 9 point stencil in 2D (19 point stencil in 3D).
     */
    void laplace(std::shared_ptr<SAMRAI::pdat::CellData<double> > dst,
                 std::shared_ptr<SAMRAI::pdat::SideData<double> > alpha,
                 double beta,
                 std::shared_ptr<SAMRAI::pdat::CellData<double> > src1,
                 double gamma,
                 std::shared_ptr<SAMRAI::pdat::CellData<double> > src2,
                 std::shared_ptr<SAMRAI::hier::Patch > patch,
                 int l = 0,
                 int m = 0,
                 int n = 0) const;

    /*!
     * \brief Computes dst_l = alpha div coef1 ((grad src1_m) + (grad src1_m)^T)
     * + beta coef2 src1_m + gamma src2_n.
     */
    void vc_laplace(std::shared_ptr<SAMRAI::pdat::SideData<double> > dst,
                    double alpha,
                    double beta,
                    std::shared_ptr<SAMRAI::pdat::NodeData<double> > coef1,
                    std::shared_ptr<SAMRAI::pdat::SideData<double> > coef2,
                    std::shared_ptr<SAMRAI::pdat::SideData<double> > src1,
                    double gamma,
                    std::shared_ptr<SAMRAI::pdat::SideData<double> > src2,
                    std::shared_ptr<SAMRAI::hier::Patch > patch,
                    bool use_harmonic_interp,
                    int l = 0,
                    int m = 0,
                    int n = 0) const;

    /*!
     * \brief Computes dst_l = alpha div coef1 ((grad src1_m) + (grad src1_m)^T)
     * + beta coef2 src1_m + gamma src2_n.
     */
    void vc_laplace(std::shared_ptr<SAMRAI::pdat::SideData<double> > dst,
                    double alpha,
                    double beta,
                    std::shared_ptr<SAMRAI::pdat::EdgeData<double> > coef1,
                    std::shared_ptr<SAMRAI::pdat::SideData<double> > coef2,
                    std::shared_ptr<SAMRAI::pdat::SideData<double> > src1,
                    double gamma,
                    std::shared_ptr<SAMRAI::pdat::SideData<double> > src2,
                    std::shared_ptr<SAMRAI::hier::Patch > patch,
                    bool use_harmonic_interp,
                    int l = 0,
                    int m = 0,
                    int n = 0) const;

    /*!
     * \brief Compute dst_i = alpha src1_j + beta src2_k, pointwise.
     */
    void pointwiseMultiply(std::shared_ptr<SAMRAI::pdat::CellData<double> > dst,
                           double alpha,
                           std::shared_ptr<SAMRAI::pdat::CellData<double> > src1,
                           double beta,
                           std::shared_ptr<SAMRAI::pdat::CellData<double> > src2,
                           std::shared_ptr<SAMRAI::hier::Patch > patch,
                           int i = 0,
                           int j = 0,
                           int k = 0) const;

    /*!
     * \brief Compute dst_i = alpha_l src1_j + beta src2_k, pointwise.
     */
    void pointwiseMultiply(std::shared_ptr<SAMRAI::pdat::CellData<double> > dst,
                           std::shared_ptr<SAMRAI::pdat::CellData<double> > alpha,
                           std::shared_ptr<SAMRAI::pdat::CellData<double> > src1,
                           double beta,
                           std::shared_ptr<SAMRAI::pdat::CellData<double> > src2,
                           std::shared_ptr<SAMRAI::hier::Patch > patch,
                           int i = 0,
                           int j = 0,
                           int k = 0,
                           int l = 0) const;

    /*!
     * \brief Compute dst_i = alpha_l src1_j + beta_m src2_k, pointwise.
     */
    void pointwiseMultiply(std::shared_ptr<SAMRAI::pdat::CellData<double> > dst,
                           std::shared_ptr<SAMRAI::pdat::CellData<double> > alpha,
                           std::shared_ptr<SAMRAI::pdat::CellData<double> > src1,
                           std::shared_ptr<SAMRAI::pdat::CellData<double> > beta,
                           std::shared_ptr<SAMRAI::pdat::CellData<double> > src2,
                           std::shared_ptr<SAMRAI::hier::Patch > patch,
                           int i = 0,
                           int j = 0,
                           int k = 0,
                           int l = 0,
                           int m = 0) const;

    /*!
     * \brief Compute dst_i = alpha src1_j + beta src2_k, pointwise.
     */
    void pointwiseMultiply(std::shared_ptr<SAMRAI::pdat::FaceData<double> > dst,
                           double alpha,
                           std::shared_ptr<SAMRAI::pdat::FaceData<double> > src1,
                           double beta,
                           std::shared_ptr<SAMRAI::pdat::FaceData<double> > src2,
                           std::shared_ptr<SAMRAI::hier::Patch > patch,
                           int i = 0,
                           int j = 0,
                           int k = 0) const;

    /*!
     * \brief Compute dst_i = alpha_l src1_j + beta src2_k, pointwise.
     */
    void pointwiseMultiply(std::shared_ptr<SAMRAI::pdat::FaceData<double> > dst,
                           std::shared_ptr<SAMRAI::pdat::FaceData<double> > alpha,
                           std::shared_ptr<SAMRAI::pdat::FaceData<double> > src1,
                           double beta,
                           std::shared_ptr<SAMRAI::pdat::FaceData<double> > src2,
                           std::shared_ptr<SAMRAI::hier::Patch > patch,
                           int i = 0,
                           int j = 0,
                           int k = 0,
                           int l = 0) const;

    /*!
     * \brief Compute dst_i = alpha_l src1_j + beta_m src2_k, pointwise.
     */
    void pointwiseMultiply(std::shared_ptr<SAMRAI::pdat::FaceData<double> > dst,
                           std::shared_ptr<SAMRAI::pdat::FaceData<double> > alpha,
                           std::shared_ptr<SAMRAI::pdat::FaceData<double> > src1,
                           std::shared_ptr<SAMRAI::pdat::FaceData<double> > beta,
                           std::shared_ptr<SAMRAI::pdat::FaceData<double> > src2,
                           std::shared_ptr<SAMRAI::hier::Patch > patch,
                           int i = 0,
                           int j = 0,
                           int k = 0,
                           int l = 0,
                           int m = 0) const;

    /*!
     * \brief Compute dst_i = alpha src1_j + beta src2_k, pointwise.
     */
    void pointwiseMultiply(std::shared_ptr<SAMRAI::pdat::NodeData<double> > dst,
                           double alpha,
                           std::shared_ptr<SAMRAI::pdat::NodeData<double> > src1,
                           double beta,
                           std::shared_ptr<SAMRAI::pdat::NodeData<double> > src2,
                           std::shared_ptr<SAMRAI::hier::Patch > patch,
                           int i = 0,
                           int j = 0,
                           int k = 0) const;

    /*!
     * \brief Compute dst_i = alpha_l src1_j + beta src2_k, pointwise.
     */
    void pointwiseMultiply(std::shared_ptr<SAMRAI::pdat::NodeData<double> > dst,
                           std::shared_ptr<SAMRAI::pdat::NodeData<double> > alpha,
                           std::shared_ptr<SAMRAI::pdat::NodeData<double> > src1,
                           double beta,
                           std::shared_ptr<SAMRAI::pdat::NodeData<double> > src2,
                           std::shared_ptr<SAMRAI::hier::Patch > patch,
                           int i = 0,
                           int j = 0,
                           int k = 0,
                           int l = 0) const;

    /*!
     * \brief Compute dst_i = alpha_l src1_j + beta_m src2_k, pointwise.
     */
    void pointwiseMultiply(std::shared_ptr<SAMRAI::pdat::NodeData<double> > dst,
                           std::shared_ptr<SAMRAI::pdat::NodeData<double> > alpha,
                           std::shared_ptr<SAMRAI::pdat::NodeData<double> > src1,
                           std::shared_ptr<SAMRAI::pdat::NodeData<double> > beta,
                           std::shared_ptr<SAMRAI::pdat::NodeData<double> > src2,
                           std::shared_ptr<SAMRAI::hier::Patch > patch,
                           int i = 0,
                           int j = 0,
                           int k = 0,
                           int l = 0,
                           int m = 0) const;

    /*!
     * \brief Compute dst_i = alpha src1_j + beta src2_k, pointwise.
     */
    void pointwiseMultiply(std::shared_ptr<SAMRAI::pdat::SideData<double> > dst,
                           double alpha,
                           std::shared_ptr<SAMRAI::pdat::SideData<double> > src1,
                           double beta,
                           std::shared_ptr<SAMRAI::pdat::SideData<double> > src2,
                           std::shared_ptr<SAMRAI::hier::Patch > patch,
                           int i = 0,
                           int j = 0,
                           int k = 0) const;

    /*!
     * \brief Compute dst_i = alpha_l src1_j + beta src2_k, pointwise.
     */
    void pointwiseMultiply(std::shared_ptr<SAMRAI::pdat::SideData<double> > dst,
                           std::shared_ptr<SAMRAI::pdat::SideData<double> > alpha,
                           std::shared_ptr<SAMRAI::pdat::SideData<double> > src1,
                           double beta,
                           std::shared_ptr<SAMRAI::pdat::SideData<double> > src2,
                           std::shared_ptr<SAMRAI::hier::Patch > patch,
                           int i = 0,
                           int j = 0,
                           int k = 0,
                           int l = 0) const;

    /*!
     * \brief Compute dst_i = alpha_l src1_j + beta_m src2_k, pointwise.
     */
    void pointwiseMultiply(std::shared_ptr<SAMRAI::pdat::SideData<double> > dst,
                           std::shared_ptr<SAMRAI::pdat::SideData<double> > alpha,
                           std::shared_ptr<SAMRAI::pdat::SideData<double> > src1,
                           std::shared_ptr<SAMRAI::pdat::SideData<double> > beta,
                           std::shared_ptr<SAMRAI::pdat::SideData<double> > src2,
                           std::shared_ptr<SAMRAI::hier::Patch > patch,
                           int i = 0,
                           int j = 0,
                           int k = 0,
                           int l = 0,
                           int m = 0) const;

    /*!
     * \brief Compute dst = |src|_1, pointwise.
     *
     * \see resetLevels
     */
    void pointwiseL1Norm(std::shared_ptr<SAMRAI::pdat::CellData<double> > dst,
                         std::shared_ptr<SAMRAI::pdat::CellData<double> > src,
                         std::shared_ptr<SAMRAI::hier::Patch > patch) const;

    /*!
     * \brief Compute dst = |src|_2, pointwise.
     *
     * \see resetLevels
     */
    void pointwiseL2Norm(std::shared_ptr<SAMRAI::pdat::CellData<double> > dst,
                         std::shared_ptr<SAMRAI::pdat::CellData<double> > src,
                         std::shared_ptr<SAMRAI::hier::Patch > patch) const;

    /*!
     * \brief Compute dst = |src|_oo, pointwise.
     *
     * \see resetLevels
     */
    void pointwiseMaxNorm(std::shared_ptr<SAMRAI::pdat::CellData<double> > dst,
                          std::shared_ptr<SAMRAI::pdat::CellData<double> > src,
                          std::shared_ptr<SAMRAI::hier::Patch > patch) const;

    /*!
     * \brief Compute dst = |src|_1, pointwise.
     *
     * \see resetLevels
     */
    void pointwiseL1Norm(std::shared_ptr<SAMRAI::pdat::NodeData<double> > dst,
                         std::shared_ptr<SAMRAI::pdat::NodeData<double> > src,
                         std::shared_ptr<SAMRAI::hier::Patch > patch) const;

    /*!
     * \brief Compute dst = |src|_2, pointwise.
     *
     * \see resetLevels
     */
    void pointwiseL2Norm(std::shared_ptr<SAMRAI::pdat::NodeData<double> > dst,
                         std::shared_ptr<SAMRAI::pdat::NodeData<double> > src,
                         std::shared_ptr<SAMRAI::hier::Patch > patch) const;

    /*!
     * \brief Compute dst = |src|_oo, pointwise.
     *
     * \see resetLevels
     */
    void pointwiseMaxNorm(std::shared_ptr<SAMRAI::pdat::NodeData<double> > dst,
                          std::shared_ptr<SAMRAI::pdat::NodeData<double> > src,
                          std::shared_ptr<SAMRAI::hier::Patch > patch) const;

    /*!
     * \brief Computes dst1 = strain src (diagonal), and dst2 = strain src (off diagonal).
     *
     * Uses centered differences.
     */
    void strain_rate(std::shared_ptr<SAMRAI::pdat::CellData<double> > dst1,
                     std::shared_ptr<SAMRAI::pdat::CellData<double> > dst2,
                     std::shared_ptr<SAMRAI::pdat::SideData<double> > src,
                     std::shared_ptr<SAMRAI::hier::Patch > patch) const;

    /*!
     * \brief Computes dst = strain src.
     *
     * Uses centered differences.
     */
    void strain_rate(std::shared_ptr<SAMRAI::pdat::CellData<double> > dst,
                     std::shared_ptr<SAMRAI::pdat::SideData<double> > src,
                     std::shared_ptr<SAMRAI::hier::Patch > patch) const;

    //\}

private:
    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    PatchMathOps(const PatchMathOps& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    PatchMathOps& operator=(const PatchMathOps& that) = delete;
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBTK_PatchMathOps
