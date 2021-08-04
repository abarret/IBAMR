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

#ifndef included_IBTK_CartCellDoubleBoundsPreservingConservativeLinearRefine
#define included_IBTK_CartCellDoubleBoundsPreservingConservativeLinearRefine

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/geom/CartesianCellDoubleConservativeLinearRefine.h"
#include "SAMRAI/pdat/CellDoubleConstantRefine.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/RefineOperator.h"


#include <string>

namespace SAMRAI
{
namespace hier
{

class Patch;

class Variable;
} // namespace hier
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class CartCellDoubleBoundsPreservingConservativeLinearRefine is a
 * concrete SAMRAI::hier::RefineOperator object which prolongs cell-centered
 * double precision patch data via conservative linear interpolation with an
 * additional bounds preservation repair step.
 */
class CartCellDoubleBoundsPreservingConservativeLinearRefine : public SAMRAI::hier::RefineOperator
{
public:
    /*!
     * \brief Default constructor.
     */
    CartCellDoubleBoundsPreservingConservativeLinearRefine() = default;

    /*!
     * \brief Destructor.
     */
    ~CartCellDoubleBoundsPreservingConservativeLinearRefine() = default;

    /*!
     * \name Implementation of SAMRAI::hier::RefineOperator interface.
     */
    //\{

    /*!
     * Return true if the refining operation matches the variable and name
     * string identifier request; false, otherwise.
     */
    bool findRefineOperator(const std::shared_ptr<SAMRAI::hier::Variable >& var,
                            const std::string& op_name) const override;

    /*!
     * Return name string identifier of the refining operation.
     */
    const std::string& getOperatorName() const override;

    /*!
     * Return the priority of this operator relative to other refining
     * operators.  The SAMRAI transfer routines guarantee that refining using
     * operators with lower priority will be performed before those with higher
     * priority.
     */
    int getOperatorPriority() const override;

    /*!
     * Return the stencil width associated with the refining operator.  The
     * SAMRAI transfer routines guarantee that the source patch will contain
     * sufficient ghost cell data surrounding the interior to satisfy the
     * stencil width requirements for each refining operator.
     */
    SAMRAI::hier::IntVector getStencilWidth() const override;

    /*!
     * Refine the source component on the fine patch to the destination
     * component on the coarse patch. The refining operation is performed on the
     * intersection of the destination patch and the coarse box.  The fine patch
     * is guaranteed to contain sufficient data for the stencil width of the
     * refining operator.
     */
    void refine(SAMRAI::hier::Patch& fine,
                const SAMRAI::hier::Patch& coarse,
                int dst_component,
                int src_component,
                const SAMRAI::hier::Box& fine_box,
                const SAMRAI::hier::IntVector& ratio) const override;

    //\}

protected:
private:
    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    CartCellDoubleBoundsPreservingConservativeLinearRefine(
        const CartCellDoubleBoundsPreservingConservativeLinearRefine& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    CartCellDoubleBoundsPreservingConservativeLinearRefine&
    operator=(const CartCellDoubleBoundsPreservingConservativeLinearRefine& that) = delete;

    /*!
     * The operator name.
     */
    static const std::string s_op_name;

    /*!
     * The basic, non-bounds preserving conservative linear refine operator.
     */
    SAMRAI::geom::CartesianCellDoubleConservativeLinearRefine d_conservative_linear_refine_op;

    /*!
     * The constant refine operator.
     */
    SAMRAI::pdat::CellDoubleConstantRefine d_constant_refine_op;
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBTK_CartCellDoubleBoundsPreservingConservativeLinearRefine
