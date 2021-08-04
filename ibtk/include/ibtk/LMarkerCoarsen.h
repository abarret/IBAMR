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

#ifndef included_IBTK_LMarkerCoarsen
#define included_IBTK_LMarkerCoarsen

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/CoarsenOperator.h"
#include "SAMRAI/hier/IntVector.h"


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
 * \brief Class LMarkerCoarsen is a concrete SAMRAI::hier::CoarsenOperator for
 * restricting IB marker data from finer levels to coarser levels in the patch
 * hierarchy.
 */
class LMarkerCoarsen : public SAMRAI::hier::CoarsenOperator
{
public:
    /*!
     * \brief Default constructor.
     */
    LMarkerCoarsen();

    /*!
     * \brief Destructor.
     */
    ~LMarkerCoarsen() = default;

    /*!
     * \name Implementation of SAMRAI::hier::CoarsenOperator interface.
     */
    //\{

    /*!
     * Return the priority of this operator relative to other coarsening
     * operators.  The SAMRAI transfer routines guarantee that coarsening using
     * operators with lower priority will be performed before those with higher
     * priority.
     */
    int getOperatorPriority() const override;

    /*!
     * Return the stencil width associated with the coarsening operator.  The
     * SAMRAI transfer routines guarantee that the source patch will contain
     * sufficient ghost cell data surrounding the interior to satisfy the
     * stencil width requirements for each coarsening operator.
     */
    SAMRAI::hier::IntVector getStencilWidth(const SAMRAI::tbox::Dimension& dim) const override;

    /*!
     * Coarsen the source component on the fine patch to the destination
     * component on the coarse patch. The coarsening operation is performed on
     * the intersection of the destination patch and the coarse box.  The fine
     * patch is guaranteed to contain sufficient data for the stencil width of
     * the coarsening operator.
     */
    void coarsen(SAMRAI::hier::Patch& coarse,
                 const SAMRAI::hier::Patch& fine,
                 int dst_component,
                 int src_component,
                 const SAMRAI::hier::Box& coarse_box,
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
    LMarkerCoarsen(const LMarkerCoarsen& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    LMarkerCoarsen& operator=(const LMarkerCoarsen& that) = delete;

    /*!
     * The operator name.
     */
    static const std::string s_op_name;
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBTK_LMarkerCoarsen
