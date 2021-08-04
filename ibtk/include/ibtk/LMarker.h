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

#ifndef included_IBTK_LMarker
#define included_IBTK_LMarker

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "ibtk/ibtk_utilities.h"

#include "SAMRAI/hier/IntVector.h"


namespace SAMRAI
{
namespace hier
{

class Index;
} // namespace hier
namespace tbox
{
class MessageStream;
} // namespace tbox
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class LMarker provides inter-processor communications functionality
 * for a Lagrangian marker.
 */
class LMarker
{
public:
    /*!
     * \brief Default constructor.
     */
    LMarker(int idx = -1,
            const Point& X = Point::Zero(),
            const Vector& U = Vector::Zero(),
            const SAMRAI::hier::IntVector& periodic_offset = SAMRAI::hier::IntVector::getZero(SAMRAI::tbox::Dimension(NDIM)));

    /*!
     * \brief Copy constructor.
     *
     * \param from The value to copy to this object.
     */
    LMarker(const LMarker& from);

    /*!
     * \brief Constructor that unpacks data from an input stream.
     */
    LMarker(SAMRAI::tbox::MessageStream& stream, const SAMRAI::hier::IntVector& offset);

    /*!
     * \brief Destructor.
     */
    virtual ~LMarker();

    /*!
     * \brief Assignment operator.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    LMarker& operator=(const LMarker& that);

    /*!
     * \return A const reference to the marker index.
     */
    const int& getIndex() const;

    /*!
     * \return A non-const reference to the marker index.
     */
    int& getIndex();

    /*!
     * \brief Set the marker index.
     */
    void setIndex(int idx);

    /*!
     * \return A const reference to the marker position.
     */
    const Point& getPosition() const;

    /*!
     * \return A non-const reference to the marker position.
     */
    Point& getPosition();

    /*!
     * \brief Set the marker position.
     */
    void setPosition(const Point& X);

    /*!
     * \return A const reference to the marker velocity.
     */
    const Vector& getVelocity() const;

    /*!
     * \return A non-const reference to the marker velocity.
     */
    Vector& getVelocity();

    /*!
     * \brief Set the marker velocity.
     */
    void setVelocity(const Vector& U);

    /*!
     * \return A const reference to the periodic offset.
     *
     * \note If the LMarker lives in cell i, the index of the source object is
     * src_index = i - offset.
     */
    const SAMRAI::hier::IntVector& getPeriodicOffset() const;

    /*!
     * \brief Set the value of the periodic offset.
     *
     * \note If the LMarker lives in cell i, the index of the source object is
     * src_index = i - offset.
     */
    void setPeriodicOffset(const SAMRAI::hier::IntVector& offset);

    /*!
     * \brief Copy data from the source.
     *
     * \note The index of the destination object is src_index + src_offset.
     */
    void copySourceItem(const SAMRAI::hier::Index& src_index,
                        const SAMRAI::hier::IntVector& src_offset,
                        const LMarker& src_item);

    /*!
     * \brief Return an upper bound on the amount of space required to pack the
     * object to a buffer.
     */
    size_t getDataStreamSize() const;

    /*!
     * \brief Pack data into the output stream.
     */
    void packStream(SAMRAI::tbox::MessageStream& stream);

    /*!
     * \brief Unpack data from the input stream.
     */
    virtual void unpackStream(SAMRAI::tbox::MessageStream& stream, const SAMRAI::hier::IntVector& offset);

private:
    /*!
     * \brief The marker index.
     */
    int d_idx;

    /*!
     * \brief The marker position and velocity.
     */
    Point d_X;
    Vector d_U;

    /*!
     * \brief The periodic offset.
     */
    SAMRAI::hier::IntVector d_offset;
};
} // namespace IBTK

/////////////////////////////// INLINE ///////////////////////////////////////

#include "ibtk/private/LMarker-inl.h" // IWYU pragma: keep

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBTK_LMarker
