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

#ifndef included_IBTK_LSetDataIterator
#define included_IBTK_LSetDataIterator

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "ibtk/LSet.h"

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/pdat/IndexData.h"


namespace SAMRAI
{
namespace hier
{

class Index;
} // namespace hier
namespace pdat
{

class CellGeometry;
} // namespace pdat
} // namespace SAMRAI

/////////////////////////////// FORWARD DECLARATIONS /////////////////////////

namespace IBTK
{
template <class T>
class LSetData;
} // namespace IBTK

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class LSetDataIterator is an iterator class which may be used to
 * iterate through LSet objects associated with a specified box in cell-centered
 * index space.
 */
template <class T>
class LSetDataIterator
{
public:
    friend class LSetData<T>;

    /*!
     * \brief Class constructor.
     */
    LSetDataIterator();

    /*!
     * \brief Class constructor.
     */
    LSetDataIterator(const LSetDataIterator& that);

    /*!
     * \brief Class destructor.
     */
    virtual ~LSetDataIterator();

    /*!
     * \brief Assignment operator.
     */
    LSetDataIterator<T>& operator=(const LSetDataIterator<T>& that);

    /*!
     * \brief Test two iterators for equality.
     */
    bool operator==(const LSetDataIterator<T>& that);

    /*!
     * \brief Test two iterators for inequality.
     */
    bool operator!=(const LSetDataIterator<T>& that);

    /*!
     * \brief Prefix increment operator.
     */
    LSetDataIterator<T>& operator++();

    /*!
     * \brief Postfix increment operator.
     */
    LSetDataIterator<T> operator++(int);

    /*!
     * \brief Return a reference to the Lagrangian data item referred to by the
     * iterator.
     */
    typename LSet<T>::value_type& operator*() const;

    /*!
     * \brief Return a reference to the Lagrangian data item referred to by the
     * iterator.
     */
    typename LSet<T>::value_type& getDataItem() const;

    /*!
     * \brief Return a const reference to the cell index referred to by the
     * iterator.
     */
    const SAMRAI::hier::Index& getCellIndex() const;

private:
    SAMRAI::hier::Box d_box = SAMRAI::hier::Box(IBTK::dim);
    std::shared_ptr<SAMRAI::pdat::IndexIterator<LSet<T>, SAMRAI::pdat::CellGeometry >> d_index_it;
    LSet<T>* d_node_set = nullptr;
    LSetData<T>* d_data = nullptr;
    typename LSet<T>::iterator d_node_it;
};

} // namespace IBTK

/////////////////////////////// INLINE ///////////////////////////////////////

#include "ibtk/private/LSetDataIterator-inl.h" // IWYU pragma: keep

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBTK_LSetDataIterator
