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

#ifndef included_IBTK_muParserCartGridFunction
#define included_IBTK_muParserCartGridFunction

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "ibtk/CartGridFunction.h"
#include "ibtk/ibtk_utilities.h"

#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/hier/PatchLevel.h"


#include "muParser.h"

#include <map>
#include <string>
#include <vector>

namespace SAMRAI
{
namespace hier
{

class Patch;

class Variable;
} // namespace hier
namespace tbox
{
class Database;
} // namespace tbox
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class muParserCartGridFunction is an implementation of the strategy
 * class CartGridFunction that allows for the run-time specification of
 * (possibly spatially- and temporally-varying) functions which are used to set
 * double precision values on standard SAMRAI SAMRAI::hier::PatchData objects.
 */
class muParserCartGridFunction : public CartGridFunction
{
public:
    /*!
     * \brief Constructor.
     */
    muParserCartGridFunction(std::string object_name,
                             std::shared_ptr<SAMRAI::tbox::Database> input_db,
                             std::shared_ptr<SAMRAI::geom::CartesianGridGeometry > grid_geom);

    /*!
     * \brief Empty destructor.
     */
    ~muParserCartGridFunction() = default;

    /*!
     * \name Methods to set patch interior data.
     */
    //\{

    /*!
     * \brief Indicates whether the concrete CartGridFunction object is
     * time-dependent.
     */
    bool isTimeDependent() const override;

    /*!
     * \brief Virtual function to evaluate the function on the patch interior.
     */
    void setDataOnPatch(int data_idx,
                        std::shared_ptr<SAMRAI::hier::Variable > var,
                        std::shared_ptr<SAMRAI::hier::Patch > patch,
                        double data_time,
                        bool initial_time = false,
                        std::shared_ptr<SAMRAI::hier::PatchLevel > level =
                            std::shared_ptr<SAMRAI::hier::PatchLevel >(NULL)) override;

    //\}

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be
     * used.
     */
    muParserCartGridFunction() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    muParserCartGridFunction(const muParserCartGridFunction& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    muParserCartGridFunction& operator=(const muParserCartGridFunction& that) = delete;

    /*!
     * The Cartesian grid geometry object provides the extents of the
     * computational domain.
     */
    std::shared_ptr<SAMRAI::geom::CartesianGridGeometry > d_grid_geom;

    /*!
     * User-provided constants specified in the input file.
     */
    std::map<std::string, double> d_constants;

    /*!
     * The strings providing the data-setting functions which are evaluated by the
     * mu::Parser objects.
     */
    std::vector<std::string> d_function_strings;

    /*!
     * The mu::Parser objects which evaluate the data-setting functions.
     */
    std::vector<mu::Parser> d_parsers;

    /*!
     * Time and position variables.
     */
    double d_parser_time = 0.0;
    Point d_parser_posn;
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBTK_muParserCartGridFunction
