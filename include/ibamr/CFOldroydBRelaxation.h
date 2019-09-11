
#ifndef included_CFOldroydBRelaxation
#define included_CFOldroydBRelaxation
#include "IBAMR_config.h"

#include "ibamr/CFRelaxationOperator.h"
#include "ibamr/ibamr_enums.h"
#include "ibamr/namespaces.h" // IWYU pragma: keep

#include "Box.h"
#include "CartesianPatchGeometry.h"
#include "CellData.h"
#include "CellIndex.h"
#include "CellVariable.h"
#include "Index.h"
#include "IntVector.h"
#include "Patch.h"
#include "PatchGeometry.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "Variable.h"
#include "tbox/Database.h"
#include "tbox/Pointer.h"
#include "tbox/Utilities.h"

#include <unsupported/Eigen/MatrixFunctions>

#include <math.h>

#include <string>

namespace SAMRAI
{
namespace hier
{
template <int DIM>
class Patch;
template <int DIM>
class PatchHierarchy;
template <int DIM>
class Variable;
} // namespace hier
namespace tbox
{
class Database;
} // namespace tbox
} // namespace SAMRAI
/////////////////////////////// INCLUDES /////////////////////////////////////

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class CFOldroydBRelaxation is a concrete CFRelaxationOperator
 * that computes the relaxation function for an Oldroyd-B fluid model.
 */
class CFOldroydBRelaxation : public IBAMR::CFRelaxationOperator
{
public:
    /*!
     * \brief This constructor reads in the parameters for the model from the input database.
     */
    CFOldroydBRelaxation(const std::string& object_name, Pointer<Database> input_db);

    /*!
     * \brief Empty destructor.
     */
    ~CFOldroydBRelaxation() = default;

    /*!
     * \name Methods to set patch data.
     */
    //\{

    /*!
     * \brief Evaluate the function on the patch interior.
     */
    void setDataOnPatch(const int data_idx,
                        SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM>> var,
                        SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM>> patch,
                        const double data_time,
                        const bool initial_time = false,
                        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> patch_level =
                            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>>(NULL)) override;

    //\}

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    CFOldroydBRelaxation();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    CFOldroydBRelaxation(const CFOldroydBRelaxation& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    CFOldroydBRelaxation& operator=(const CFOldroydBRelaxation& that);

    double d_lambda;
    bool d_sqr_root, d_log_conform;
}; // Private

} // Namespace IBAMR
#endif
