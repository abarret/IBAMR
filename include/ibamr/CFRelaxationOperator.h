
#ifndef included_CFRelaxationOperator
#define included_CFRelaxationOperator
#include "IBAMR_config.h"

#include "ibamr/ibamr_enums.h"
#include "ibamr/namespaces.h" // IWYU pragma: keep
#include <ibamr/app_namespaces.h>

#include "ibtk/CartGridFunction.h"
#include "ibtk/ibtk_utilities.h"

#include "Patch.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "tbox/Pointer.h"
#include "tbox/Utilities.h"

#include <string>
#include <vector>

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
 * \brief Class CFRelaxationOperator is an abstract class that provides an interface
 * for specifying a relaxation function for the extra stress for Oldroyd-B type viscoelastic fluid models.
 */
class CFRelaxationOperator : public IBTK::CartGridFunction
{
public:
    /*!
     * \brief This constructor does nothing interesting.
     */
    CFRelaxationOperator(const std::string& object_name);

    /*!
     * \brief Empty destructor.
     */
    virtual ~CFRelaxationOperator() = default;

    /*!
     * \name Methods to set patch data.
     */
    //\{

    /*!
     * \brief Indicates whether the concrete CFRelaxationOperator object is
     * time-dependent. Returns true.
     */
    bool isTimeDependent() const override;

    virtual void setPatchDataIndex(const int data_idx);

    //\}

protected:
    /*!
     * The object name is used for error/warning reporting.
     */
    std::string d_object_name;
    int d_W_cc_idx = IBTK::invalid_index;

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    CFRelaxationOperator();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    CFRelaxationOperator(const CFRelaxationOperator& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    CFRelaxationOperator& operator=(const CFRelaxationOperator& that);
}; // Private

} // Namespace IBAMR
#endif
