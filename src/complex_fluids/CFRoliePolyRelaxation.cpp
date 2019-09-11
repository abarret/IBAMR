#include "ibamr/CFRoliePolyRelaxation.h"

// Namespace
namespace IBAMR
{
// Constructor
CFRoliePolyRelaxation::CFRoliePolyRelaxation(const std::string& object_name, Pointer<Database> input_db)
    : CFRelaxationOperator(object_name)
{
    d_lambda_d = input_db->getDouble("lambda_d");
    d_lambda_R = input_db->getDouble("lambda_R");
    d_beta = input_db->getDouble("beta");
    d_delta = input_db->getDouble("delta");
    d_sqr_root = input_db->getBoolWithDefault("square_root_evolve", false);
    d_log_conform = input_db->getBoolWithDefault("log_conform_evolve", false);
    return;
} // End Constructor

void
CFRoliePolyRelaxation::setDataOnPatch(const int data_idx,
                                      Pointer<Variable<NDIM>> /*var*/,
                                      Pointer<Patch<NDIM>> patch,
                                      const double /*data_time*/,
                                      const bool initial_time,
                                      Pointer<PatchLevel<NDIM>> /*patch_level*/)
{
    const Box<NDIM>& patch_box = patch->getBox();
    const Pointer<CartesianPatchGeometry<NDIM>> p_geom = patch->getPatchGeometry();
    Pointer<CellData<NDIM, double>> ret_data = patch->getPatchData(data_idx);
    Pointer<CellData<NDIM, double>> in_data = patch->getPatchData(d_W_cc_idx);
    ret_data->fillAll(0.0);
    double tr = 0.0;
    if (initial_time) return;
    for (CellIterator<NDIM> i(patch_box); i; i++)
    {
        const CellIndex<NDIM>& idx = i();
#if (NDIM == 2)
        double Q_xx, Q_yy, Q_xy;
        if (d_sqr_root)
        {
            Q_xx = (*in_data)(idx, 0) * (*in_data)(idx, 0) + (*in_data)(idx, 2) * (*in_data)(idx, 2);
            Q_yy = (*in_data)(idx, 2) * (*in_data)(idx, 2) + (*in_data)(idx, 1) * (*in_data)(idx, 1);
            Q_xy = (*in_data)(idx, 0) * (*in_data)(idx, 2) + (*in_data)(idx, 1) * (*in_data)(idx, 2);
        }
        else if (d_log_conform)
        {
            Eigen::Matrix2d mat;
            mat(0, 0) = (*in_data)(idx, 0);
            mat(1, 1) = (*in_data)(idx, 1);
            mat(0, 1) = mat(1, 0) = (*in_data)(idx, 2);
            mat = mat.exp();
            Q_xx = mat(0, 0);
            Q_yy = mat(1, 1);
            Q_xy = mat(0, 1);
        }
        else
        {
            Q_xx = (*in_data)(idx, 0);
            Q_yy = (*in_data)(idx, 1);
            Q_xy = (*in_data)(idx, 2);
        }
        tr = Q_xx + Q_yy;
        (*ret_data)(idx, 0) =
            -1.0 / d_lambda_d * (Q_xx - 1.0) -
            2.0 * (1.0 - sqrt(2.0 / tr)) / d_lambda_R * (Q_xx + d_beta * pow(tr / 2.0, d_delta) * (Q_xx - 1.0));
        (*ret_data)(idx, 1) =
            -1.0 / d_lambda_d * (Q_yy - 1.0) -
            2.0 * (1.0 - sqrt(2.0 / tr)) / d_lambda_R * (Q_yy + d_beta * pow(tr / 2.0, d_delta) * (Q_yy - 1.0));
        (*ret_data)(idx, 2) = -1.0 / d_lambda_d * (Q_xy)-2.0 * (1.0 - sqrt(2.0 / tr)) / d_lambda_R *
                              (Q_xy + d_beta * pow(tr / 2.0, d_delta) * (Q_xy));
#endif
#if (NDIM == 3)
        double Q_xx, Q_yy, Q_zz, Q_yz, Q_xz, Q_xy;
        if (d_sqr_root)
        {
            Q_xx = (*in_data)(idx, 0) * (*in_data)(idx, 0) + (*in_data)(idx, 5) * (*in_data)(idx, 5) +
                   (*in_data)(idx, 4) * (*in_data)(idx, 4);
            Q_yy = (*in_data)(idx, 5) * (*in_data)(idx, 5) + (*in_data)(idx, 1) * (*in_data)(idx, 1) +
                   (*in_data)(idx, 3) * (*in_data)(idx, 3);
            Q_zz = (*in_data)(idx, 4) * (*in_data)(idx, 4) + (*in_data)(idx, 3) * (*in_data)(idx, 3) +
                   (*in_data)(idx, 2) * (*in_data)(idx, 2);
            Q_yz = (*in_data)(idx, 5) * (*in_data)(idx, 4) + (*in_data)(idx, 1) * (*in_data)(idx, 3) +
                   (*in_data)(idx, 3) * (*in_data)(idx, 2);
            Q_xz = (*in_data)(idx, 0) * (*in_data)(idx, 4) + (*in_data)(idx, 5) * (*in_data)(idx, 3) +
                   (*in_data)(idx, 4) * (*in_data)(idx, 2);
            Q_xy = (*in_data)(idx, 0) * (*in_data)(idx, 5) + (*in_data)(idx, 5) * (*in_data)(idx, 1) +
                   (*in_data)(idx, 4) * (*in_data)(idx, 3);
        }
        else if (d_log_conform)
        {
            Eigen::Matrix3d mat;
            mat(0, 0) = (*in_data)(idx, 0);
            mat(1, 1) = (*in_data)(idx, 1);
            mat(2, 2) = (*in_data)(idx, 2);
            mat(1, 2) = mat(2, 1) = (*in_data)(idx, 3);
            mat(0, 2) = mat(2, 0) = (*in_data)(idx, 4);
            mat(0, 1) = mat(1, 0) = (*in_data)(idx, 5);
            mat = mat.exp();
            Q_xx = mat(0, 0);
            Q_yy = mat(1, 1);
            Q_zz = mat(2, 2);
            Q_xy = mat(0, 1);
            Q_xz = mat(0, 2);
            Q_yz = mat(1, 2);
        }

        else
        {
            Q_xx = (*in_data)(idx, 0);
            Q_yy = (*in_data)(idx, 1);
            Q_zz = (*in_data)(idx, 2);
            Q_yz = (*in_data)(idx, 3);
            Q_xz = (*in_data)(idx, 4);
            Q_xy = (*in_data)(idx, 5);
        }
        tr = Q_xx + Q_yy + Q_zz;
        (*ret_data)(idx, 0) =
            -1.0 / d_lambda_d * (Q_xx - 1.0) -
            2.0 * (1.0 - sqrt(3.0 / tr)) / d_lambda_R * (Q_xx + d_beta * pow(tr / 3.0, d_delta) * (Q_xx - 1.0));
        (*ret_data)(idx, 1) =
            -1.0 / d_lambda_d * (Q_yy - 1.0) -
            2.0 * (1.0 - sqrt(3.0 / tr)) / d_lambda_R * (Q_yy + d_beta * pow(tr / 3.0, d_delta) * (Q_yy - 1.0));
        (*ret_data)(idx, 2) =
            -1.0 / d_lambda_d * (Q_zz - 1.0) -
            2.0 * (1.0 - sqrt(3.0 / tr)) / d_lambda_R * (Q_zz + d_beta * pow(tr / 3.0, d_delta) * (Q_zz - 1.0));
        (*ret_data)(idx, 3) = -1.0 / d_lambda_d * (Q_yz)-2.0 * (1.0 - sqrt(3.0 / tr)) / d_lambda_R *
                              (Q_yz + d_beta * pow(tr / 3.0, d_delta) * (Q_yz));
        (*ret_data)(idx, 4) = -1.0 / d_lambda_d * (Q_xz)-2.0 * (1.0 - sqrt(3.0 / tr)) / d_lambda_R *
                              (Q_xz + d_beta * pow(tr / 3.0, d_delta) * (Q_xz));
        (*ret_data)(idx, 5) = -1.0 / d_lambda_d * (Q_xy)-2.0 * (1.0 - sqrt(3.0 / tr)) / d_lambda_R *
                              (Q_xy + d_beta * pow(tr / 3.0, d_delta) * (Q_xy));
#endif
    }
} // end setDataOnPatch

} // namespace IBAMR
