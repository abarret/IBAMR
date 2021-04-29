#include "ibamr/cut_cells/LSFromMesh.h"
#include "ibamr/cut_cells/ls_functions.h"

#include "ibtk/DebuggingUtilities.h"
#include "ibtk/IBTK_MPI.h"

#include <CGAL/Alpha_shape_2.h>
#include <CGAL/Alpha_shape_face_base_2.h>
#include <CGAL/Alpha_shape_vertex_base_2.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Quotient.h>
#include <CGAL/algorithm.h>
#include <CGAL/assertions.h>
#include <CGAL/draw_polygon_2.h>
#include <CGAL/draw_triangulation_2.h>

#include <algorithm>

// FORTRAN ROUTINES
#if (NDIM == 2)
#define SIGN_SWEEP_FC IBAMR_FC_FUNC(signsweep2dn, SIGNSWEEP2D)
#endif

extern "C"
{
    void SIGN_SWEEP_FC(double* U,
                       const int& U_gcw,
                       const int& ilower0,
                       const int& iupper0,
                       const int& ilower1,
                       const int& iupper1,
                       const double& large_dist,
                       int& n_updates);
}

namespace
{
static Timer* t_updateVolumeAreaSideLS = nullptr;
static Timer* t_findIntersection = nullptr;
} // namespace

namespace LS
{
const double LSFromMesh::s_eps = 0.5;

LSFromMesh::LSFromMesh(std::string object_name,
                       Pointer<PatchHierarchy<NDIM> > hierarchy,
                       MeshBase* mesh,
                       FEDataManager* fe_data_manager,
                       bool use_inside /* = true*/)
    : LSFindCellVolume(std::move(object_name), hierarchy),
      d_mesh(mesh),
      d_fe_data_manager(fe_data_manager),
      d_use_inside(use_inside),
      d_sgn_var(new CellVariable<NDIM, double>(d_object_name + "SGN"))
{
    IBAMR_DO_ONCE(t_updateVolumeAreaSideLS =
                      TimerManager::getManager()->getTimer("LS::LSFromMesH::updateVolumeAreaSideLS()");
                  t_findIntersection = TimerManager::getManager()->getTimer("LS::LSFromMesh::findIntersection()"););
    d_cut_cell_mesh_mapping = libmesh_make_unique<CutCellMeshMapping>(
        d_object_name + "::CutCellMapping", nullptr, static_cast<Mesh*>(d_mesh), d_fe_data_manager);
    return;
} // Constructor

void
LSFromMesh::updateVolumeAreaSideLS(int vol_idx,
                                   Pointer<CellVariable<NDIM, double> > /*vol_var*/,
                                   int area_idx,
                                   Pointer<CellVariable<NDIM, double> > /*area_var*/,
                                   int side_idx,
                                   Pointer<SideVariable<NDIM, double> > /*side_var*/,
                                   int phi_idx,
                                   Pointer<NodeVariable<NDIM, double> > phi_var,
                                   double /*data_time*/,
                                   bool extended_box)
{
    LS_TIMER_START(t_updateVolumeAreaSideLS);
    TBOX_ASSERT(phi_idx != IBTK::invalid_index);
    TBOX_ASSERT(phi_var);
    HierarchyNodeDataOpsReal<NDIM, double> hier_nc_data_ops(d_hierarchy, 0, d_hierarchy->getFinestLevelNumber());
    hier_nc_data_ops.setToScalar(phi_idx, s_eps, false);
    HierarchyCellDataOpsReal<NDIM, double> hier_cc_data_ops(d_hierarchy, 0, d_hierarchy->getFinestLevelNumber());
    hier_cc_data_ops.setToScalar(vol_idx, 0.0, false);

    d_cut_cell_mesh_mapping->initializeObjectState(d_hierarchy);
    d_cut_cell_mesh_mapping->generateCutCellMappings();
    const std::vector<std::map<IndexList, std::vector<CutCellElems> > >& idx_cut_cell_map_vec =
        d_cut_cell_mesh_mapping->getIdxCutCellElemsMap(d_hierarchy->getFinestLevelNumber());

    Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(d_hierarchy->getFinestLevelNumber());
    unsigned int local_patch_num = 0;
    for (PatchLevel<NDIM>::Iterator p(level); p; p++, ++local_patch_num)
    {
        Pointer<Patch<NDIM> > patch = level->getPatch(p());
        const std::map<IndexList, std::vector<CutCellElems> >& idx_cut_cell_map = idx_cut_cell_map_vec[local_patch_num];
        if (idx_cut_cell_map.size() == 0) continue;

        Pointer<CellData<NDIM, double> > vol_data = patch->getPatchData(vol_idx);
        Pointer<CellData<NDIM, double> > area_data = patch->getPatchData(area_idx);
        Pointer<SideData<NDIM, double> > side_data = patch->getPatchData(side_idx);
        Pointer<NodeData<NDIM, double> > phi_data = patch->getPatchData(phi_idx);

        Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
        const double* const x_low = pgeom->getXLower();
        const double* const dx = pgeom->getDx();
        const hier::Index<NDIM>& idx_low = patch->getBox().lower();

        const IntVector<NDIM>& lower = patch->getBox().lower();
        const IntVector<NDIM>& upper = patch->getBox().upper();

        // Loop through all cut cells in map
        for (const auto& idx_elem_vec_pair : idx_cut_cell_map)
        {
            const CellIndex<NDIM>& idx = idx_elem_vec_pair.first.d_idx;
            plog << "On index: " << idx << "\n";
            const std::vector<CutCellElems>& cut_cell_elem_vec = idx_elem_vec_pair.second;

            // Determine area contributions
            if (area_data)
            {
                double area = 0.0;
                for (const auto& cut_cell_elem : cut_cell_elem_vec)
                {
                    area += (cut_cell_elem.d_elem->point(1) - cut_cell_elem.d_elem->point(0)).norm();
                }
                (*area_data)(idx) = area;
            }

            for (const auto& cut_cell_elem : cut_cell_elem_vec)
            {
                const std::unique_ptr<Elem>& elem = cut_cell_elem.d_elem;
                const Elem* const parent_elem = cut_cell_elem.d_parent_elem;
                const std::vector<std::pair<libMesh::Point, int> >& pt_side_vec = cut_cell_elem.d_intersection_side_vec;

                plog << "Cut cell has points: " << elem->point(0) << " and " << elem->point(1) << "\n";
                plog << "Parent elem has points: " << parent_elem->point(0) << " and " << parent_elem->point(1) << "\n";
                plog << "We have the following intersections:\n";
                for (const auto& pt_side : pt_side_vec)
                    plog << "Pt: " << pt_side.first << " and side: " << pt_side.second << "\n";
            }
            plog << "\n\n";

            // Let's find distances from background cell nodes to structure
            // Determine normal for elements
            // Warning, we need the normal to be consistent between parent and child elements.
            std::vector<IBTK::Vector3d> elem_normals;
            int i = 0;
            for (const auto& cut_cell_elem : cut_cell_elem_vec)
            {
                // Note we use the parent element to calculate normals to preserve directions
                Vector3d v, w;
                const Elem* const elem = cut_cell_elem.d_parent_elem;
                v << elem->point(0)(0), elem->point(0)(1), elem->point(0)(2);
                w << elem->point(1)(0), elem->point(1)(1), elem->point(1)(2);
                const unsigned int domain_id = cut_cell_elem.d_parent_elem->subdomain_id();
                Vector3d e3 = Vector3d::UnitZ();
                if (!d_use_inside) e3 *= -1.0;
                if (d_norm_reverse_domain_id.find(domain_id) != d_norm_reverse_domain_id.end() ||
                    d_norm_reverse_elem_id.find(cut_cell_elem.d_parent_elem->id()) != d_norm_reverse_elem_id.end())
                    e3 *= -1.0;
                Vector3d n = (w - v).cross(e3);
                plog << "Parent normal is: " << n.transpose() << "\n";
                elem_normals.push_back(n);
            }
            // Determine distances to nodes
            for (int x = 0; x < 2; ++x)
            {
                for (int y = 0; y < 2; ++y)
                {
                    Vector3d P = Vector3d::Zero();
                    for (int d = 0; d < NDIM; ++d) P(d) = static_cast<double>(idx(d) - idx_low(d)) + (d == 0 ? x : y);
                    // Project P onto element
                    Vector3d avg_proj, avg_unit_normal;
                    avg_proj.setZero();
                    avg_unit_normal.setZero();
                    double min_dist = std::numeric_limits<double>::max();
                    int num_min = 0;
                    // Loop through all elements and calculate the smallest distance
                    for (unsigned int i = 0; i < elem_normals.size(); ++i)
                    {
                        const std::unique_ptr<Elem>& elem = cut_cell_elem_vec[i].d_elem;
                        const Vector3d& n = elem_normals[i];
                        Vector3d v, w;
                        v << (elem->point(0)(0) - x_low[0]) / dx[0], (elem->point(0)(1) - x_low[1]) / dx[1], 0.0;
                        w << (elem->point(1)(0) - x_low[0]) / dx[0], (elem->point(1)(1) - x_low[1]) / dx[1], 0.0;
                        const double t = std::max(0.0, std::min(1.0, (P - v).dot(w - v) / (v - w).squaredNorm()));
                        const Vector3d proj = v + t * (w - v);
                        VectorNd x_proj;
                        for (int d = 0; d < NDIM; ++d) x_proj[d] = x_low[d] + dx[d] * (proj(d) - idx_low(d));
                        const double dist = (proj - P).norm();
                        if (dist < min_dist)
                        {
                            min_dist = dist;
                            avg_proj = proj;
                            avg_unit_normal = n;
                            num_min = 1;
                        }
                        else if (MathUtilities<double>::equalEps(dist, min_dist))
                        {
                            avg_proj += proj;
                            avg_unit_normal += n;
                            ++num_min;
                        }
                    }
                    avg_proj /= static_cast<double>(num_min);
                    avg_unit_normal /= static_cast<double>(num_min);
                    avg_unit_normal.normalize();

                    const double dist = (P - avg_proj).norm();
                    Vector3d phys_vec;
                    for (unsigned int d = 0; d < NDIM; ++d) phys_vec(d) = dx[d] * (P - avg_proj)[d];
                    double dist_phys = phys_vec.norm();
                    double sgn = (avg_unit_normal.dot(P - avg_proj) <= 0.0 ? -1.0 : 1.0);
                    NodeIndex<NDIM> n_idx(idx, IntVector<NDIM>(x, y));
                    (*phi_data)(n_idx) =
                        dist_phys < std::abs((*phi_data)(n_idx)) ? (dist_phys * sgn) : (*phi_data)(n_idx);
                }
            }

            // Find cell volumes. We use a triangulation algorithm from CGAL
            // First find all vertices and put them in one vector
            using K = CGAL::Exact_predicates_inexact_constructions_kernel;
            using Triangulation = CGAL::Triangulation_2<K>;
            bool draw = false;
            if (idx(0) == 128 && idx(1) == 67) draw = true;
            Triangulation t;
            t.clear();
            std::vector<K::Point_2> points;
            // First vertices of cell nodes
            VectorNd x_pt;
            for (int d = 0; d < NDIM; ++d) x_pt[d] = x_low[d] + dx[d] * static_cast<double>(idx(d) - idx_low(d));
            NodeIndex<NDIM> idx_ll(idx, NodeIndex<NDIM>::LowerLeft);
            NodeIndex<NDIM> idx_lr(idx, NodeIndex<NDIM>::LowerRight);
            NodeIndex<NDIM> idx_ur(idx, NodeIndex<NDIM>::UpperRight);
            NodeIndex<NDIM> idx_ul(idx, NodeIndex<NDIM>::UpperLeft);
            if ((*phi_data)(idx_ll) < 0.0) points.push_back(K::Point_2(x_pt[0], x_pt[1]));
            if ((*phi_data)(idx_lr) < 0.0) points.push_back(K::Point_2(x_pt[0] + dx[0], x_pt[1]));
            if ((*phi_data)(idx_ur) < 0.0) points.push_back(K::Point_2(x_pt[0] + dx[0], x_pt[1] + dx[1]));
            if ((*phi_data)(idx_ul) < 0.0) points.push_back(K::Point_2(x_pt[0], x_pt[1] + dx[1]));
            // Now loop through all intersections
            // Only count interior points once
            std::vector<libMesh::Point> int_pts;
            std::map<int, std::vector<libMesh::Point> > side_pt_map;
            for (const auto& cut_cell : cut_cell_elem_vec)
            {
                const std::vector<std::pair<libMesh::Point, int> >& pt_side_vec = cut_cell.d_intersection_side_vec;
                for (const auto& pt_side : pt_side_vec)
                {
                    const libMesh::Point& pt = pt_side.first;
                    if (pt_side.second == -1 && std::find(int_pts.begin(), int_pts.end(), pt) != int_pts.end())
                        continue; // Already found point.
                    else
                        int_pts.push_back(pt); // Haven't found point yet.
                    side_pt_map[pt_side.second].push_back(pt);
                    points.push_back(K::Point_2(pt(0), pt(1)));
                }
            }
            // If we have multiple intersections on a side, we need to check whether that side has INTERIOR cell nodes.
            // If there are INTERIOR cell nodes, our structure is concave, and we need to use a different algorithm
            bool concave_structure = false;
            for (const auto& side_pt_vec_pair : side_pt_map)
            {
                int side = side_pt_vec_pair.first;
                if (side == -1) continue;
                if (side_pt_vec_pair.second.size() > 1)
                {
                    NodeIndex<NDIM> idx_l, idx_u;
                    switch (side)
                    {
                    case 0:
                        idx_l = NodeIndex<NDIM>(idx, NodeIndex<NDIM>::LowerLeft);
                        idx_u = NodeIndex<NDIM>(idx, NodeIndex<NDIM>::UpperLeft);
                        break;
                    case 1:
                        idx_l = NodeIndex<NDIM>(idx, NodeIndex<NDIM>::LowerRight);
                        idx_u = NodeIndex<NDIM>(idx, NodeIndex<NDIM>::UpperRight);
                        break;
                    case 2:
                        idx_l = NodeIndex<NDIM>(idx, NodeIndex<NDIM>::LowerLeft);
                        idx_u = NodeIndex<NDIM>(idx, NodeIndex<NDIM>::LowerRight);
                        break;
                    case 3:
                        idx_l = NodeIndex<NDIM>(idx, NodeIndex<NDIM>::UpperLeft);
                        idx_u = NodeIndex<NDIM>(idx, NodeIndex<NDIM>::UpperRight);
                        break;
                    default:
                        TBOX_ERROR("Unknown side!\n");
                    }
                    if ((*phi_data)(idx_l) < 0.0 && (*phi_data)(idx_u) < 0.0) concave_structure = true;
                }
            }
            // We now have all the points in the vector. Use the triangulation algorithm.
            t.insert(points.begin(), points.end());
            if (draw) CGAL::draw(t);
            Triangulation::Finite_faces_iterator face_iter = t.finite_faces_begin();
            double vol = 0.0;
            for (; face_iter != t.finite_faces_end(); ++face_iter)
            {
                Triangulation::Triangle tri = t.triangle(face_iter);
                vol += tri.area();
            }
            (*vol_data)(idx) = vol / (dx[0] * dx[1]);
        }
    }
    LS_TIMER_STOP(t_updateVolumeAreaSideLS);
}

bool
LSFromMesh::findIntersection(libMesh::Point& p, Elem* elem, libMesh::Point r, libMesh::VectorValue<double> q)
{
    LS_TIMER_START(t_findIntersection);
    bool found_intersection = false;
    switch (elem->type())
    {
    case libMesh::EDGE2:
    {
        // Use linear interpolation
        // Plane through r in q direction:
        // p = r + t * q
        // Plane through two element points p0, p1
        // p = 0.5*(1+u)*p0 + 0.5*(1-u)*p1
        // Set equal and solve for u and t.
        // Note that since q is aligned with a grid axis, we can solve for u first, then find t later
        // Solve for u via a * u + b = 0
        // with a = 0.5 * (p0 - p1)
        //      b = 0.5 * (p0 + p1) - r
        const libMesh::Point& p0 = elem->point(0);
        const libMesh::Point& p1 = elem->point(1);
        const int search_dir = q(0) == 0.0 ? 1 : 0;
        const int trans_dir = (search_dir + 1) % NDIM;
        double a = 0.5 * (p0(trans_dir) - p1(trans_dir));
        double b = 0.5 * (p0(trans_dir) + p1(trans_dir)) - r(trans_dir);
        const double u = -b / a;
        // Determine if this intersection is on the interior of the element
        // This means that u is between -1 and 1
        if (u >= -1.0 && u <= 1.0)
        {
            // Now determine if intersection occurs on axis
            // This amounts to t being between -0.5 and 0.5
            double p_search = 0.5 * p0(search_dir) * (1.0 + u) + 0.5 * (1.0 - u) * p1(search_dir);
            double t = (p_search - r(search_dir)) / q(search_dir);
            if (t >= -0.5 && t <= 0.5)
            {
                // We've found an intersection on this axis
                p = 0.5 * (1.0 + u) * p0 + 0.5 * (1.0 - u) * p1;
                found_intersection = true;
            }
        }
        break;
    }
    default:
        TBOX_ERROR("Unknown element.\n");
    }
    LS_TIMER_STOP(t_findIntersection);
    return found_intersection;
}
} // namespace LS
