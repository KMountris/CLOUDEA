/*
 * CLOUDEA - Software for solving PDEs using explicit methods.
 * Copyright (C) 2017  <Konstantinos A. Mountris> <konstantinos.mountris@gmail.com>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef CLOUDEA_APPROXIMANTS_FPM_FLUX_CORRECTOR_TPP_
#define CLOUDEA_APPROXIMANTS_FPM_FLUX_CORRECTOR_TPP_


#include "CLOUDEA/engine/approximants/fpm_flux_corrector.hpp"


namespace CLOUDEA {


////!  PUBLIC  ////


template<short DIM>
FpmFluxCorrector<DIM>::FpmFluxCorrector() : facet_ids_(), facet_measures_(),
    facet_centroids_(), facet_dir_(), facet_normals_()
{}


template<short DIM>
FpmFluxCorrector<DIM>::~FpmFluxCorrector()
{}


template<short DIM>
void FpmFluxCorrector<DIM>::ComputeCorrectorData(const IMP::Voronoi<DIM> &voro)
{
    switch (DIM) {
        case 2 :
            this->ComputeCorrectorData2D(voro); break;
        case 3 :
            this->ComputeCorrectorData3D(voro); break;
        default:
            auto err_msg = "Could not compute FPM flux corrector data. Domain dimensions are not supported.";
            throw std::invalid_argument(Logger::Error(err_msg));
        break;
    }
}


template<short DIM>
void FpmFluxCorrector<DIM>::ComputeCorrectorData2D(const IMP::Voronoi<DIM> &voro)
{
    // Initialize storage containers.
    auto fnum = voro.FacetsNum();
    this->facet_ids_.clear();        this->facet_ids_.reserve(fnum);
    this->facet_measures_.clear();   this->facet_measures_.reserve(fnum);
    this->facet_centroids_.clear();  this->facet_centroids_.reserve(fnum);
    this->facet_dir_.clear();        this->facet_dir_.reserve(fnum);
    this->facet_normals_.clear();    this->facet_normals_.reserve(fnum);

    // Iterate over voronoi tesselation facets.
    auto fid = int{0}, n1 = int{0}, n2 = int{0};
    auto facet_centroid = IMP::Vec<DIM,double>{};
    auto facet_normal = IMP::Vec<DIM,double>{};
    auto facet_measure = double{};
    for (const auto &facet : voro.Facets()) {
        // Store the facet index.
        this->facet_ids_.emplace_back(fid);

        // Coordinates of the facet points.
        auto v1 = voro.Points(facet.C(0));
        auto v2 = voro.Points(facet.C(1));

        // The facet centroid
        facet_centroid = 0.5 * (v1+v2);
        this->facet_centroids_.emplace_back(facet_centroid.CopyToEigen());

        // The facet normal.
        // Extract the coordinates of the points of the polyfacet.
        facet_normal = IMP::Vec<DIM,double>({v1[1]-v2[1], v2[0]-v1[0]});

        // Make sure that normal unit vector is outward.
        n1 = facet.ParentCellId();
        n2 = facet.NeighCellId();
        if (n2 != -1) {
            if (facet_normal.CwiseMul(voro.Nodes(n2)-voro.Nodes(n1)).Sum() < 0.)
                facet_normal = -facet_normal;
        }
        this->facet_dir_.emplace_back( facet_normal.CopyToEigen() );

        // Make normal unit vector.
        facet_normal /= facet_normal.Norm();
        this->facet_normals_.emplace_back( facet_normal.CopyToEigen() );

        // Compute the measure of the facet. For 2D domains is the length of the facet.
        facet_measure = voro.Facets(fid).Length(voro.Points());
        this->facet_measures_.emplace_back(facet_measure);

        fid++;
    }

    this->facet_ids_.shrink_to_fit();
    this->facet_measures_.shrink_to_fit();
    this->facet_centroids_.shrink_to_fit();
    this->facet_dir_.shrink_to_fit();
    this->facet_normals_.shrink_to_fit();
}


template<short DIM>
void FpmFluxCorrector<DIM>::ComputeCorrectorData3D(const IMP::Voronoi<DIM> &voro)
{
    // Initialize storage containers.
    auto fnum = voro.FacetsNum();
    this->facet_ids_.clear();        this->facet_ids_.reserve(fnum);
    this->facet_measures_.clear();   this->facet_measures_.reserve(fnum);
    this->facet_centroids_.clear();  this->facet_centroids_.reserve(fnum);
    this->facet_dir_.clear();    this->facet_dir_.reserve(fnum);
    this->facet_normals_.clear();    this->facet_normals_.reserve(fnum);

    // Iterate over voronoi tesselation facets.
    auto fid = int{0}, n1 = int{0}, n2 = int{0};
    auto facet_centroid = IMP::Vec<DIM,double>{};
    auto facet_normal = IMP::Vec<DIM,double>{};
    auto facet_measure = IMP::Vec<DIM, double>{};
    auto oa = IMP::Vec<DIM, double>{};
    auto ob = IMP::Vec<DIM, double>{};
    auto tri_normal = IMP::Vec<DIM, double>{};
    for (const auto &facet : voro.Facets()) {
        // Store the facet index.
        this->facet_ids_.emplace_back(fid);

        // Get the number of the facet points.
        auto facet_pnum = static_cast<int>(facet.Connectivity().size());

        // Compute centroid.
        facet_centroid.SetZero();
        for (auto i = 0; i != facet_pnum; ++i) {
            facet_centroid += voro.Points(facet.C(i));
        }
        facet_centroid /= static_cast<double>(facet_pnum);
        this->facet_centroids_.emplace_back(facet_centroid.CopyToEigen());

        // Compute the measure of the facet. For 3D domains is the area of the facet.
        facet_measure.SetZero();
        facet_normal.SetZero();
        for (auto i = 0; i != facet_pnum; ++i) {
            auto i_next = (i + 1) % facet_pnum;

            // Collect the area contributions of the triangles making up the surface.
            oa = voro.Points(facet.C(i)) - facet_centroid;
            ob = voro.Points(facet.C(i_next)) - facet_centroid;
            tri_normal.Set(
                { IMP::ALGORITHMS::ProdsDiff(oa[1], ob[2], oa[2], ob[1]),
                  IMP::ALGORITHMS::ProdsDiff(oa[2], ob[0], oa[0], ob[2]),
                  IMP::ALGORITHMS::ProdsDiff(oa[0], ob[1], oa[1], ob[0]) }
            );
            facet_measure += 0.5*tri_normal;
            facet_normal += tri_normal;
        }
        this->facet_measures_.emplace_back(facet_measure.Norm());

        // Make sure that normal unit vector is outward.
        n1 = facet.ParentCellId();
        n2 = facet.NeighCellId();
        if (n2 != -1) {
            if (facet_normal.CwiseMul(voro.Nodes(n2)-voro.Nodes(n1)).Sum() < 0.)
                facet_normal = -facet_normal;
        }
        this->facet_dir_.emplace_back( facet_normal.CopyToEigen() );

        // Make normal unit vector.
        facet_normal /= facet_normal.Norm();
        this->facet_normals_.emplace_back( facet_normal.CopyToEigen() );

        fid++;
    }

    this->facet_ids_.shrink_to_fit();
    this->facet_measures_.shrink_to_fit();
    this->facet_centroids_.shrink_to_fit();
    this->facet_dir_.shrink_to_fit();
    this->facet_normals_.shrink_to_fit();

}
// void FpmFluxCorrector<DIM>::ComputeCorrectorData3D(const IMP::Voronoi<DIM> &voro)
// {
//     // Initialize storage containers.
//     auto fnum = voro.FacetsNum();
//     this->facet_ids_.clear();        this->facet_ids_.reserve(fnum);
//     this->facet_measures_.clear();   this->facet_measures_.reserve(fnum);
//     this->facet_centroids_.clear();  this->facet_centroids_.reserve(fnum);
//     this->facet_dir_.clear();    this->facet_dir_.reserve(fnum);
//     this->facet_normals_.clear();    this->facet_normals_.reserve(fnum);

//     // Iterate over voronoi tesselation facets.
//     auto fid = int{0}, n1 = int{0}, n2 = int{0};
//     auto facet_centroid = IMP::Vec<DIM,double>{};
//     auto facet_normal = IMP::Vec<DIM,double>{};
//     auto facet_measure = IMP::Vec<DIM, double>{};
//     auto facet_triangles = std::vector<IMP::Cell<DIM,3>>{};
//     auto tri_edge1 = IMP::Vec<DIM, double>{};
//     auto tri_edge2 = IMP::Vec<DIM, double>{};
//     auto tri_normal = IMP::Vec<DIM, double>{};
//     for (const auto &facet : voro.Facets()) {
//         // Store the facet index.
//         this->facet_ids_.emplace_back(fid);

//         // Get the number of the facet points.
//         auto facet_pnum = static_cast<int>(facet.Connectivity().size());

//         // Initialize centroid of the polygon summing first and last vertices.
//         // rest will be treated during triangle extraction.
//         facet_centroid.SetZero();
//         facet_centroid += voro.Points(facet.C(0));

//         // Initialize facet normal computation with Newell's method
//         // rest will be treated during triangle extraction.
//         auto v1 = facet.C(0);
//         auto v2 = facet.C(1);
//         facet_normal.SetZero();
//         facet_normal[0] += (voro.Points(v1)[1] - voro.Points(v2)[1]) * (voro.Points(v1)[2] + voro.Points(v2)[2]);
//         facet_normal[1] += (voro.Points(v1)[2] - voro.Points(v2)[2]) * (voro.Points(v1)[0] + voro.Points(v2)[0]);
//         facet_normal[2] += (voro.Points(v1)[0] - voro.Points(v2)[0]) * (voro.Points(v1)[1] + voro.Points(v2)[1]);

//         // Extract facet triangles.
//         auto tid = int{0};
//         const auto v0 = facet.C(0);
//         facet_triangles.clear();
//         facet_triangles.resize(facet_pnum-2);
//         for (auto i = 1; i != facet_pnum-1; ++i) {
//             v1 = facet.C(i);
//             v2 = facet.C(i+1);
//             facet_triangles[tid++].SetConnectivity({v0, v1, v2});

//             // Add to the centroid computation.
//             facet_centroid += voro.Points(v1);

//             // Add to the normal computation.
//             facet_normal[0] += (voro.Points(v1)[1] - voro.Points(v2)[1]) * (voro.Points(v1)[2] + voro.Points(v2)[2]);
//             facet_normal[1] += (voro.Points(v1)[2] - voro.Points(v2)[2]) * (voro.Points(v1)[0] + voro.Points(v2)[0]);
//             facet_normal[2] += (voro.Points(v1)[0] - voro.Points(v2)[0]) * (voro.Points(v1)[1] + voro.Points(v2)[1]);
//         }

//         // Complete centroid computation with last point.
//         facet_centroid += voro.Points(v2);
//         facet_centroid /= static_cast<double>(facet_pnum);
//         this->facet_centroids_.emplace_back(facet_centroid.CopyToEigen());

//         // Complete normal computation with pair of last and first points.
//         facet_normal[0] += (voro.Points(v2)[1] - voro.Points(v0)[1]) * (voro.Points(v2)[2] + voro.Points(v0)[2]);
//         facet_normal[1] += (voro.Points(v2)[2] - voro.Points(v0)[2]) * (voro.Points(v2)[0] + voro.Points(v0)[0]);
//         facet_normal[2] += (voro.Points(v2)[0] - voro.Points(v0)[0]) * (voro.Points(v2)[1] + voro.Points(v0)[1]);

//         // Make sure that normal unit vector is outward.
//         n1 = facet.ParentCellId();
//         n2 = facet.NeighCellId();
//         if (n2 != -1) {
//             if (facet_normal.CwiseMul(voro.Nodes(n2)-voro.Nodes(n1)).Sum() < 0.)
//                 facet_normal = -facet_normal;
//         }
//         this->facet_dir_.emplace_back( facet_normal.CopyToEigen() );

//         // Make normal unit vector.
//         facet_normal /= facet_normal.Norm();
//         this->facet_normals_.emplace_back( facet_normal.CopyToEigen() );

//         // Compute the measure of the facet. For 3D domains is the area of the facet.
//         facet_measure.SetZero();
//         for (const auto &tri : facet_triangles) {
//             // Collect the area contributions of the triangles making up the surface.
//             tri_edge1 = voro.Points(tri.N(1)) - voro.Points(tri.N(0));
//             tri_edge2 = voro.Points(tri.N(2)) - voro.Points(tri.N(0));
//             tri_normal.Set({tri_edge1[1]*tri_edge2[2] - tri_edge1[2]*tri_edge2[1],
//                             tri_edge1[2]*tri_edge2[0] - tri_edge1[0]*tri_edge2[2],
//                             tri_edge1[0]*tri_edge2[1] - tri_edge1[1]*tri_edge2[0]});
//             facet_measure += 0.5*tri_normal;
//         }
//         this->facet_measures_.emplace_back(facet_measure.Norm());

//         fid++;
//     }

//     this->facet_ids_.shrink_to_fit();
//     this->facet_measures_.shrink_to_fit();
//     this->facet_centroids_.shrink_to_fit();
//     this->facet_dir_.shrink_to_fit();
//     this->facet_normals_.shrink_to_fit();

// }

} // End of namespace CLOUDEA

#endif //CLOUDEA_APPROXIMANTS_FPM_FLUX_CORRECTOR_TPP_