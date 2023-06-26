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

#ifndef CLOUDEA_APPROXIMANTS_FPM_TPP_
#define CLOUDEA_APPROXIMANTS_FPM_TPP_

#include "CLOUDEA/engine/approximants/fpm.hpp"

namespace CLOUDEA {


////!  PUBLIC  ////


template <short DIM>
Fpm<DIM>::Fpm() : support_(), flux_corrector_(), phi_grad_(), phi_(), penalty_(0.)
{}


template <short DIM>
Fpm<DIM>::~Fpm()
{}


template <short DIM>
void Fpm<DIM>::SetFluxCorrector(const IMP::Voronoi<DIM> &voro)
{
    this->flux_corrector_.ComputeCorrectorData(voro);
}


template <short DIM>
void Fpm<DIM>::Compute(const IMP::Voronoi<DIM> &voro_tess)
{
    this->phi_.clear();       this->phi_.reserve(voro_tess.NodesNum());
    this->phi_grad_.clear();  this->phi_grad_.reserve(voro_tess.NodesNum());

    // Iterate over field nodes.
    auto nid = int{0}, neighs_num = int{0}, cell_cnt = int{0};
    auto cell_center = IMP::Vec<DIM, double>{};
    auto neigh_nodes = std::vector<IMP::Vec<DIM, double>>{};
    for (const auto &node : voro_tess.Nodes()) {
        // Number of influence nodes to the field node.
        // Skip the first influence node as it is the current field node.
        neighs_num = this->support_.InfluenceNodeIds(nid).size() - 1;

        // Get the coordinates of the actual influence nodes.
        neigh_nodes.resize(neighs_num);
        for (auto i = std::size_t{0}; i != neighs_num; ++i) {
            // Get the coordinates of the supporting field node.
            auto neigh_id = this->support_.InfluenceNodeIds(nid)[i+1];
            neigh_nodes[i] = voro_tess.Nodes(neigh_id);
        }

        // Compute the A matrix.
        arma::mat A(neighs_num, DIM);
        for (std::size_t i = 0; i != neighs_num; ++i) {
            for (short j = 0; j != DIM; ++j)
                A(i,j) = neigh_nodes[i][j] - node[j];
        }

        // Compute the C matrix.
        arma::mat C = arma::solve(A.t()*A, A.t());

        // Compute the FPM gradient matrix.
        arma::mat aux = arma::join_cols(-1.*arma::ones(1,neighs_num), arma::eye(neighs_num, neighs_num));

        arma::mat B = aux*C.t();
        this->phi_grad_.emplace_back( this->CastToEigen(B) );

        // Compute cell center to be used as integration point.
        cell_center.SetZero();
        cell_cnt = 0;
        for (const auto &fid : voro_tess.Cells(nid).Connectivity()) {
            for (const auto &pid : voro_tess.Facets(fid).Connectivity()) {
                cell_center += voro_tess.Points(pid);
                cell_cnt++;
            }
        }
        cell_center /= static_cast<double>(cell_cnt);

        // Compute the FPM shape function.
        arma::vec xd(DIM);
        for (short i = 0; i != DIM; ++i)  xd(i) = cell_center[i] - node[i];

        arma::vec fpm_phi = B*xd;
        fpm_phi(0) += 1.;
        this->phi_.emplace_back( this->CastToEigen(fpm_phi) );

        // Increase node index.
        nid++;
    }
}


} // End of namespace CLOUDEA

#endif //CLOUDEA_APPROXIMANTS_FPM_TPP_
