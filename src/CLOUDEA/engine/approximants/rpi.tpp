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

#ifndef CLOUDEA_APPROXIMANTS_RPI_TPP_
#define CLOUDEA_APPROXIMANTS_RPI_TPP_

#include "CLOUDEA/engine/approximants/rpi.hpp"

namespace CLOUDEA {


template <short DIM>
Rpi<DIM>::Rpi()
{
    this->type_ = MfreeType::rpi;
}


template <short DIM>
Rpi<DIM>::~Rpi()
{}


template <short DIM>
void Rpi<DIM>::ComputeEnrichedRbf(const IMP::Vec<DIM, double> &eval_point, const std::vector<IMP::Vec<DIM, double>> &neigh_nodes,
                                  const std::vector<int> &neigh_ids, arma::vec &rpi_values, arma::mat &rpi_grad)
{

    // Compute monomial basis for eval_point.
    this->mbasis_->Compute(eval_point);
    short m = this->mbasis_->Rank();

    // Get the number of neighbor field nodes in the support domain.
    std::size_t neighs_num = neigh_nodes.size();

    // Reset values of the radial basis function and its gradient.
    rpi_values = arma::vec(neighs_num + m, arma::fill::zeros);
    rpi_grad   = arma::mat(neighs_num + m, DIM);

    // Iterate over neighbor nodes.
    for (std::size_t i = 0; i != neighs_num; ++i) {

        // Compute the weight function of neighbor node at eval_point.
        this->weight_->Compute(eval_point, neigh_nodes[i], this->support_.Radius(neigh_ids[i]), this->support_.DilateCoeff(neigh_ids[i]));

        // Assign the rbf value.
        rpi_values(i) = this->weight_->Val();

        // Assign the gradient of the rbf value.
        for (std::size_t d = 0; d != DIM; ++d) {
            rpi_grad(i, d) = this->weight_->Grad(d);
        }
    }
    
    // Set the values and gradient for the monomial basis enrichment.
    for (short j = 0; j != m; ++j) {
        rpi_values(neighs_num+j) = this->mbasis_->Basis(0)(j);

        for (short d = 0; d != DIM; ++d) {
            rpi_grad(neighs_num+j, d) = this->mbasis_->Grad(0)(d, j);
        }
    }
}


template <short DIM>
void Rpi<DIM>::Compute(const std::vector<IMP::Vec<DIM, double>> &eval_points,
                       const std::vector<IMP::Vec<DIM, double>> &field_nodes)
{
    std::size_t eval_points_num = eval_points.size();
    this->phi_.clear();       this->phi_.reserve(eval_points_num);
    this->phi_grad_.clear();  this->phi_grad_.reserve(eval_points_num);

    // Get monomial basis rank.
    short m = this->mbasis_->Rank();

    // Iterate over eval_points.
    for (const auto &point : eval_points) {

        // Get the index of the current evaluation point.
        auto eval_point_id = &point - &eval_points[0];

        // Number of neighbor field nodes to the evaluation point.
        auto neighs_num = this->support_.InfluenceNodeIds(eval_point_id).size();

        // Get the coordinates of the neighbor nodes in the support domain.
        std::vector<IMP::Vec<DIM, double>> neigh_nodes(neighs_num, IMP::Vec<DIM, double>());
        for (const auto &neigh_id : this->support_.InfluenceNodeIds(eval_point_id)) {
        
            // Get the coordinates of the supporting field node.
            auto i = &neigh_id - &this->support_.InfluenceNodeIds(eval_point_id)[0];
            neigh_nodes[i] = field_nodes[neigh_id];
        }

        // Initialize radial basis function values and gradient.
        arma::vec rho(neighs_num + m);
        arma::mat drho(neighs_num + m, DIM);

        // Compute the Vandermonde matrix for the neighbor nodes.
        this->mbasis_->Compute(neigh_nodes);

        // Iterate over the neighbor nodes to assemble the RPI matrix G.
        arma::mat G(neighs_num + m, neighs_num + m, arma::fill::zeros);
        for (const auto &neigh_id : this->support_.InfluenceNodeIds(eval_point_id)) {

            // Get position of the iterator.
            auto i = &neigh_id - &this->support_.InfluenceNodeIds(eval_point_id)[0];

            // Compute the values and gradient of the rbf for the current field node in the support domain.
            this->ComputeEnrichedRbf(field_nodes[neigh_id], neigh_nodes, this->support_.InfluenceNodeIds(eval_point_id), rho, drho);

            // Fill the G matrix with the rbf values.
            for (std::size_t j = 0; j != neighs_num; ++j) { G(j, i) = rho(j); }

            // Fill the G matrix with the monomial basis values.
            for (short j = 0; j != m; ++j) {
                G(i, neighs_num+j) = this->mbasis_->Basis(i)(j);
                G(neighs_num+j, i) = this->mbasis_->Basis(i)(j);
            }
        } // End of Iterate over the support nodes to assemble the RPI matrix G.

        // Compute the values and gradient of the rbf for the evaluation point.
        this->ComputeEnrichedRbf(point, neigh_nodes, this->support_.InfluenceNodeIds(eval_point_id), rho, drho);

        // Get the RPI basis function.
        arma::vec rpi_phi = arma::solve(G, rho);
        rpi_phi = rpi_phi.head_rows(neighs_num);
        this->phi_.emplace_back( this->CastToEigen(rpi_phi) );

        // Get the RPI basis function gradient.
        arma::mat rpi_phi_grad(neighs_num + m, DIM);
        for (int d = 0; d != DIM; ++d) {
            rpi_phi_grad.col(d) = arma::solve(G, drho.col(d));
        }
        rpi_phi_grad = rpi_phi_grad.head_rows(neighs_num);
        this->phi_grad_.emplace_back( this->CastToEigen(rpi_phi_grad) );
    }

}


template <short DIM>
void Rpi<DIM>::Compute(const IMP::Vec<DIM, double> &eval_point, int eval_point_id,
                       const std::vector<IMP::Vec<DIM, double>> &field_nodes)
{
    // Number of neighbor field nodes to the evaluation point.
    std::size_t neighs_num = this->support_.InfluenceNodeIds(eval_point_id).size();

    // Reset RPI basis function and gradient to zero for a single evaluation point.
    this->phi_.clear();
    this->phi_.resize(1, Eigen::VectorXd::Zero(neighs_num));
    this->phi_grad_.clear();
    this->phi_grad_.resize(1, Eigen::MatrixXd::Zero(neighs_num, DIM));

    // Get the monomial basis rank.
    short m = this->mbasis_->Rank();
    
    // Get the coordinates of the neighbor nodes in the support domain.
    std::vector<IMP::Vec<DIM, double>> neigh_nodes(neighs_num, IMP::Vec<DIM, double>());
    for (const auto &neigh_id : this->support_.InfluenceNodeIds(eval_point_id)) {
    
        // Get the coordinates of the supporting field node.
        auto i = &neigh_id - &this->support_.InfluenceNodeIds(eval_point_id)[0];
        neigh_nodes[i] = field_nodes[neigh_id];
    }

    // Initialize radial basis function values and gradient.
    arma::vec rho(neighs_num + m);
    arma::mat drho(neighs_num + m, DIM);

    // Compute the Vandermonde matrix for the neighbor nodes.
    this->mbasis_->Compute(neigh_nodes);

    // Iterate over the support nodes to assemble the RPI matrix G.
    arma::mat G(neighs_num + m, neighs_num + m, arma::fill::zeros);
    for (const auto &neigh_id : this->support_.InfluenceNodeIds(eval_point_id)) {

        // Get position of the iterator.
        auto i = &neigh_id - &this->support_.InfluenceNodeIds(eval_point_id)[0];

        // Compute the values and gradient of the rbf for the current field node in the support domain.
        this->ComputeEnrichedRbf(field_nodes[neigh_id], neigh_nodes, this->support_.InfluenceNodeIds(eval_point_id), rho, drho);

        // Fill the G matrix with the rbf values.
        for (std::size_t j = 0; j != neighs_num; ++j) { G(j, i) = rho(j); }

        // Fill the G matrix with the monomial basis values.
        for (short j = 0; j != m; ++j) {
            G(i, neighs_num+j) = this->mbasis_->Basis(i)(j);
            G(neighs_num+j, i) = this->mbasis_->Basis(i)(j);
        }
    } // End of Iterate over the support nodes to assemble the RPI matrix G.

    // Compute the values and gradient of the rbf for the evaluation point.
    this->ComputeEnrichedRbf(eval_point, neigh_nodes, this->support_.InfluenceNodeIds(eval_point_id), rho, drho);

    // Get the RPI basis function values.
    arma::vec rpi_phi = arma::solve(G, rho);
    rpi_phi = rpi_phi.head_rows(neighs_num);
    this->phi_[0] = this->CastToEigen(rpi_phi);

    // Get the RPI basis function gradient.
    arma::mat rpi_phi_grad(neighs_num + m, DIM);
    for (int d = 0; d != DIM; ++d) {
        rpi_phi_grad.col(d) = arma::solve(G, drho.col(d));
    }
    rpi_phi_grad = rpi_phi_grad.head_rows(neighs_num);
    this->phi_grad_[0] = this->CastToEigen(rpi_phi_grad);

}



} // End of namespace CLOUDEA

#endif //CLOUDEA_APPROXIMANTS_RPI_TPP_
