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

#ifndef CLOUDEA_APPROXIMANTS_MKI_TPP_
#define CLOUDEA_APPROXIMANTS_MKI_TPP_

#include "CLOUDEA/engine/approximants/mki.hpp"

namespace CLOUDEA {

template <short DIM>
void Mki<DIM>::ComputeCorrelationFunction(const IMP::Vec<DIM, double> &eval_point, const std::vector<IMP::Vec<DIM, double>> &neigh_nodes, 
                                          const std::vector<int> &neigh_ids, arma::vec &mki_values, arma::mat &mki_grad)
{
    // Get the number of neighbor field nodes in the support domain.
    std::size_t neighs_num = neigh_nodes.size();

    // Reset values of the radial basis function and its gradient.
    mki_values = arma::vec(neighs_num, arma::fill::zeros);
    mki_grad   = arma::mat(neighs_num, DIM, arma::fill::zeros);

    // Iterate over neighbor nodes.
    for (std::size_t i = 0; i != neighs_num; ++i) {

        // Compute the weight function of neighbor node at eval_point.
        this->weight_->Compute(eval_point, neigh_nodes[i], this->support_.Radius(neigh_ids[i]), this->support_.DilateCoeff(neigh_ids[i]));

        // Compute the gaussian rbf value at the current neighbor.
        mki_values(i) = this->weight_->Val();

        // Compute the gradient of the gaussian rbf value at the current neighbor.
        for (std::size_t d = 0; d != DIM; ++d) {
            mki_grad(i, d) = this->weight_->Grad(d);
        }
    }

}


template <short DIM>
Mki<DIM>::Mki()
{
    this->type_ = MfreeType::mki;
}


template <short DIM>
Mki<DIM>::~Mki()
{}


template <short DIM>
void Mki<DIM>::Compute(const std::vector<IMP::Vec<DIM, double>> &eval_points,
                       const std::vector<IMP::Vec<DIM, double>> &field_nodes)
{
    std::size_t eval_points_num = eval_points.size();
    this->phi_.clear();       this->phi_.reserve(eval_points_num);
    this->phi_grad_.clear();  this->phi_grad_.reserve(eval_points_num);

    // Iterate over eval_points.
    for (const auto &point : eval_points) {

        // Get the index of the current evaluation point.
        auto eval_point_id = &point - &eval_points[0];

        // Number of neighbor field nodes to the evaluation point.
        auto neighs_num = this->support_.InfluenceNodeIds(eval_point_id).size();


        // Get the coordinates of the neighbor nodes in the support domain.
        std::vector<IMP::Vec<DIM, double>> neigh_nodes(neighs_num, IMP::Vec<DIM, double>());
        for (std::size_t i = 0; i != neighs_num; ++i) {
            // Get the coordinates of the supporting field node.
            auto neigh_id = this->support_.InfluenceNodeIds(eval_point_id)[i];
            neigh_nodes[i] = field_nodes[neigh_id];
        }

        // Initialize moving kriging function values and gradient.
        arma::vec rho(neighs_num);
        arma::mat drho(neighs_num, DIM);

        // Compute the monomial basis Vandermonde matrix for neighbors.
        this->mbasis_->Compute(neigh_nodes);
        
        // Compute R moment matrix.
        arma::mat R(neighs_num, neighs_num);
        for (std::size_t i = 0; i != neighs_num; ++i) {

            // Get position of the iterator.
            auto neigh_id = this->support_.InfluenceNodeIds(eval_point_id)[i];

            // Compute the values and gradient of the mki for the current field node in the support domain.
            this->ComputeCorrelationFunction(field_nodes[neigh_id], neigh_nodes, this->support_.InfluenceNodeIds(eval_point_id), rho, drho);

            // Fill the R matrix with the mki values.
            for (std::size_t j = 0; j != neighs_num; ++j) { R(i,j) = rho(j); }

        } // End Compute R moment matrix.

        // Compute inverse R and transpose P matrices.
        arma::mat I = arma::eye(neighs_num, neighs_num);
        arma::mat Rinv = arma::solve(R,I);
        arma::mat P = this->mbasis_->Vandermonde();
        arma::mat Pt = P.t();
        
        // Compute A and B matrices.
        arma::mat A = arma::solve(Pt*Rinv*P, Pt*Rinv);
        arma::mat B = Rinv*(I - P*A);

        // Compute monomial basis for evaluation point.
        this->mbasis_->Compute(point);

        // Compute the values and gradient of the mki for the evaluation point.
        this->ComputeCorrelationFunction(point, neigh_nodes, this->support_.InfluenceNodeIds(eval_point_id), rho, drho);

        // Get the basis function values and gradient on the support domain nodes.
        arma::vec mki_phi = (this->mbasis_->Basis(0)*A + rho.t()*B).t();
        this->phi_.emplace_back( this->CastToEigen(mki_phi) );

        arma::mat mki_phi_grad = (this->mbasis_->Grad(0)*A + drho.t()*B).t();
        this->phi_grad_.emplace_back( this->CastToEigen(mki_phi_grad) );

    }

}


template <short DIM>
void Mki<DIM>::Compute(const IMP::Vec<DIM, double> &eval_point, int eval_point_id,
                       const std::vector<IMP::Vec<DIM, double>> &field_nodes)
{
    // Number of neighbor field nodes to the evaluation point.
    std::size_t neighs_num = this->support_.InfluenceNodeIds(eval_point_id).size();

    // Reset MKI basis function and gradient to zero for a single evaluation point.
    this->phi_.clear();
    this->phi_.resize(1, Eigen::VectorXd::Zero(neighs_num));
    this->phi_grad_.clear();
    this->phi_grad_.resize(1, Eigen::MatrixXd::Zero(neighs_num, DIM));

    // Get the coordinates of the neighbor nodes in the support domain.
    std::vector<IMP::Vec<DIM, double>> neigh_nodes(neighs_num, IMP::Vec<DIM, double>());
    for (std::size_t i = 0; i != neighs_num; ++i) {
        // Get the coordinates of the supporting field node.
        auto neigh_id = this->support_.InfluenceNodeIds(eval_point_id)[i];
        neigh_nodes[i] = field_nodes[neigh_id];
    }

    // Initialize moving kriging function values and gradient.
    arma::vec rho(neighs_num);
    arma::mat drho(neighs_num, DIM);

    // Compute the monomial basis Vandermonde matrix for neighbors.
    this->mbasis_->Compute(neigh_nodes);

    // Compute R moment matrix.
    arma::mat R(neighs_num, neighs_num);
    for (const auto &neigh_id : this->support_.InfluenceNodeIds(eval_point_id)) {

        // Get position of the iterator.
        auto i = &neigh_id - &this->support_.InfluenceNodeIds(eval_point_id)[0];

        // Compute the values and gradient of the mki for the current field node in the support domain.
        this->ComputeCorrelationFunction(field_nodes[neigh_id], neigh_nodes, this->support_.InfluenceNodeIds(eval_point_id), rho, drho);

        // Fill the R matrix with the mki values.
        for (std::size_t j = 0; j != neighs_num; ++j) { R(i,j) = rho(j); }

    } // End Compute R moment matrix.

    // Compute inverse R and transpose P matrices.
    arma::mat I = arma::eye(neighs_num, neighs_num);
    arma::mat Rinv = arma::solve(R,I);
    arma::mat P = this->mbasis_->Vandermonde();
    arma::mat Pt = P.t();
    
    // Compute A and B matrices.
    arma::mat A = arma::solve(Pt*Rinv*P, Pt*Rinv);
    arma::mat B = Rinv*(I - P*A);

    // Compute monomial basis for evaluation point.
    this->mbasis_->Compute(eval_point);

    // Compute the values and gradient of the mki for the evaluation point.
    this->ComputeCorrelationFunction(eval_point, neigh_nodes, this->support_.InfluenceNodeIds(eval_point_id), rho, drho);

    // Get the basis function values and gradient on the support domain nodes.
    arma::vec mki_phi = (this->mbasis_->Basis(0)*A + rho.t()*B).t();
    this->phi_[0] = this->CastToEigen(mki_phi);

    arma::mat mki_phi_grad = (this->mbasis_->Grad(0)*A + drho.t()*B).t();
    this->phi_grad_[0] = this->CastToEigen(mki_phi_grad);

}



} // End of namespace CLOUDEA

#endif //CLOUDEA_APPROXIMANTS_MKI_TPP_
