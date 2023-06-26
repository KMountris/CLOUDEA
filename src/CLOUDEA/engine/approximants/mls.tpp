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


#ifndef CLOUDEA_APPROXIMANTS_MLS_TPP_
#define CLOUDEA_APPROXIMANTS_MLS_TPP_


#include "CLOUDEA/engine/approximants/mls.hpp"

namespace CLOUDEA {


template<short DIM>
Mls<DIM>::Mls()
{
    // Set the approximation type to MLS.
    this->type_ = MfreeType::mls;
}


template<short DIM>
Mls<DIM>::~Mls()
{}


template<short DIM>
void Mls<DIM>::Compute(const std::vector<IMP::Vec<DIM, double>> &eval_points,
                       const std::vector<IMP::Vec<DIM, double>> &field_nodes)
{
    std::size_t eval_points_num = eval_points.size();
    this->phi_.clear();       this->phi_.reserve(eval_points_num);
    this->phi_grad_.clear();  this->phi_grad_.reserve(eval_points_num);

    // Get the rank of the monomial basis.
    short m = this->mbasis_->Rank();

    // Iterate over evaluation points.
    for (const auto &eval_point : eval_points) {

        // Get the index of the current evaluation eval_point.
        auto eval_point_id = &eval_point - &eval_points[0];
        
        // Number of neighbor field nodes to the evaluation eval_point.
        auto neighs_num = this->support_.InfluenceNodeIds(eval_point_id).size();

        // Initialize matrix A and gradient.
        arma::mat A(m, m, arma::fill::zeros);
        std::vector<arma::mat> dA(DIM, A);

        // Initialize matrix B and gradient.
        arma::mat B(m, neighs_num, arma::fill::zeros);
        std::vector<arma::mat> dB(DIM, B);

        // Compute A,B matrices and their gradients.
        for (const auto &neigh_id : this->support_.InfluenceNodeIds(eval_point_id)) {
        
            // Get position of the iterator.
            auto i = &neigh_id - &this->support_.InfluenceNodeIds(eval_point_id)[0];

            // Compute the weight function of the supporting field node at the evaluation eval_point.    
            this->weight_->Compute(eval_point, field_nodes[neigh_id], this->support_.Radius(neigh_id), this->support_.DilateCoeff(neigh_id));

            // Compute the monomial basis of the supporting field node.
            this->mbasis_->Compute(field_nodes[neigh_id]);

            A += this->weight_->Val() * ( this->mbasis_->Basis(0).t()*this->mbasis_->Basis(0) );
            B.col(i) = this->weight_->Val() * this->mbasis_->Basis(0).t();

            // Set weight function and update A matrix gradients.
            for (short d = 0; d != DIM; ++d) {
                dA[d] += this->weight_->Grad(d) * ( this->mbasis_->Basis(0).t()*this->mbasis_->Basis(0) );
                dB[d].col(i) += this->weight_->Grad(d) * this->mbasis_->Basis(0).t();
            }

        }

        // Apply correction in A for high rank monomial basis. 
        if (m > (DIM+1)) {
            for (short i = DIM+1; i != m; ++i) {
                A(i,i) += 1.e-7;
            }
        }

        // Compute the C matrix and its transpose.
        arma::mat C_t = arma::solve(A, B);
        arma::mat C = C_t.t();

        // Compute the monomial basis of the evaluation eval_point.
        this->mbasis_->Compute(eval_point);

        // Compute and store the basis function.
        arma::vec mls_phi = C * this->mbasis_->Basis(0).t();
        this->phi_.emplace_back( this->CastToEigen(mls_phi) );

        // Compute the basis function gradient.
        arma::mat mls_phi_grad = arma::mat(neighs_num, DIM, arma::fill::zeros);
        for (short d = 0; d != DIM; ++d) {

            arma::mat B1 = dB[d] - dA[d]*C_t;
            arma::mat C1 = ( arma::solve(A, B1) ).t();

            mls_phi_grad.col(d) = C*this->mbasis_->Grad(0).row(d).t() + C1*this->mbasis_->Basis(0).t();
        }

        // Store the basis function gradient.
        this->phi_grad_.emplace_back( this->CastToEigen(mls_phi_grad) );
    }
}


template<short DIM>
void Mls<DIM>::Compute(const IMP::Vec<DIM, double> &eval_point, int eval_point_id,
                       const std::vector<IMP::Vec<DIM, double>> &field_nodes)
{

    // Number of neighbor field nodes to the evaluation eval_point.
    std::size_t neighs_num = this->support_.InfluenceNodeIds(eval_point_id).size();

    // Reset MLS basis function and gradient to zero.
    this->phi_.clear();
    this->phi_.resize(1, Eigen::VectorXd::Zero(neighs_num));
    this->phi_grad_.clear();
    this->phi_grad_.resize(1, Eigen::MatrixXd::Zero(neighs_num, DIM));

    // Get the rank of the monomial basis.
    short m = this->mbasis_->Rank();

    // Initialize matrix A and gradient.
    arma::mat A(m, m, arma::fill::zeros);
    std::vector<arma::mat> dA(DIM, A);

    // Initialize matrix B and gradient.
    arma::mat B(m, neighs_num, arma::fill::zeros);
    std::vector<arma::mat> dB(DIM, B);

    // Compute A,B matrices and their gradients.
    for (const auto &neigh_id : this->support_.InfluenceNodeIds(eval_point_id)) {
        
        // Get position of the iterator.
        auto i = &neigh_id - &this->support_.InfluenceNodeIds(eval_point_id)[0];

        // Compute the weight function of the supporting field node at the evaluation eval_point.    
        this->weight_->Compute(eval_point, field_nodes[neigh_id], this->support_.Radius(neigh_id), this->support_.DilateCoeff(neigh_id));

        // Compute the monomial basis of the supporting field node.
        this->mbasis_->Compute(field_nodes[neigh_id]);

        A += this->weight_->Val() * ( this->mbasis_->Basis(0).t()*this->mbasis_->Basis(0) );
        B.col(i) = this->weight_->Val() * this->mbasis_->Basis(0).t();

        // Set weight function and update A matrix gradients.
        for (short d = 0; d != DIM; ++d) {
            dA[d] += this->weight_->Grad(d) * ( this->mbasis_->Basis(0).t()*this->mbasis_->Basis(0) );
            dB[d].col(i) += this->weight_->Grad(d) * this->mbasis_->Basis(0).t();
        }

    }

    // Apply correction in A for high rank monomial basis. 
    if (m > (DIM+1)) {
        for (short i = DIM+1; i != m; ++i) {
            A(i,i) += 1.e-7;
        }
    }

    // Compute the C matrix and its transpose.
    arma::mat C_trans = arma::solve(A, B);
    arma::mat C = C_trans.t();

    // Compute the monomial basis of the evaluation eval_point.
    this->mbasis_->Compute(eval_point);

    // Compute and store the basis function.
    arma::vec mls_phi = C * this->mbasis_->Basis(0).t(); 
    this->phi_[0] = this->CastToEigen(mls_phi);

    // Compute the basis function gradient.
    arma::mat mls_phi_grad = arma::mat(neighs_num, DIM, arma::fill::zeros);
    for (short d = 0; d != DIM; ++d) {

        arma::mat B1 = dB[d] - dA[d]*C_trans;
        arma::mat C1 = ( arma::solve(A, B1) ).t();

        mls_phi_grad.col(d) = C*this->mbasis_->Grad(0).row(d).t() + C1*this->mbasis_->Basis(0).t();
    }

    this->phi_grad_[0] = this->CastToEigen(mls_phi_grad);

}


} // End of namespace CLOUDEA

#endif //CLOUDEA_APPROXIMANTS_MLS_TPP_