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


#ifndef CLOUDEA_MONOMIAL_LINEAR_BASIS_TPP_
#define CLOUDEA_MONOMIAL_LINEAR_BASIS_TPP_

#include "CLOUDEA/engine/monomial/linear_basis.hpp"


namespace CLOUDEA {


template<short DIM>
LinearBasis<DIM>::LinearBasis()
{
    this->rank_ = DIM+1;
    this->type_ = MonomialType::linear;
}


template<short DIM>
LinearBasis<DIM>::~LinearBasis()
{}


template<short DIM>
void LinearBasis<DIM>::Compute(const IMP::Vec<DIM, double> &point)
{
    // Initialize the vandermonde matrix of the linear monomial basis for the point.
    this->vandermonde_ = arma::mat(1, this->rank_, arma::fill::zeros);

    // Initialize the linear monomial basis for the point.
    this->basis_.clear();
    this->basis_.resize(1, arma::rowvec(this->rank_, arma::fill::zeros));

    // Initialize the gradient of the linear monomial basis for the point.
    this->grad_.clear();
    this->grad_.resize(1, arma::mat(DIM, this->rank_, arma::fill::zeros));

    // Set linear basis and gradient.
    this->basis_[0](0) = 1.;
    for (short d = 0; d != DIM; ++d) {
        this->basis_[0](d+1) = point[d];
        this->grad_[0](d, d+1) = 1.;
    }

    // Assign the monomial basis to the Vandermonde matrix.
    this->vandermonde_.row(0) = this->basis_[0];
}


template<short DIM>
void LinearBasis<DIM>::Compute(const std::vector<IMP::Vec<DIM, double>> &points)
{
    // Initialize the vandermonde matrix of the linear monomial basis for the points.
    this->vandermonde_ = arma::mat(points.size(), this->Rank(), arma::fill::zeros);

    // Initialize the container of the linear monomial basis for the points.
    this->basis_.clear();
    this->basis_.resize(points.size(), arma::rowvec(this->Rank(), arma::fill::zeros));

    // Initialize the gradient of the linear monomial basis for the points.
    this->grad_.clear();
    this->grad_.resize(points.size(), arma::mat(DIM, this->Rank(), arma::fill::zeros));

    // Compute the basis and its gradient of the linear monomial basis for each point.
    for (const auto &point : points) {
        auto id = &point - &points[0];

        this->basis_[id](0) = 1.;
        for (short d = 0; d != DIM; ++d) {
            this->basis_[id](d+1) = point[d];
            this->grad_[id](d, d+1) = 1.;
        }

        // Assign the monomial basis to the Vandermonde matrix.
        this->vandermonde_.row(id) = this->basis_[id];
    }
}


} //End of namespace CLOUDEA

#endif //CLOUDEA_MONOMIAL_LINEAR_BASIS_TPP_