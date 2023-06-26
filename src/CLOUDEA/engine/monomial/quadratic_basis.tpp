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


#ifndef CLOUDEA_MONOMIAL_QUADRATIC_BASIS_TPP_
#define CLOUDEA_MONOMIAL_QUADRATIC_BASIS_TPP_

#include "CLOUDEA/engine/monomial/quadratic_basis.hpp"


namespace CLOUDEA {


template<short DIM>
QuadraticBasis<DIM>::QuadraticBasis()
{
    if (DIM == 1) { this->rank_ = 3; }
    else if (DIM == 2) { this->rank_ = 6; }
    else if (DIM == 3) { this->rank_ = 10; }
        
    this->type_ = MonomialType::quadratic;
}


template<short DIM>
QuadraticBasis<DIM>::~QuadraticBasis()
{}


template<short DIM>
void QuadraticBasis<DIM>::Compute(const IMP::Vec<DIM, double> &point)
{
    // Initialize the vandermonde matrix of the quadratic monomial basis for the point.
    this->vandermonde_ = arma::mat(1, this->rank_, arma::fill::zeros);

    // Initialize the cubic quadratic basis for the point.
    this->basis_.clear();
    this->basis_.resize(1, arma::rowvec(this->rank_, arma::fill::zeros));

    // Initialize the gradient of the cubic quadratic basis for the point.
    this->grad_.clear();
    this->grad_.resize(1, arma::mat(DIM, this->rank_, arma::fill::zeros));

    // Set quadratic basis and gradient.
    if (DIM == 1) {
        double x = point[0];
        this->basis_[0] = { 1.,  x,  x*x };  

        this->grad_[0]  = {{ 0., 1., 2.*x }};
    }
    else if (DIM == 2) {
        double x = point[0];
        double y = point[1];

        this->basis_[0] = { 1., x,  y,  x*x,  y*y, x*y };

        this->grad_[0] = {{ 0., 1., 0., 2.*x,   0.,   y },
                          { 0., 0., 1.,   0., 2.*y,   x }};
    }
    else if (DIM == 3) {
        double x = point[0];
        double y = point[1];
        double z = point[2];

        this->basis_[0] = { 1.,  x,  y,  z,  x*x,  y*y,  z*z, x*y, y*z, x*z };

        this->grad_[0] = {{ 0., 1., 0., 0., 2.*x,   0.,   0.,   y,  0.,   z },
                          { 0., 0., 1., 0.,   0.,  2*y,   0.,   x,   z,  0. },
                          { 0., 0., 0., 1.,   0.,   0., 2.*z,  0.,   y,   x }};
    }

    // Assign the monomial basis to the Vandermonde matrix.
    this->vandermonde_.row(0) = this->basis_[0];
}


template<short DIM>
void QuadraticBasis<DIM>::Compute(const std::vector<IMP::Vec<DIM, double>> &points)
{
    // Initialize the vandermonde matrix of the quadratic monomial basis for the points.
    this->vandermonde_ = arma::mat(points.size(), this->Rank(), arma::fill::zeros);

    // Initialize the container of the quadratic monomial basis for the points.
    this->basis_.clear();
    this->basis_.resize(points.size(), arma::rowvec(this->Rank(), arma::fill::zeros));

    // Initialize the gradient of the quadratic monomial basis for the points.
    this->grad_.clear();
    this->grad_.resize(points.size(), arma::mat(DIM, this->Rank(), arma::fill::zeros));

    // Compute the basis and its gradient of the quadratic monomial basis for each point.
    for (const auto &point : points) {
        auto id = &point - &points[0];

        if (DIM == 1) {
            double x = point[0];

            this->basis_[id] = { 1.,  x,  x*x };  

            this->grad_[id] = {{ 0., 1., 2.*x }};
        }
        else if (DIM == 2) {
            double x = point[0];
            double y = point[1];

            this->basis_[id] = { 1., x,  y,  x*x,  y*y, x*y };

            this->grad_[id] = {{ 0., 1., 0., 2.*x,   0.,   y },
                               { 0., 0., 1.,   0., 2.*y,   x }};
        }
        else if (DIM == 3) {
            double x = point[0];
            double y = point[1];
            double z = point[2];

            this->basis_[id] = { 1.,  x,  y,  z,  x*x,  y*y,  z*z, x*y, y*z, x*z };

            this->grad_[id] = {{ 0., 1., 0., 0., 2.*x,   0.,   0.,   y,  0.,   z },
                               { 0., 0., 1., 0.,   0.,  2*y,   0.,   x,   z,  0. },
                               { 0., 0., 0., 1.,   0.,   0., 2.*z,  0.,   y,   x }};
        }

        // Assign the monomial basis to the Vandermonde matrix.
        this->vandermonde_.row(id) = this->basis_[id];
    }
}


} //End of namespace CLOUDEA

#endif //CLOUDEA_MONOMIAL_QUADRATIC_BASIS_TPP_