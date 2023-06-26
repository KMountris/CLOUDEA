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


#ifndef CLOUDEA_WEIGHT_FUNCTIONS_MULTIQUADRIC_RBF_TPP_
#define CLOUDEA_WEIGHT_FUNCTIONS_MULTIQUADRIC_RBF_TPP_


#include "CLOUDEA/engine/weight_functions/multiquadric_rbf.hpp"

namespace CLOUDEA {


template<short DIM>
MultiquadricRbf<DIM>::MultiquadricRbf()
{
    this->type_ = WeightType::multiquadric;
}


template<short DIM>
MultiquadricRbf<DIM>::~MultiquadricRbf()
{}


template<short DIM>
void MultiquadricRbf<DIM>::Compute(const IMP::Vec<DIM, double> &point, const IMP::Vec<DIM, double> &center, double radius, double dilate_coeff)
{
    // Dilatation coefficient not applicable for multiquadric radial basis function.
    dilate_coeff = 1.0;
    radius *= dilate_coeff;

    // Squared distance of the evaluation point from the support center point of the multiquadric radial basis function.
    double r2 = point.Distance2(center);

    // Compute shape parameter.
    double a = this->theta_ * radius;

    // Compute the multiquadric radial basis function value w.
    this->val_ = std::pow(r2 + a*a, this->beta_);

    // Compute 1/r * dw/dr.
    double dval_dr = 2*this->beta_ * std::pow(r2 + a*a, this->beta_-1.);

    // Compute the multiquadric radial basis function gradient.
    this->grad_.SetZero();
    for (short d = 0; d != DIM; ++d) {
        // Compute dw/dx = (1/r * dw/dr) * (dr/dx * r)
        this->grad_[d] = dval_dr * (point[d]-center[d]);
    }
}


} // End of namespace CLOUDEA

#endif //CLOUDEA_WEIGHT_FUNCTIONS_MULTIQUADRIC_RBF_TPP_
