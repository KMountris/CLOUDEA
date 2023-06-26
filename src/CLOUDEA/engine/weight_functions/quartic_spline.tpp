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


#ifndef CLOUDEA_WEIGHT_FUNCTIONS_QUARTIC_SPLINE_TPP_
#define CLOUDEA_WEIGHT_FUNCTIONS_QUARTIC_SPLINE_TPP_


#include "CLOUDEA/engine/weight_functions/quartic_spline.hpp"

namespace CLOUDEA {


template<short DIM>
QuarticSpline<DIM>::QuarticSpline()
{
    this->type_ = WeightType::quartic;
}


template<short DIM>
QuarticSpline<DIM>::~QuarticSpline()
{}


template<short DIM>
void QuarticSpline<DIM>::Compute(const IMP::Vec<DIM, double> &point, const IMP::Vec<DIM, double> &center, double radius, double dilate_coeff)
{
    double dilated_radius = dilate_coeff * radius;

    // Normalized distance of the evaluation point from the quartic spline center point.
    double r = std::sqrt(point.Distance2(center))  / dilated_radius;

    // Compute quartic spline value w.
    this->val_ = 0.;
    this->grad_.SetZero();
    if (r < 1.0 + 2*std::numeric_limits<double>::epsilon()) {
        this->val_ = 1. - 6.*r*r + 8.*r*r*r - 3.*r*r*r*r;

        // Compute 1/r * dw/dr.
        double dval_dr = -12 + 24*r - 12*r*r;

        // Compute the quartic spline gradient.
        for (short d = 0; d != DIM; ++d) {
            // Compute dw/dx = (1/r * dw/dr) * (dr/dx * r)
            this->grad_[d] = dval_dr * (point[d]-center[d]) / (dilated_radius*dilated_radius);
        }
    }
    
}


} // End of namespace CLOUDEA

#endif //CLOUDEA_WEIGHT_FUNCTIONS_QUARTIC_SPLINE_TPP_
