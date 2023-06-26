 /*
 * CLOUDEA - Software for solving PDEs using explicit methods.
 * Copyright (C) 2017  <Konstantinos A. Mountris> <konstantinos.mountris\gmail.com>
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

/**
   \file quartic_spline.hpp
   \brief QuarticSpline class header file.
   \author Konstantinos A. Mountris
   \date 25/06/2020
*/

#ifndef CLOUDEA_WEIGHT_FUNCTIONS_QUARTIC_SPLINE_HPP_
#define CLOUDEA_WEIGHT_FUNCTIONS_QUARTIC_SPLINE_HPP_

#include "CLOUDEA/engine/weight_functions/weight_function.hpp"
#include "CLOUDEA/engine/utilities/logger.hpp"

#include <cmath>
#include <limits>

namespace CLOUDEA {

/** \addtogroup WeightFunctions \{ */


/**
 * \class QuarticSpline
 * \brief Class implemmenting a quartic spline weight function.
 * \tparam DIM The spatial dimension of the quartic spline.
 */
template <short DIM>
class QuarticSpline : public WeightFunction<DIM>
{

public:

    inline QuarticSpline();


    inline ~QuarticSpline();


    /**
     * \brief Compute the quartic spline. 
     * \param [in] point The point to evaluate the quartic spline.
     * \param [in] center The center of the quartic spline. 
     * \param [in] radius The radius of the quartic spline.
     * \param [in] dilate_coeff The dilatation coefficient for the quartic spline's radius.
     * \return [void]
     */
    void Compute(const IMP::Vec<DIM, double> &point, const IMP::Vec<DIM, double> &center, double radius, double dilate_coeff);

};

/** \} End of Doxygen Groups */

} // End of namespace CLOUDEA

#endif //CLOUDEA_WEIGHT_FUNCTIONS_QUARTIC_SPLINE_HPP_

#include "CLOUDEA/engine/weight_functions/quartic_spline.tpp"
