
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

/**
   \file weight_factory.hpp
   \brief WeightFactory class header file.
   \author Konstantinos A. Mountris
   \date 25/06/2020
*/

#ifndef CLOUDEA_WEIGHT_FUNCTIONS_WEIGHT_FACTORY_HPP_
#define CLOUDEA_WEIGHT_FUNCTIONS_WEIGHT_FACTORY_HPP_


#include "CLOUDEA/engine/weight_functions/weight_function.hpp"
#include "CLOUDEA/engine/weight_functions/cubic_spline.hpp"
#include "CLOUDEA/engine/weight_functions/quartic_spline.hpp"
#include "CLOUDEA/engine/weight_functions/gaussian_rbf.hpp"
#include "CLOUDEA/engine/weight_functions/multiquadric_rbf.hpp"
#include "CLOUDEA/engine/weight_functions/polyharmonic_rbf.hpp"
#include "CLOUDEA/engine/utilities/logger.hpp"

#include <memory>
#include <string>


namespace CLOUDEA {

/** \addtogroup WeightFunctions \{ */

/**
 * \class WeightFactory
 * \author Konstantinos A. Mountris
 * \brief Factory for weight functions generation.
 * \tparam DIM The spatial dimensions of the generated weight function.
 */

template <short DIM>
class WeightFactory
{   
public:

    /**
     * \brief Creates a weight function according to the provided type.
     * \param [in] weight_type The type of the desired weight function.
     * \return [std::unique_ptr<WeightFunction>] Unique pointer to the created weight function.
     */
    inline static std::unique_ptr<WeightFunction<DIM>> Create(WeightType weight_type) {
        
        std::unique_ptr<WeightFunction<DIM>> weight_ptr;

        switch (weight_type) {
        case WeightType::cubic :
            weight_ptr = std::make_unique<CubicSpline<DIM>>();
            break;
        case WeightType::quartic :
            weight_ptr = std::make_unique<QuarticSpline<DIM>>();
            break;
        case WeightType::gaussian :
            weight_ptr = std::make_unique<GaussianRbf<DIM>>();
            break;
        case WeightType::multiquadric :
            weight_ptr = std::make_unique<MultiquadricRbf<DIM>>();
            break;
        case WeightType::polyharmonic :
            weight_ptr = std::make_unique<PolyharmonicRbf<DIM>>();
            break;
        default:
            std::string error_str = "Could not create weight function. Not supported type.";
            throw std::invalid_argument(Logger::Error(error_str, __FILE__, __LINE__));
            break;
        }

        return weight_ptr;
    }
};


/** \} End of Doxygen Groups */

} // End of namespace ELECTRA

#endif // CLOUDEA_WEIGHT_FUNCTIONS_WEIGHT_FACTORY_HPP_