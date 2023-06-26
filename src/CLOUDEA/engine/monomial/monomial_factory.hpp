
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
   \file monomial_factory.hpp
   \brief MonomialFactory class header file.
   \author Konstantinos A. Mountris
   \date 25/06/2020
*/

#ifndef CLOUDEA_MONOMIAL_MONOMIAL_FACTORY_HPP_
#define CLOUDEA_MONOMIAL_MONOMIAL_FACTORY_HPP_


#include "CLOUDEA/engine/monomial/monomial_basis.hpp"
#include "CLOUDEA/engine/monomial/linear_basis.hpp"
#include "CLOUDEA/engine/monomial/quadratic_basis.hpp"
#include "CLOUDEA/engine/monomial/cubic_basis.hpp"
#include "CLOUDEA/engine/utilities/logger.hpp"

#include <memory>
#include <string>


namespace CLOUDEA {

/** \addtogroup Monomial \{ */

/**
 * \class MonomialFactory
 * \author Konstantinos A. Mountris
 * \brief Factory for monomial basis function generation.
 * \tparam DIM The spatial dimensions of the generated monomial basis function.
 */

template <short DIM>
class MonomialFactory
{   
    
public:

    /**
     * \brief Creates a monomial basis function according to the provided type.
     * \param [in] weight_type The type of the desired monomial basis function.
     * \return [std::unique_ptr<MonomialBasis>] Unique pointer to the created monomial basis function.
     */
    inline static std::unique_ptr<MonomialBasis<DIM>> Create(MonomialType monomial_type) {
        
        std::unique_ptr<MonomialBasis<DIM>> monomial_ptr;

        switch (monomial_type) {
        case MonomialType::linear :
            monomial_ptr = std::make_unique<LinearBasis<DIM>>();
            break;
        case MonomialType::quadratic :
            monomial_ptr = std::make_unique<QuadraticBasis<DIM>>();
            break;
        case MonomialType::cubic :
            monomial_ptr = std::make_unique<CubicBasis<DIM>>();
            break;
        default:
            std::string error_str = "Could not create monomial basis function. Not supported type.";
            throw std::invalid_argument(Logger::Error(error_str, __FILE__, __LINE__));
            break;
        }

        return monomial_ptr;
    }
};


/** \} End of Doxygen Groups */

} // End of namespace ELECTRA

#endif // CLOUDEA_WEIGHT_FUNCTIONS_WEIGHT_FACTORY_HPP_