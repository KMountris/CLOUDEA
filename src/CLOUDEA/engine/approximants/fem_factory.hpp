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
   \file fem_factory.hpp
   \brief FemFactory class header file.
   \author Konstantinos A. Mountris
   \date 11/07/2021
*/

#ifndef CLOUDEA_APPROXIMANTS_FEM_FACTORY_HPP_
#define CLOUDEA_APPROXIMANTS_FEM_FACTORY_HPP_

#include "CLOUDEA/engine/approximants/approximant_props.hpp"
#include "CLOUDEA/engine/approximants/base_fem.hpp"
#include "CLOUDEA/engine/approximants/fem_d1v2.hpp"
#include "CLOUDEA/engine/approximants/fem_d2v3.hpp"
#include "CLOUDEA/engine/approximants/fem_d2v4.hpp"
#include "CLOUDEA/engine/approximants/fem_d3v4.hpp"
#include "CLOUDEA/engine/approximants/fem_d3v8.hpp"

#include <IMP/IMP>

#include <memory>
#include <string>


namespace CLOUDEA {

/** \addtogroup Approximants \{ */
/** \addtogroup Finite-Elements \{ */

/**
 * \class FemFactory
 * \brief Class implementing a factory for finite element generation.
 * \tparam DIM dimension of the generated finite element.
 * \author Konstantinos A. Mountris
 */

template<short DIM>
class FemFactory
{
public:

    /**
     * \brief Creates a finite element according to the provided type.
     * \param [in] ap_model_type The type of the desired finite element.
     * \return [std::unique_ptr<BaseFem>] Unique pointer to the created finite element.
     */
    inline static std::unique_ptr<BaseFem<DIM>> Create(IMP::CellShape cell_type) {

        std::unique_ptr<BaseFem<DIM>> fem_ptr;

        switch (cell_type) {
        case IMP::CellShape::edge :
            fem_ptr = std::make_unique<FemD1v2<DIM>>();
            break;
        case IMP::CellShape::tri :
            fem_ptr = std::make_unique<FemD2v3<DIM>>();
            break;
        case IMP::CellShape::quad :
            fem_ptr = std::make_unique<FemD2v4<DIM>>();
            break;
        case IMP::CellShape::tet :
            fem_ptr = std::make_unique<FemD3v4<DIM>>();
            break;
        case IMP::CellShape::hex :
            fem_ptr = std::make_unique<FemD3v8<DIM>>();
            break;
        default:
            std::string error_str = Logger::Error("Could not create finite element. Not supported type.");
            throw std::invalid_argument(error_str);
            break;
        }

        return fem_ptr;
    }
};


/** \} End of Doxygen Groups */
/** \} End of Doxygen Groups */

} // End of namespace CLOUDEA

#endif // CLOUDEA_APPROXIMANTS_FEM_FACTORY_HPP_