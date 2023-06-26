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


/*!
   \file fem_d3v8.hpp
   \brief FemD3v8 class header file.
   \author Konstantinos A. Mountris
   \date 05/05/2019
*/

#ifndef CLOUDEA_APPROXIMANTS_FEM_D3V8_HPP_
#define CLOUDEA_APPROXIMANTS_FEM_D3V8_HPP_

#include "CLOUDEA/engine/approximants/base_fem.hpp"
#include "CLOUDEA/engine/utilities/logger.hpp"

#include <vector>
#include <iostream>
#include <iomanip>
#include <iterator>
#include <string>
#include <initializer_list>

namespace CLOUDEA {

/** \addtogroup Approximants \{ */
/** \addtogroup Finite-Elements \{ */

/**
 * \class FemD3v8
 * \brief Class implemmenting a trilinear hexahedral finite element in 3 dimensions.
 */
template<short DIM>
class FemD3v8 : public BaseFem<DIM>
{

private:

    const FemType type_ = FemType::D3V8;                /**< The type of the element */

public:

    /**
     * \brief The FemD3v8 constructor.
     */
    FemD3v8();


    /**
     * \brief The FemD3v8 destructor.
     */
    virtual ~FemD3v8();


    /**
     * \brief Set the quadrature rule of the element.
     * quad_points_per_axis The number of the quadrature points per axis X, Y, Z for the element.
     * \return [void]
     */
    virtual void SetQuadrature(std::initializer_list<short> quad_points_per_axis);


    /**
     * \brief Compute the Jacobian of the element for each quadrature point.
     * Internally computes and stores in addition the determinant and the inverse of the Jacobian.
     * \param [in] nodes The nodes composing the element.
     */
    virtual void ComputeJacobians(const std::vector<IMP::Vec<DIM, double> > &nodes);


    /**
     * \brief Compute the shape function values of the element for each quadrature point.
     * \return [void]
     */
    virtual void ComputeShapeFunctions();


    /**
     * \brief Compute the derivatives of the element's shape function in the natural coordinates (ξ, η, ζ) for each quadrature point.
     * \return [void]
     */
    virtual void ComputeDerivsNatural();


    /**
     * \brief Compute the derivatives of the element's shape function in the physical coordinates (x, y, z) for each quadrature point.
     * \return [void]
     */
    virtual void ComputeDerivs();


    /**
     * \brief Get the type of the element.
     * \return [FemType::D3V8] The type of the element.
     */
    inline virtual FemType Type() const { return this->type_; }

};


/** \} End of Doxygen Groups */
/** \} End of Doxygen Groups */

} // End of namespace CLOUDEA

#include "CLOUDEA/engine/approximants/fem_d3v8.tpp"

#endif //CLOUDEA_APPROXIMANTS_FEM_D3V8_HPP_