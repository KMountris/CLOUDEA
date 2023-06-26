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
   \file fem_d1v2.hpp
   \brief FemD1v2 class header file.
   \author Konstantinos A. Mountris
   \date 11/07/2021
*/

#ifndef CLOUDEA_APPROXIMANTS_FEM_D1V2_HPP_
#define CLOUDEA_APPROXIMANTS_FEM_D1V2_HPP_

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
 * \class FemD1v2
 * \brief Class implemmenting a linear finite element in 1 dimension.
 */

template<short DIM>
class FemD1v2 : public BaseFem<DIM>
{

private:

    const FemType type_ = FemType::D1V2;              /**< The type of the element */


public:

    /**
     * \brief The FemD1v2 constructor.
     */
    FemD1v2();


    /**
     * \brief The FemD1v2 destructor.
     */
    virtual ~FemD1v2();


    /**
     * \brief Set the quadrature rule of the element.
     * \param [in] quad_points_num The number of the quadrature points.
     * \return [void]
     */
    virtual void SetQuadrature(std::initializer_list<short> quad_points_num);


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
     * \brief Compute the derivatives of the element's shape function in the natural coordinate (ξ) for each quadrature point.
     * \return [void]
     */
    virtual void ComputeDerivsNatural();


    /**
     * \brief Compute the derivatives of the element's shape function in the physical coordinate (x) for each quadrature point.
     * \return [void]
     */
    virtual void ComputeDerivs();


    /**
     * \brief Get the type of the element.
     * \return [FemType::D1V2] The type of the element.
     */
    inline virtual FemType Type() const { return this->type_; }


};


/** \} End of Doxygen Groups */
/** \} End of Doxygen Groups */

} // End of namespace CLOUDEA

#endif //CLOUDEA_APPROXIMANTS_FEM_D1V2_HPP_

#include "CLOUDEA/engine/approximants/fem_d1v2.tpp"