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
   \file cubic_basis.hpp
   \brief CubicBasis class header file.
   \author Konstantinos A. Mountris
   \date 25/06/2020
*/

#ifndef CLOUDEA_MONOMIAL_CUBIC_BASIS_HPP_
#define CLOUDEA_MONOMIAL_CUBIC_BASIS_HPP_

#include "CLOUDEA/engine/monomial/monomial_basis.hpp"


namespace CLOUDEA {

/** \addtogroup Monomial \{ */

/**
 * \class CubicBasis
 * \author Konstantinos A. Mountris
 * \brief Class implemmenting a cubic monomial basis.
 * \tparam DIM The spatial dimension of the cubic monomial basis.
 */
template <short DIM>
class CubicBasis : public MonomialBasis<DIM>
{

public:

    /**
    * \brief The CubicBasis default constructor.
    */
    inline CubicBasis();


    /**
    * \brief The CubicBasis destructor.
    */
    inline virtual ~CubicBasis();
    
    
    /**
    * \brief Compute the cubic monomial basis for a given point. 
    * \param [in] point The coordinates of the point.
    * \return [void]
    */
    inline void Compute(const IMP::Vec<DIM, double> &point);
  
  
    /**
    * \brief Compute the cubic monomial basis for a vector of points. 
    * \param point The coordinates of the points in the vector.
    * \return [void]
    */
    inline void Compute(const std::vector<IMP::Vec<DIM, double>> &points);

};

} //End of namespace CLOUDEA

#endif //CLOUDEA_MONOMIAL_CUBIC_BASIS_HPP_

#include "CLOUDEA/engine/monomial/cubic_basis.tpp"