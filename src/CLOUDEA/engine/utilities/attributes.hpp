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
   \file attributes.hpp
   \brief Collection of enumerations and types of various attributes used in CLOUDEA.
   \author Konstantinos A. Mountris
   \date 22/05/2019
*/

#ifndef CLOUDEA_UTILITIES_ATTRIBUTES_HPP_
#define CLOUDEA_UTILITIES_ATTRIBUTES_HPP_


namespace CLOUDEA {

/** \addtogroup Utilities \{ */



/**
 * \enum Collection of radial basis function types.
 * \author Konstantinos A. Mountris
 */
enum class RbfType { unknown,          /**< Unknown radial basis function type */
                     multiquadric,     /**< Multiquadratic radial basis function type */
                     gaussian,         /**< Gaussian radial basis function type */
                     polyharmonic      /**< Polyharmonic radial basis function type */
                   };



                              




/*! \} End of Doxygen Groups*/

} //end of namespace CLOUDEA

#endif //CLOUDEA_UTILITIES_ATTRIBUTES_HPP_