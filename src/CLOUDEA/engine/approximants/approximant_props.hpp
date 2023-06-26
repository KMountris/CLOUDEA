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
   \file approximant_props.hpp
   \brief Collection of approximant properties header file.
   \author Konstantinos A. Mountris
   \date 11/06/2021
*/

#ifndef CLOUDEA_APPROXIMANTS_APPROXIMANT_PROPS_HPP_
#define CLOUDEA_APPROXIMANTS_APPROXIMANT_PROPS_HPP_


namespace CLOUDEA {

/** \addtogroup Approximant-Properties \{ */


/**
 * \enum Enumeration for aproximant type declaration.
 * \author Konstantinos A. Mountris
 */
enum class ApproxType { fem,    /**< Finite Element Method approximation type */
                        fpm,    /**< Fragile Points Method approximation type */
                        efg,    /**< Element Free Galerkin Method approximation type */
                        mcm     /**< Mixed Collocation Method approximation type */
                      };

/**
 * \enum Type of finite element aproximants.
 * \author Konstantinos A. Mountris
 */
enum class FemType { UNKNOWN,   /**< Unknown element */
                     D1V2,      /**< Linear 1D line element with two vertices */
                     D2V3,      /**< Linear 2D triangle element with three vertices */
                     D2V4,      /**< Bilinear 2D quadrilateral element with four vertices */
                     D3V4,      /**< Linear 3D tetrahedral element with four vertices */
                     D3V8       /**< Trilinear 3D hexahedral element with eight vertices */
                   };


/**
 * \enum Type of meshfree approximants.
 * \author Konstantinos A. Mountris
 */
enum class MfreeType { unknown, /**< Unknown meshfree approximant type */
                       mls,     /**< Moving least squares meshfree approximant type */
                       rpi,     /**< Radial point interpolation meshfree approximant type */
                       mki      /**< Moving Kriging interpolation meshfree approximant type */
                     };



/** \} End of Doxygen Groups */

} // End of namespace CLOUDEA

#endif //CLOUDEA_APPROXIMANTS_APPROXIMANT_PROPS_HPP_