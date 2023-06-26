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
   \file quadrature_rule.hpp
   \brief QuadratureRule class header file.
   \author Konstantinos A. Mountris
   \date 11/07/2021
*/

#ifndef CLOUDEA_QUADRATURE_QUADRATURE_RULE_HPP_
#define CLOUDEA_QUADRATURE_QUADRATURE_RULE_HPP_


#include "CLOUDEA/engine/utilities/logger.hpp"

#include <IMP/Vectors>

#include <vector>
#include <stdexcept>
#include <exception>


namespace CLOUDEA {

/** \addtogroup QuadratureRule \{ */


/**
 * \class QuadratureRule
 * \brief Class implemmenting quadrature rules for various element types.
 *
 * Quadrature is used for the evaluation of integrals that arise in the Galerkin weak form of Partial Differential Equations.
 * It is required in the Finite Elements and global weak form Meshless (e.g. Element Free Galerkin) methods. Quadrature rules
 * in the natural coordinates system are available in 1-, 2-, and 3 dimensions for line, square, cube, triangle, and tetrahedron elements.
 *
 * \section  sec_quad_on_line Quadrature on line element
 * For line elements the Gauss Legendre quadrature rules described in \cite abbott2005 are used. The available quadrature accuracy is up to order 12.
 * The number of quadrature points and the corresponding order for the triangle quadrature rule is given in the table below:
 *
 * Points | Order
 * ------ | ------
 *   1    |   1
 *   2    |   2
 *   3    |   3
 *   4    |   4
 *   5    |   5
 *   6    |   6
 *   7    |   7
 *   8    |   8
 *   9    |   9
 *  10    |  10
 *  11    |  11
 *  12    |  12
 *
 * \section  sec_quad_on_square Quadrature on square element
 * For square elements the quadrature rules described in \cite abbott2005 are used for both X and Y axis. The available quadrature accuracy per axis is up to order 12.
 * The minimum quadrature order of the two axis is chosen as the quadrature order of the element. The number of quadrature points and the corresponding order for each axis
 * of the square element is given in Section \ref sec_quad_on_line.
 *
 * \section  sec_quad_on_triangle Quadrature on triangle element
 * For triangle elements the quadrature rules described in \cite dunavant1985 are used. The available quadrature accuracy is up to order 17.
 * The number of quadrature points and the corresponding order for the triangle quadrature rule is given in the table below:
 *
 * Points | Order
 * ------ | ------
 *   1    |   1
 *   3    |   2
 *   4    |   3
 *   6    |   4
 *   7    |   5
 *  12    |   6
 *  13    |   7
 *  16    |   8
 *  19    |   9
 *  25    |  10
 *  27    |  11
 *  33    |  12
 *  37    |  13
 *  42    |  14
 *  48    |  15
 *  52    |  16
 *  61    |  17
 *
 * \section  sec_quad_on_tetrahedron Quadrature on tetrahedron element
 * For tetrahedral elements the quadrature rules described in \cite keast1986 are used. The available quadrature accuracy is up to order 17.
 * The number of quadrature points and the corresponding order for the tetrahedron quadrature rule is given in the table below:
 *
 * Points | Order
 * ------ | ------
 *   1    |   0
 *   4    |   1
 *   5    |   2
 *  10    |   3
 *  11    |   4
 *  15    |   5
 *  24    |   6
 *  31    |   7
 *  45    |   8
 *
 * \section  sec_quad_on_cube Quadrature on cubic element
 * For cubic elements the quadrature rules described in \cite abbott2005 are used for the X, Y and Z axis. The available quadrature accuracy per axis is up to order 12.
 * The minimum quadrature order of the three axis is chosen as the quadrature order of the element. The number of quadrature points and the corresponding order for each axis
 * of the cubic element is given in Section \ref sec_quad_on_line.
 *
 * \tparam DIM The dimensions number of the quadrature points coordinates. Supported: 1 | 2 | 3
 */

template<short DIM>
class QuadratureRule
{


private:

    std::vector<IMP::Vec<DIM, double>> points_;   /**< The quadrature points natural coordinates */

    std::vector<double> weights_;                 /**< The weights of the quadrature points */

    short order_;                                 /**< The order of the quadrature */


public:

    /**
     * \brief The QuadratureRule constructor.
     */
    inline QuadratureRule();


    /**
     * \brief The QuadratureRule destructor.
     */
    inline virtual ~QuadratureRule();


    /**
     * \brief Set the quadrature points and their weights for a line in natural coordinates system using the rules described in \cite abbott2005.
     * \param [in] quad_points_num The number of the required quadrature points.
     * \return [void]
     */
    inline void SetForLineNat(short quad_points_num);


    /**
     * \brief Set the quadrature points and their weights for a square in natural coordinates system using the rules described in \cite abbott2005.
     * The quadrature order in this case corresponds to the minimum quadrature order over the X, Y axis.
     * \param [in] quad_points_num_x The number of the required quadrature points over the X-axis.
     * \param [in] quad_points_num_y The number of the required quadrature points over the Y-axis.
     */
    inline void SetForQuadrilateralNat(short quad_points_num_x, short quad_points_num_y);


    /**
     * \brief Set the quadrature points and their weights for a triangle in natural coordinates system using the rules described in \cite dunavant1985.
     * \param [in] quad_points_num The number of the required quadrature points.
     * \return [void]
     */
    inline void SetForTriangleNat(short quad_points_num);


    /**
     * \brief Set the quadrature points and their weights for a tetrahedron in natural coordinates system using the rules described in \cite keast1986.
     * \param [in] quad_points_num The number of the required quadrature points.
     * \return [void]
     */
    inline void SetForTetrahedronNat(short quad_points_num);


    /**
     * \brief Set the quadrature points and their weights for a cube in natural coordinates system using the rules described in \cite abbott2005.
     * The quadrature order in this case corresponds to the minimum quadrature order over the X, Y, Z axis.
     * \param [in] quad_points_num_x The number of the required quadrature points over the X-axis.
     * \param [in] quad_points_num_y The number of the required quadrature points over the Y-axis.
     * \param [in] quad_points_num_z The number of the required quadrature points over the Z-axis.
     */
    inline void SetForHexahedronNat(short quad_points_num_x, short quad_points_num_y, short quad_points_num_z);


    /**
     * \brief Get the order of the quadrature.
     * \return [short int] The order of the quadrature.
     */
    inline short Order() const { return this->order_; }


    /**
     * \brief Get the number of the quadrature points.
     * \return [int] The number of the quadrature points.
     */
    inline int PointsNum() const { return static_cast<int>(this->points_.size()); }


    /**
     * \brief Get the natural coordinates of the quadrature points.
     * \return [const std::vector<IMP::Vec<DIM, double>>&] The natural coordinates of the quadrature points.
     */
    inline const std::vector<IMP::Vec<DIM, double>> & Points() const { return this->points_; }


    /**
     * \brief Get the natural coordinates of the quadrature point with index id.
     * Fast access, no range check.
     * \param [in] id The index of the quadrature point.
     * \return [const IMP::Vec<DIM, double>&] The natural coordinates of the quadrature points.
     */
    inline const IMP::Vec<DIM, double> & Points(std::size_t id) const { return this->points_[id]; }


    /**
     * \brief Get the natural coordinates of the quadrature point with index id.
     * Slow access with range check.
     * \param [in] id The index of the quadrature point.
     * \return [const IMP::Vec<DIM, double>&] The natural coordinates of the quadrature point with index id.
     */
    inline const IMP::Vec<DIM, double> & PointsAt(std::size_t id) const { return this->points_.at(id); }


    /**
     * \brief Get the weights of the quadrature points.
     * \return [const std::vector<double>&] The weights of the quadrature points.
     */
    inline const std::vector<double> & Weights() const { return this->weights_; }


    /**
     * \brief Get the weights of the quadrature point with index id.
     * Fast access, no range check.
     * \param [in] id The index of the quadrature point.
     * \return [double] The weights of the quadrature point with index id.
     */
    inline double Weights(std::size_t id) const { return this->weights_[id]; }


    /**
     * \brief Get the weights of the quadrature point with index id.
     * Slow access with range check
     * \param [in] id The index of the quadrature point.
     * \return [double] The weights of the quadrature point with index id.
     */
    inline double WeightsAt(std::size_t id) const { return this->weights_.at(id); }

};


/** \} End of Doxygen Groups */

} // End of namespace CLOUDEA


#include "CLOUDEA/engine/quadrature/quadrature_rule.tpp"


#endif //CLOUDEA_QUADRATURE_QUADRATURE_RULE_HPP_
