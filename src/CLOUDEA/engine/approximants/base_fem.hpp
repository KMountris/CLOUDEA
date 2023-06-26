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
   \file base_fem.hpp
   \brief BaseFem class header file.
   \author Konstantinos A. Mountris
   \date 11/07/2021
*/

#ifndef CLOUDEA_APPROXIMANTS_BASE_FEM_HPP_
#define CLOUDEA_APPROXIMANTS_BASE_FEM_HPP_


#include "CLOUDEA/engine/approximants/approximant_props.hpp"
#include "CLOUDEA/engine/quadrature/quadrature_rule.hpp"

#include <IMP/Vectors>

#include <Eigen/Dense>


namespace CLOUDEA {

/** \addtogroup Approximants \{ */
/** \addtogroup Finite-Elements \{ */

/**
 * \class BaseFem
 * \brief Class implemmenting an abstract class for a generic finite element.
 * \tparam DIM The dimensions of the physical space that the element lies in.
 */
template<short DIM>
class BaseFem
{
private:

    QuadratureRule<DIM> quadrature_;              /*!< The quadrature of the element */

    std::vector<Eigen::MatrixXd> jacobians_;          /*!< The jacobian of the element at each quadrature point */

    std::vector<Eigen::MatrixXd> inv_jacobians_;       /*!< The inverse of the jacobian of the element at each quadrature point */

    std::vector<double> det_jacobians_;         /*!< The determinant of the jacobian of the element at each quadrature point */

    std::vector<Eigen::VectorXd> shape_functions_;    /*!< The shape function values of the element at each quadrature point */

    std::vector<Eigen::MatrixXd> derivs_nat_;         /*!< The element's derivatives at each quadrature point in the natural (ξ, η) coordinates */

    std::vector<Eigen::MatrixXd> derivs_;             /*!< The element's derivatives at each quadrature point in the physical (x, y) coordinates */


public:

    /**
     * \brief The default constructor of the BaseFem class.
     */
    explicit BaseFem() : quadrature_(), jacobians_(), inv_jacobians_(),
            det_jacobians_(), shape_functions_(), derivs_nat_(), derivs_() {}


    /**
     * \brief The destructor of the BaseFem class.
     */
    virtual ~BaseFem() {}


    /**
     * \brief To be used by the inheriting finite element to set the quadrature rule of the element.
     * \param [in] quad_points_num The number of the quadrature points for the element.
     * \return [void]
     */
    virtual void SetQuadrature(std::initializer_list<short> quad_points_num) = 0;


    /**
     * \brief To be used by the inheriting finite element to compute the Jacobian of the element for each quadrature point.
     * Internally should additionally compute and store the determinant and the inverse of the Jacobian.
     * \param [in] nodes The nodes composing the element.
     */
    virtual void ComputeJacobians(const std::vector<IMP::Vec<DIM, double> > &nodes) = 0;


    /**
     * \brief To be used by the inheriting finite element to compute the shape function values of the element for each quadrature point.
     * \return [void]
     */
    virtual void ComputeShapeFunctions() = 0;


    /**
     * \brief To be used by the inheriting finite element to compute the derivatives of the element's shape function in the natural coordinates (ξ, η) for each quadrature point. 
     * \return [void]
     */
    virtual void ComputeDerivsNatural() = 0;


    /**
     * \brief To be used by the inheriting finite element to compute the derivatives of the element's shape function in the physical coordinates (x, y) for each quadrature point. 
     * \return [void]
     */
    virtual void ComputeDerivs() = 0;


    /**
     * \brief To be used by the inheriting finite element to get the type of the element.
     * \return [FemType::UNKNOWN] The type of the element.
     */
    virtual FemType Type() const = 0;


    /**
     * \brief To be used by the inheriting finite element to get the quadrature points and weights of the element. 
     * \return [const QuadratureRule<DIM>&] The quadrature points and weights of the element.
     */
    inline auto & Quadrature() const { return this->quadrature_; }
    inline auto & Quadrature() { return this->quadrature_; }


    /**
     * \brief To be used by the inheriting finite element to get the Jacobian matrix of the element for each quadrature point.
     * \return [const std::vector<Eigen::MatrixXd>&] The Jacobian matrix of the element for each quadrature point.
     */
    inline auto & Jacobians() const { return this->jacobians_; }
    inline auto & Jacobians() { return this->jacobians_; }


    inline auto & Jacobians(std::size_t id) const { return this->jacobians_[id]; }
    inline auto & Jacobians(std::size_t id) { return this->jacobians_[id]; }


    /**
     * \brief To be used by the inheriting finite element to get the inverse of the Jacobian matrix of the element for each quadrature point.
     * \return [const std::vector<Eigen::MatrixXd>&] The inverse of the Jacobian matrix of the element for each quadrature point.
     */
    inline auto & InvJacobians() const { return this->inv_jacobians_; }
    inline auto & InvJacobians() { return this->inv_jacobians_; }


    inline auto & InvJacobians(std::size_t id) const { return this->inv_jacobians_[id]; }
    inline auto & InvJacobians(std::size_t id) { return this->inv_jacobians_[id]; }


    /**
     * \brief To be used by the inheriting finite element to get the determinant of the Jacobian matrix of the element for each quadrature point.
     * \return [const std::vector<double>&] The determinant of the Jacobian matrix of the element for each quadrature point.
     */
    inline auto & DetJacobians() const { return this->det_jacobians_; }
    inline auto & DetJacobians() { return this->det_jacobians_; }


    inline auto & DetJacobians(std::size_t id) const { return this->det_jacobians_[id]; }
    inline auto & DetJacobians(std::size_t id) { return this->det_jacobians_[id]; }


    /**
     * \brief To be used by the inheriting finite element to get the shape function values of the element for each quadrature point. 
     * \return [const std::vector<Eigen::VectorXd>&] The shape function values of the element for each quadrature point.
     */
    inline auto & ShapeFunctions() const { return this->shape_functions_; }
    inline auto & ShapeFunctions() { return this->shape_functions_; }


    inline auto & ShapeFunctions(std::size_t id) const { return this->shape_functions_[id]; }
    inline auto & ShapeFunctions(std::size_t id) { return this->shape_functions_[id]; }


    /**
     * \brief To be used by the inheriting finite element to get the derivatives of the element's shape function in the natural coordinates (ξ, η) for each quadrature point. 
     * \return [const std::vector<Eigen::MatrixXd>&] The derivatives of the element's shape function in the physical coordinates (ξ, η) for each quadrature point.
     */
    inline auto & DerivsNatural() const { return this->derivs_nat_; }
    inline auto & DerivsNatural() { return this->derivs_nat_; }


    inline auto & DerivsNatural(std::size_t id) const { return this->derivs_nat_[id]; }
    inline auto & DerivsNatural(std::size_t id) { return this->derivs_nat_[id]; }


    /**
     * \brief To be used by the inheriting finite element to get the derivatives of the element's shape function in the physical coordinates (x, y) for each quadrature point. 
     * \return [const std::vector<Eigen::MatrixXd>&] The derivatives of the element's shape function in the physical coordinates (x, y) for each quadrature point.
     */
    inline auto & Derivs() const { return this->derivs_; }
    inline auto & Derivs() { return this->derivs_; }


    inline auto & Derivs(std::size_t id) const { return this->derivs_[id]; }
    inline auto & Derivs(std::size_t id) { return this->derivs_[id]; }

};

/** \} End of Doxygen Groups */
/** \} End of Doxygen Groups */

} // End of namespace CLOUDEA

#endif // CLOUDEA_APPROXIMANTS_BASE_FEM_HPP_