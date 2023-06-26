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
   \file fem_mats.hpp
   \brief FemMats class header file.
   \author Konstantinos A. Mountris
   \date 28/02/2022
*/

#ifndef CLOUDEA_APPROXIMANTS_FEM_MATS_HPP_
#define CLOUDEA_APPROXIMANTS_FEM_MATS_HPP_


#include "CLOUDEA/engine/approximants/approximant_props.hpp"
#include "CLOUDEA/engine/approximants/fem_factory.hpp"

#include <IMP/IMP>
#include <Eigen/Dense>


namespace CLOUDEA {

/** \addtogroup Approximants \{ */
/** \addtogroup Finite-Elements \{ */

/**
 * \class FemMats
 * \brief Class implemmenting a collection of Fem matrices for the nodes of a mesh.
 * \tparam DIM The dimensions of the physical space that the element lies in.
 */
template<short DIM, short CELL_VERTS>
class FemMats
{
private:

    std::vector<std::vector<int>> fbar_patches_;

    std::vector<Eigen::VectorXd> shape_funcs_;

    std::vector<Eigen::MatrixXd> grads_;

    std::vector<double> jacobian_dets_;

    std::vector<double> gauss_weights_;

    int gauss_per_cell_;


public:

    /**
     * \brief The default constructor of the FemMats class.
     */
    explicit FemMats();


    /**
     * \brief The destructor of the FemMats class.
     */
    virtual ~FemMats();


    /**
     * \brief Compute non-overlapping patches where each elemnent of the mesh belong to be used in F-bar simplex elements for
     * volumetric locking correction at near incompressibility.
     * \return [void]
     */
    inline void ComputeFbarPatches(IMP::Mesh<DIM,CELL_VERTS> mesh);


    /**
     * \brief To be used by the inheriting finite element to set the quadrature rule of the element.
     * \param [in] quad_points_num The number of the quadrature points for the element.
     * \return [void]
     */
    inline void ComputeMatrices(const IMP::Mesh<DIM,CELL_VERTS> &mesh);


    inline auto & FbarPatches() const { return this->fbar_patches_; }


    inline auto & FbarPatches(std::size_t id) const { return this->fbar_patches_[id]; }


    /**
     * \brief Get the finite element shape functions.
     * \return [std::vector<Eigen::VectorXd>&] The finite element shape functions.
     */
    inline auto & ShapeFuncs() const { return this->shape_funcs_; }


    /**
     * \brief Get the finite element shape functions for gauss point with index id.
     * Fast access, no range check.
     * \return [Eigen::VectorXd&] The finite element shape functions for gauss point with index id.
     */
    inline auto & ShapeFuncs(std::size_t id) const { return this->shape_funcs_[id]; }


    /**
     * \brief Get the finite element shape functions for gauss point with index id.
     * Slow acces, with range check.
     * \return [Eigen::VectorXd&] The finite element shape functions for gauss point with index id.
     */
    inline auto & ShapeFuncsAt(std::size_t id) const { return this->shape_funcs_.at(id); }


    /**
     * \brief Get the finite element gradients.
     * \return [std::vector<Eigen::MatrixXd>&] The finite element gradients.
     */
    inline auto & Grads() const { return this->grads_; }


    /**
     * \brief Get the finite element gradients for gauss point with index id.
     * Fast access, no range check.
     * \return [Eigen::VectorXd&] The finite element gradients for gauss point with index id.
     */
    inline auto & Grads(std::size_t id) const { return this->grads_[id]; }


    /**
     * \brief Get the finite element gradients for gauss point with index id.
     * Slow acces, with range check.
     * \return [Eigen::VectorXd&] The finite element gradients for gauss point with index id.
     */
    inline auto & GradsAt(std::size_t id) const { return this->grads_.at(id); }


    /**
     * \brief Get the finite element jacobian determinants.
     * \return [std::vector<double>&] The finite element jacobian determinants.
     */
    inline auto & JacobianDets() const { return this->jacobian_dets_; }


    /**
     * \brief Get the finite element jacobian determinants for gauss point with index id.
     * Fast access, no range check.
     * \return [Eigen::VectorXd&] The finite element jacobian determinants for gauss point with index id.
     */
    inline auto & JacobianDets(std::size_t id) const { return this->jacobian_dets_[id]; }


    /**
     * \brief Get the finite element jacobian determinants for gauss point with index id.
     * Slow acces, with range check.
     * \return [Eigen::VectorXd&] The finite element jacobian determinants for gauss point with index id.
     */
    inline auto & JacobianDetsAt(std::size_t id) const { return this->jacobian_dets_.at(id); }


    /**
     * \brief Get the finite element gaussian weights.
     * \return [std::vector<double>&] The finite element gaussian weights.
     */
    inline auto & GaussWeights() const { return this->gauss_weights_; }


    /**
     * \brief Get the finite element gaussian weight for gauss point with index id.
     * Fast access, no range check.
     * \return [Eigen::VectorXd&] The finite element gaussian weight for gauss point with index id.
     */
    inline auto & GaussWeights(std::size_t id) const { return this->gauss_weights_[id]; }


    /**
     * \brief Get the finite element gaussian weight for gauss point with index id.
     * Slow acces, with range check.
     * \return [Eigen::VectorXd&] The finite element gaussian weight for gauss point with index id.
     */
    inline auto & GaussWeightsAt(std::size_t id) const { return this->gauss_weights_.at(id); }


    /**
     * \brief Get the number of gauss points per cell in the mesh.
     * \return [int] The number of gauss points per cell in the mesh.
     */
    inline auto GaussPerCell() const { return this->gauss_per_cell_; }


};

/** \} End of Doxygen Groups */
/** \} End of Doxygen Groups */

} // End of namespace CLOUDEA

#endif // CLOUDEA_APPROXIMANTS_FEM_MATS_HPP_

#include "CLOUDEA/engine/approximants/fem_mats.tpp"