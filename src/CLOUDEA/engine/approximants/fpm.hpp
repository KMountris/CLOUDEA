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
   \file fpm.hpp
   \brief Fpm class header file.
   \author Konstantinos A. Mountris
   \date 12/10/2020
*/

#ifndef CLOUDEA_APPROXIMANTS_FPM_HPP_
#define CLOUDEA_APPROXIMANTS_FPM_HPP_

#include "CLOUDEA/engine/support_domain/support_domain.hpp"
#include "CLOUDEA/engine/approximants/fpm_flux_corrector.hpp"

#include <IMP/IMP>
#include <Eigen/Dense>
#include <armadillo>

#include <cstddef>
#include <vector>
#include <memory>


namespace CLOUDEA {

/**  \addtogroup Approximants \{ */
/**  \addtogroup Fragile-Points \{ */

/**
 * \class Fpm
 * \author Konstantinos A. Mountris
 * \brief A class implementing an approximant using the Fragile Points Method.
 * \tparam DIM the spatial dimensions of the approximant.
 */
template<short DIM>
class Fpm
{

private:

    SupportDomain<DIM> support_;                  /**< The support domain of the approximant */

    FpmFluxCorrector<DIM> flux_corrector_;        /**< The numerical flux corrector of the fpm approximant */

    std::vector<Eigen::MatrixXd> phi_grad_;       /**< The basis function gradient values of the approximant on the evaluation points */

    std::vector<Eigen::VectorXd> phi_;            /**< The basis function values of the approximant on the evaluation points */

    double penalty_;                              /**< The penalty coefficienc of the fpm approximant */

public:
    /**
     * \brief Fpm constructor.
     */
    inline Fpm();


    /**
     * \brief Fpm destructor.
     */
    inline virtual ~Fpm();


    /**
     * \brief Set up the numerical flux corrector of the fpm approximant.
     * \param voro The voronoi tesselation to set up the data of the numerical flux corrector.
     * \return [void]
     */
    inline void SetFluxCorrector(const IMP::Voronoi<DIM> &voro);


    /**
     * \brief Set the penalty coefficient of the fpm approximant.
     * \param [in] penalty The penalty coefficient of the fpm approximant.
     * \return [void]
     */
    inline void SetPenalty(double penalty) { this->penalty_ = penalty; }


    /**
     * \brief Compute the fpm approximant for the field nodes of a voronoi tesselation.
     * The voronoi tesselation is composed by field nodes, voronoi points, voronoi cells, and optionally voronoi facets for 3D tesselations.
     * \param [in] voro The voronoi tesselation of which field nodes will be used for the fpm approximant evaluation.
     * \return [void]
     */
    void Compute(const IMP::Voronoi<DIM> &voro);


    /**
     * \brief Edit the support domain for each field node of the approximant.
     * \return [SupportDomain<DIM>&] The support domain for each field node of the approximant.
     */
    inline auto & Support() { return this->support_; }


    /**
     * \brief Get the support domain for each field node of the approximant.
     * \return [const SupportDomain<DIM>&] The support domain for each field node of the approximant.
     */
    inline auto & Support() const { return this->support_; }


    /**
     * \brief Get the numerical flux corrector of the approximant.
     * \return [FluxCorrector<DIM> &] The numerical flux corrector of the approximant.
     */
    inline auto & FluxCorrector() const { return this->flux_corrector_; }


    /**
     * \brief Get the basis functions for each evaluation point of the approximant.
     * \return [const std::vector<Eigen::VectorXd>&] The basis functions for each evaluation point of the approximant.
     */
    inline auto & Phi() const { return this->phi_; }


    /**
     * \brief Get the basis function for a specified evaluation point of the approximant.
     * Fast access, no range check.
     * \param [in] id The index of the specified evaluation point of the approximant.
     * \return [const std::vector<Eigen::VectorXd>&] The basis functions for each evaluation point of the approximant.
     */
    inline auto & Phi(std::size_t id) const { return this->phi_[id]; }


    /**
     * \brief Get the basis function for a specified evaluation point of the approximant.
     * Slower access with range check.
     * \param [in] id The index of the specified evaluation point of the approximant.
     * \return [const std::vector<Eigen::VectorXd>&] The basis functions for each evaluation point of the approximant.
     */
    inline auto & PhiAt(std::size_t id) const { return this->phi_.at(id); }


    /**
     * \brief Get the basis functions gradients for each evaluation point of the approximant.
     * \return [const std::vector<Eigen::MatrixXd>&] The basis functions gradients for each evaluation point of the approximant.
     */
    inline auto & PhiGrad() const { return this->phi_grad_; }


    /**
     * \brief Get the basis functions gradients for a specified evaluation point of the approximant.
     * Fast access, no range check.
     * \param [in] id The index of the specified evaluation point of the approximant.
     * \return [const std::vector<Eigen::MatrixXd>&] The basis functions gradients for each evaluation point of the approximant.
     */
    inline auto & PhiGrad(std::size_t id) const { return this->phi_grad_[id]; }


    /**
     * \brief Get the basis functions gradients for a specified evaluation point of the approximant.
     * Slower access with range check.
     * \param [in] id The index of the specified evaluation point of the approximant.
     * \return [const std::vector<Eigen::MatrixXd>&] The basis functions gradients for each evaluation point of the approximant.
     */
    inline auto & PhiGradAt(std::size_t id) const { return this->phi_grad_.at(id); }


    /**
     * \brief Get the penalty coefficient of the fpm approximant.
     * \return [double] The penalty coefficient of the fpm approximant.
     */
    inline auto Penalty() const { return this->penalty_; }



protected:

    /**
     * \brief Cast an armadillo matrix to Eigen.
     * \param arma_A The armadillo matrix to be casted.
     * \return Eigen::MatrixXd The casted matrix to Eigen.
     */
    inline Eigen::MatrixXd CastToEigen(arma::mat arma_A) {
      Eigen::MatrixXd eigen_B = Eigen::Map<Eigen::MatrixXd>(arma_A.memptr(), arma_A.n_rows, arma_A.n_cols);
      return eigen_B;
    }


    /**
     * \brief Cast an armadillo vector to Eigen.
     * \param arma_A The armadillo vector to be casted.
     * \return Eigen::VectorXd The casted vector to Eigen.
     */
    inline Eigen::VectorXd CastToEigen(arma::vec arma_A) {
      Eigen::VectorXd eigen_B = Eigen::Map<Eigen::VectorXd>(arma_A.memptr(), arma_A.n_rows, arma_A.n_cols);
      return eigen_B;
    }


    /**
     * \brief Cast an armadillo row vector to Eigen.
     * \param arma_A The armadillo row vector to be casted.
     * \return Eigen::RowVectorXd The casted row vector to Eigen.
     */
    inline Eigen::RowVectorXd CastToEigen(arma::rowvec arma_A) {
      Eigen::RowVectorXd eigen_B = Eigen::Map<Eigen::RowVectorXd>(arma_A.memptr(), arma_A.n_rows, arma_A.n_cols);
      return eigen_B;
    }


};


/** \} End of Doxygen Groups */
/** \} End of Doxygen Groups */

} //end of namespace CLOUDEA

#endif //CLOUDEA_APPROXIMANTS_FPM_HPP_

#include "CLOUDEA/engine/approximants/fpm.tpp"