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
   \file mfree.hpp
   \brief Mfree class header file.
   \author Konstantinos A. Mountris
   \date 26/11/2019
*/

#ifndef CLOUDEA_APPROXIMANTS_MFREE_HPP_
#define CLOUDEA_APPROXIMANTS_MFREE_HPP_

#include "CLOUDEA/engine/weight_functions/weight_factory.hpp"
#include "CLOUDEA/engine/monomial/monomial_factory.hpp"
#include "CLOUDEA/engine/support_domain/support_domain.hpp"

#include <IMP/IMP>
#include <Eigen/Dense>
#include <armadillo>

#include <cstddef>
#include <vector>
#include <memory>


namespace CLOUDEA {

/**  \addtogroup Approximants \{ */
/** \addtogroup Meshfree \{ */


/**
 * \class Mfree
 * \author Konstantinos A. Mountris
 * \brief A basic meshfree approximant.
 * \tparam DIM the spatial dimensions of the approximant.
 */
template<short DIM>
class Mfree
{

protected:

    SupportDomain<DIM> support_;                        /**< The support domain of the approximant */

    std::vector<Eigen::MatrixXd> phi_grad_;             /**< The basis function gradient values of the approximant on the evaluation points */

    std::vector<Eigen::VectorXd> phi_;                  /**< The basis function values of the approximant on the evaluation points */

    std::unique_ptr<WeightFunction<DIM>> weight_;       /**< The weight function of the approximant. */

    std::unique_ptr<MonomialBasis<DIM>> mbasis_;        /**< The monomial basis of the approximant. */

    MfreeType type_;                                    /**< The type of the approximant */


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


public:
    /**
     * \brief Mfree constructor.
     */
    inline Mfree() : support_(), phi_grad_(), phi_(), weight_(), mbasis_(), type_(MfreeType::unknown) {}


    /**
     * \brief Mfree destructor.
     */
    inline virtual ~Mfree() {}


    /**
    * \brief Set the weight function of the approximant.
    * \param type The type of the weight function.
    * \param theta The shape parameter of the weight function (Applicable at specific types of weight function).
    * \param beta The exponent of the weight function (Applicable at specific types of weight function).
    */
    inline void SetWeightFunction(WeightType type, double theta = 1.0, double beta = 1.0) { this->weight_ = WeightFactory<DIM>::Create(type);
                                                                                            this->weight_->SetParameters(theta, beta);
                                                                                          }


    /**
    * \brief Set the monomial basis for the approximant.
    * \param type The type of the monomial basis.
    */
    inline void SetMonomialBasis(MonomialType type) { this->mbasis_ = MonomialFactory<DIM>::Create(type); }


    /**
    * \brief Set the dilatation coefficient of the support domains. It can vary amongst different support domains.
    * \param [in] dilate_coeffs The dilatation coefficient values per support domain.
    * \return [void]
    */
    inline void SetSupportDilatation(std::vector<double> dilate_coeffs) { this->support_.SetDilateCoeff(dilate_coeffs); }


    /**
    * \brief Set the dilatation coefficient of the support domains. Same for all support domains.
    * \param [in] dilate_coeff The dilatation coefficient value.
    * \return [void]
    */
    inline void SetSupportDilatation(double dilate_coeff) { this->support_.SetDilateCoeff(dilate_coeff); }


    /**
     * \brief Compute the approximant for a group of evaluation points. Both basis functions and their gradients are computed.
     * \param [in] eval_points The points where the approximant is evaluated.
     * \param [in] field_nodes The field nodes of the domain where the approximant is applied.
     * \return [void]
     */
    virtual void Compute(const std::vector<IMP::Vec<DIM, double>> &eval_points,
                         const std::vector<IMP::Vec<DIM, double>> &field_nodes) = 0;


    /**
     * \brief Compute the approximant for a single evaluation point. Both basis functions and their gradients are computed.
     * \param [in] eval_point The single point where the approximant is evaluated.
     * \param [in] eval_point_id The index of the corresponding support domain to the evaluation point.
     * \param [in] field_nodes The field nodes of the domain where the approximant is applied.
     * \return [void]
     */
    virtual void Compute(const IMP::Vec<DIM, double> &eval_point, int eval_point_id,
                         const std::vector<IMP::Vec<DIM, double>> &field_nodes) = 0;


    /**
     * \brief Get the type of the meshfree approximant.
     * \return [CLOUDEA::MfreeType] The type of the meshfree approximant.
     */
    inline MfreeType Type() const { return this->type_; }


    /**
     * \brief Edit the support domain for each field node of the approximant.
     * \return [SupportDomain<DIM>&] The support domain for each field node of the approximant.
     */
    inline SupportDomain<DIM> & EditSupport() { return this->support_; }


    /**
     * \brief Get the support domain for each field node of the approximant.
     * \return [const SupportDomain<DIM>&] The support domain for each field node of the approximant.
     */
    inline const SupportDomain<DIM> & Support() const { return this->support_; }


    /**
     * \brief Get the basis functions for each evaluation point of the approximant.
     * \return [const std::vector<Eigen::VectorXd>&] The basis functions for each evaluation point of the approximant.
     */
    inline const std::vector<Eigen::VectorXd> & Phi() const { return this->phi_; }


    /**
     * \brief Get the basis function for a specified evaluation point of the approximant.
     * Fast access, no range check.
     * \param [in] id The index of the specified evaluation point of the approximant.
     * \return [const std::vector<Eigen::VectorXd>&] The basis functions for each evaluation point of the approximant.
     */
    inline const Eigen::VectorXd & Phi(std::size_t id) const { return this->phi_[id]; }


    /**
     * \brief Get the basis function for a specified evaluation point of the approximant.
     * Slower access with range check.
     * \param [in] id The index of the specified evaluation point of the approximant.
     * \return [const std::vector<Eigen::VectorXd>&] The basis functions for each evaluation point of the approximant.
     */
    inline const Eigen::VectorXd & PhiAt(std::size_t id) const { return this->phi_.at(id); }


    /**
     * \brief Get the basis functions gradients for each evaluation point of the approximant.
     * \return [const std::vector<Eigen::MatrixXd>&] The basis functions gradients for each evaluation point of the approximant.
     */
    inline const std::vector<Eigen::MatrixXd> & PhiGrad() const { return this->phi_grad_; }


    /**
     * \brief Get the basis functions gradients for a specified evaluation point of the approximant.
     * Fast access, no range check.
     * \param [in] id The index of the specified evaluation point of the approximant.
     * \return [const std::vector<Eigen::MatrixXd>&] The basis functions gradients for each evaluation point of the approximant.
     */
    inline const Eigen::MatrixXd & PhiGrad(std::size_t id) const { return this->phi_grad_[id]; }


    /**
     * \brief Get the basis functions gradients for a specified evaluation point of the approximant.
     * Slower access with range check.
     * \param [in] id The index of the specified evaluation point of the approximant.
     * \return [const std::vector<Eigen::MatrixXd>&] The basis functions gradients for each evaluation point of the approximant.
     */
    inline const Eigen::MatrixXd & PhiGradAt(std::size_t id) const { return this->phi_grad_.at(id); }

};


/** \} End of Doxygen Groups */
/** \} End of Doxygen Groups */

} //end of namespace CLOUDEA

#endif //CLOUDEA_APPROXIMANTS_MFREE_HPP_
