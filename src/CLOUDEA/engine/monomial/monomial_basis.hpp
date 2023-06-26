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
   \file monomial_basis.hpp
   \brief MonomialBasis class header file.
   \author Konstantinos A. Mountris
   \date 29/07/2019
*/

#ifndef CLOUDEA_MONOMIAL_MONOMIAL_BASIS_HPP_
#define CLOUDEA_MONOMIAL_MONOMIAL_BASIS_HPP_

#include <IMP/IMP>
#include <armadillo>

#include <vector>

namespace CLOUDEA {

/** \addtogroup Monomial \{ */

/**
 * \enum MonomialType
 * \author Konstantinos A. Mountris
 * \brief Types of monomials.
 */
enum class MonomialType { unknown,      /**< Unknown monomial type */
                          linear,       /**< Linear monomial type */
                          quadratic,    /**< Quadratic monomial type */
                          cubic         /**< Cubic monomial type */
                        };

/**
 * \class MonomialBasis
 * \author Konstantinos A. Mountris
 * \brief Class implemmenting a monomial basis.
 * \tparam DIM The spatial dimension of the monomial basis.
 */
template <short DIM>
class MonomialBasis
{

protected:

  arma::mat vandermonde_;             /**< The Vandermonde matrix of the monomial basis */

  std::vector<arma::mat> grad_;       /**< The gradient matrix for each row entry of the Vandermonde matrix */  

  std::vector<arma::rowvec> basis_;   /**< The values for each row entry of the Vandermonde matrix */

  MonomialType type_;                       /**< The type of the monomial basis */

  short rank_;                                /**< The rank of the monomial basis */



public:

  /**
   * \brief The MonomialBasis default constructor.
   */
  inline MonomialBasis() :  vandermonde_(), grad_(), basis_(), type_(MonomialType::unknown), rank_(0) {}


  /**
   * \brief The MonomialBasis destructor.
   */
  inline virtual ~MonomialBasis() {}
    
    
  /**
   * \brief Compute monomial basis for a given point. 
   * \param [in] point The coordinates of the point.
   * \return [void]
   */
  inline virtual void Compute(const IMP::Vec<DIM, double> &point) = 0;
  
  
  /**
   * \brief Compute monomial basis for a vector of points. 
   * \param point The coordinates of the point.
   * \return [void]
   */
  inline virtual void Compute(const std::vector<IMP::Vec<DIM, double>> &points) = 0;


  /**
   * \brief Get the Vandermonde matrix of the monomial basis.
   * \return [const arma::mat&] The Vandermonde matrix of the monomial basis.
   */
  inline const arma::mat & Vandermonde() const { return this->vandermonde_; }
  
  
  /**
   * \brief Get the gradient matrix of the monomial basis for a row entry in the Vandemonde matrix. 
   * \param [in] id The row index of the Vandermonde matrix pointing to the basis for which we get the gradient. 
   *                Default: Get the gradient for the basis at the first row of the Vandermonde matrix (id=0).
   * \return [const arma::mat&] The gradient matrix of the monomial basis for a row entry in the Vandemonde matrix.
   */
  inline const arma::mat & Grad(std::size_t id) const { return this->grad_[id]; }
  
  
  /**
   * \brief Get the monomial basis from the Vandermonde matrix.
   * \param [in] id The row index of the Vandermonde matrix. Default: get the basis at the first row of the Vandermonde matrix (id=0). 
   * \return [const arma::rowvec&] The monomial basis from the Vandermonde matrix.
   */
  inline const arma::rowvec & Basis(std::size_t id) const { return this->basis_[id]; }


  /**
   * \brief Get the type of the monomial basis.  
   * \return [MonomialType] The type of the monomial basis. 
   */
  inline MonomialType Type() const { return this->type_; }


  /**
   * \brief Return the rank of the monomial basis.
   * The rank denotes the number of monomials in the basis. 
   * \return [short] The rank of the monomial basis.
   */
  inline short Rank() const { return this->rank_; }


};


} //End of namespace CLOUDEA

#endif //CLOUDEA_MONOMIAL_MONOMIAL_BASIS_HPP_