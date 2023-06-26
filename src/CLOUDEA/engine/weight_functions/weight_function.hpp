 /*
 * CLOUDEA - Software for solving PDEs using explicit methods.
 * Copyright (C) 2017  <Konstantinos A. Mountris> <konstantinos.mountris\gmail.com>
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
   \file weight_function.hpp
   \brief WeightFunction class header file.
   \author Konstantinos A. Mountris
   \date 10/07/2019
*/

#ifndef CLOUDEA_WEIGHT_FUNCTIONS_WEIGHT_FUNCTION_HPP_
#define CLOUDEA_WEIGHT_FUNCTIONS_WEIGHT_FUNCTION_HPP_


#include <IMP/IMP>


namespace CLOUDEA {

/**  \addtogroup WeightFunctions \{ */


/**
 * \enum Various types of weight functions.
 * \author Konstantinos A. Mountris
 */
enum class WeightType { unknown,            /**< Uknown weight function type */
                        cubic,              /**< Cubic spline weight function type */
                        quartic,            /**< Quartic spline weight function type */
                        gaussian,           /**< Gaussian radial basis weight function type */
                        multiquadric,       /**< Multiquadric radial basis weight function type */
                        polyharmonic        /**< Polyharmonic weight function type */
                      };


/**
 * \class WeightFunction
 * \brief Class implemmenting various weight interpolation functions.
 * \tparam DIM The spatial dimension of the weight function.
 */
template <short DIM>
class WeightFunction
{

protected:

    IMP::Vec<DIM, double> grad_;        /**< Gradient of the weight function at a given evaluation point */

    double val_;                        /**< Value of the weight function at a given evaluation point */

    double theta_;                      /**< The weight function shape parameter */

    double beta_;                       /**< The weight function exponent */

    WeightType type_;                   /**< Type of the weight function */


public:

    /**
     * \brief The weight function default constructor. 
     */
    inline WeightFunction() : grad_(), val_(0.), theta_(1.), beta_(1.), type_(WeightType::unknown) {}


    /**
     * \brief The weight function destructor. 
     */
    inline virtual ~WeightFunction() {}


    /**
     * \brief Set the weight functions shape parameter and exponent.
     * \param [in] theta The shape parameter of the weight function. 
     * \param [in] beta The exponent of the weight function.
     * \return [void]
     */
    inline void SetParameters(double theta, double beta) { this->theta_ = theta; this->beta_ = beta; }


    /**
     * \brief Compute the weight function. 
     * \param [in] point The point to evaluate the weight function.
     * \param [in] center The center of the weight function. 
     * \param [in] radius The radius of the weight function.
     * \param [in] dilate_coeff The dilatation coefficient for the weight function's radius.
     * \return [void]
     */
    virtual void Compute(const IMP::Vec<DIM, double> &point, const IMP::Vec<DIM, double> &center, double radius, double dilate_coeff) = 0;


    /**
     * \brief Get the value of the weight function.
     * \return [double] The value of the weight function.
     */
    inline double Val() const { return this->val_; }


    /**
     * \brief Get the gradient of the weight function.
     * \return [const IMP::Vec<DIM, double>&] The gradient of the weight function. 
     */
    inline const IMP::Vec<DIM, double> & Grad() const { return this->grad_; } 


    /**
     * \brief Get the value of the gradient of the weight function for a dimension.
     * \param id The index of the required dimension. [X:0, Y:1, Z:2]
     * \return [double] The value of the gradient of the weight function at the required dimension. 
     */
    inline double Grad(std::size_t id) const { return this->grad_[id]; } 


};

/** \} End of Doxygen Groups */

} // End of namespace CLOUDEA

#endif //CLOUDEA_WEIGHT_FUNCTIONS_WEIGHT_FUNCTION_HPP_