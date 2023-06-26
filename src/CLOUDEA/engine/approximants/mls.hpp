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
   \file mls.hpp
   \brief Mls class header file.
   \author Konstantinos A. Mountris
   \date 22/05/2019
*/

#ifndef CLOUDEA_APPROXIMANTS_MLS_HPP_
#define CLOUDEA_APPROXIMANTS_MLS_HPP_

#include "CLOUDEA/engine/approximants/mfree.hpp"
#include "CLOUDEA/engine/utilities/logger.hpp"

#include <vector>
#include <stdexcept>
#include <exception>
#include <iterator>
#include <set>


namespace CLOUDEA {

/**  \addtogroup Approximants \{ */
/** \addtogroup Meshfree \{ */


/**
 * \class Mls
 * \author Konstantinos A. Mountris
 * \brief The Moving Least Squares (MLS) meshfree approximant.
 * \tparam DIM The spatial dimensions of the approximant.
 */
template <short DIM>
class Mls : public Mfree<DIM>
{

public:


   /**
    * \brief The constructor of the moving least squares (MLS) approximants.
    */
   inline Mls();


   /**
    * \brief The destructor of the moving least squares (MLS) approximants.
    */
   inline ~Mls();


   /**
    * \brief Compute the MLS approximant for a group of evaluation points. Both basis functions and their gradients are computed.
    * \param [in] eval_points The points where the approximant is evaluated.
    * \param [in] field_nodes The field nodes of the domain where the approximant is applied.
    * \return [void]
    */
   virtual void Compute(const std::vector<IMP::Vec<DIM, double>> &eval_points,
                        const std::vector<IMP::Vec<DIM, double>> &field_nodes);


   /**
    * \brief Compute the MLS approximant for a single evaluation point. Both basis functions and their gradients are computed.
    * \param [in] eval_point The single point where the approximant is evaluated.
    * \param [in] eval_point_id The index of the corresponding support domain to the evaluation point.
    * \param [in] field_nodes The field nodes of the domain where the approximant is applied.
    * \return [void]
    */
   virtual void Compute(const IMP::Vec<DIM, double> &eval_point, int eval_point_id,
                        const std::vector<IMP::Vec<DIM, double>> &field_nodes);

};

/** \} End of Doxygen Groups */
/** \} End of Doxygen Groups */

} //End of namespace CLOUDEA

#endif //CLOUDEA_APPROXIMANTS_MLS_HPP_

#include "CLOUDEA/engine/approximants/mls.tpp"