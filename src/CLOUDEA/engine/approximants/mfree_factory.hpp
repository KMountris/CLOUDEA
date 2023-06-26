
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
   \file mfree_factory.hpp
   \brief MfreeFactory class header file.
   \author Konstantinos A. Mountris
   \date 28/11/2019
*/

#ifndef CLOUDEA_APPROXIMANTS_MFREE_FACTORY_HPP_
#define CLOUDEA_APPROXIMANTS_MFREE_FACTORY_HPP_


#include "CLOUDEA/engine/approximants/mfree.hpp"
#include "CLOUDEA/engine/approximants/mls.hpp"
#include "CLOUDEA/engine/approximants/rpi.hpp"
#include "CLOUDEA/engine/approximants/mki.hpp"
#include "CLOUDEA/engine/utilities/attributes.hpp"
#include "CLOUDEA/engine/utilities/logger.hpp"

#include <memory>
#include <string>


namespace CLOUDEA {

/** \addtogroup Approximants \{ */
/** \addtogroup Meshfree \{ */

/**
 * \class MfreeFactory
 * \author Konstantinos A. Mountris
 * \brief Factory for meshfree (Mfree) approximants generation.
 * \tparam DIM The spatial dimensions of the generated approximant.
 */

template <short DIM>
class MfreeFactory
{
public:

    /**
     * \brief Creates a meshfree approximant according to the provided type.
     * \param [in] mfree_type The type of the desired meshfree approximant.
     * \return [std::unique_ptr<Mfree>] Unique pointer to the created meshfree approximant.
     */
    inline static std::unique_ptr<Mfree<DIM>> Create(MfreeType mfree_type) {

        std::unique_ptr<Mfree<DIM>> mfree_ptr;

        switch (mfree_type) {
        case MfreeType::mls :
            mfree_ptr = std::make_unique<Mls<DIM>>();
            break;
        case MfreeType::rpi :
            mfree_ptr = std::make_unique<Rpi<DIM>>();
            break;
        case MfreeType::mki :
            mfree_ptr = std::make_unique<Mki<DIM>>();
            break;
        default:
            std::string error_str = "Could not create Mfree approximant. Not supported type.";
            throw std::invalid_argument(Logger::Error(error_str, __FILE__, __LINE__));
            break;
        }
        return mfree_ptr;
    }
};


/** \} End of Doxygen Groups */
/** \} End of Doxygen Groups */

} // End of namespace ELECTRA

#endif // CLOUDEA_APPROXIMANTS_MFREE_FACTORY_HPP_