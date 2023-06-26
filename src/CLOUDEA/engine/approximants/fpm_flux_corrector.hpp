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
   \file fpm_flux_corrector.hpp
   \brief FpmFluxCorrector class header file.
   \author Konstantinos A. Mountris
   \date 18/01/2022
*/

#ifndef CLOUDEA_APPROXIMANTS_FPM_FLUX_CORRECTOR_HPP_
#define CLOUDEA_APPROXIMANTS_FPM_FLUX_CORRECTOR_HPP_

#include <IMP/IMP>
#include <Eigen/Dense>

#include <vector>
#include <iostream>
#include <stdexcept>
#include <exception>
#include <cmath>


namespace CLOUDEA {

/** \addtogroup Approximants \{ */
/**  \addtogroup Fragile-Points \{ */


/**
 * \class FpmFluxCorrector
 * \brief Class storing data for efficient numerical flux correction application on the internal facets
 * of a Voronoi tesselation during simulations using the Fragile Points Method approximation.
 * \tparam DIM the dimensions of the domain that the facets belong.
 * \author Konstantinos A. Mountris
 */
template<short DIM>
class FpmFluxCorrector
{
private:

    std::vector<int> facet_ids_;                            /**< The indices of the internal facets */

    std::vector<double> facet_measures_;                    /**< The measures of the internal facests. 1D: length, 2D: area */

    std::vector<Eigen::VectorXd> facet_centroids_;          /**< The centroids of the internal facets */

    std::vector<Eigen::VectorXd> facet_dir_;                /**< The direction vectors of the internal facets */

    std::vector<Eigen::VectorXd> facet_normals_;            /**< The normal unit vectors of the internal facets */

public:


    /**
     * \brief Default constructor for FpmFluxCorrector object.
     */
    FpmFluxCorrector();


    /**
     * \brief Destructor for FpmFluxCorrector object.
     */
    virtual ~FpmFluxCorrector();


    void ComputeCorrectorData(const IMP::Voronoi<DIM> &voro);


    inline auto FacetsNum() const { return static_cast<int>(this->facet_ids_.size()); }


    inline auto & FacetIds() const { return this->facet_ids_; }


    inline auto & FacetIds(std::size_t id) const { return this->facet_ids_[id]; }


    inline auto & FacetMeasures() const { return this->facet_measures_; }


    inline auto & FacetMeasures(std::size_t id) const { return this->facet_measures_[id]; }


    inline auto & FacetCentroids() const { return this->facet_centroids_; }


    inline auto & FacetCentroids(std::size_t id) const { return this->facet_centroids_[id]; }


    inline auto & FacetDir() const { return this->facet_dir_; }


    inline auto & FacetDir(std::size_t id) const { return this->facet_dir_[id]; }


    inline auto & FacetNormals() const { return this->facet_normals_; }


    inline auto & FacetNormals(std::size_t id) const { return this->facet_normals_[id]; }


protected:

    void ComputeCorrectorData2D(const IMP::Voronoi<DIM> &voro);


    void ComputeCorrectorData3D(const IMP::Voronoi<DIM> &voro);

};


/** \} End of Doxygen Groups */
/** \} End of Doxygen Groups */

} // End of namespace CLOUDEA

#endif // CLOUDEA_APPROXIMANTS_FPM_FLUX_CORRECTOR_HPP_

#include "CLOUDEA/engine/approximants/fpm_flux_corrector.tpp"