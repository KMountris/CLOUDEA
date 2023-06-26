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


#ifndef CLOUDEA_INTEGRATION_INTEG_POINTS_HPP_
#define CLOUDEA_INTEGRATION_INTEG_POINTS_HPP_

/*!
   \file integ_points.hpp
   \brief IntegPoints class header file.
   \author Konstantinos A. Mountris
   \date 16/05/2017
*/


#include "CLOUDEA/engine/approximants/mmls_3d.hpp"
#include "CLOUDEA/engine/support_domain/inf_support_domain.hpp"
#include "CLOUDEA/engine/vectors/vec3.hpp"
#include "CLOUDEA/engine/elements/element_properties.hpp"
#include "CLOUDEA/engine/elements/node.hpp"
#include "CLOUDEA/engine/elements/tetrahedron.hpp"
#include "CLOUDEA/engine/mesh/mesh_properties.hpp"
#include "CLOUDEA/engine/mesh/tetramesh.hpp"
#include "CLOUDEA/engine/integration/integ_options.hpp"

#include <string>
#include <vector>
#include <utility>
#include <algorithm>
#include <iostream>


namespace CLOUDEA {

/*!
 *  \addtogroup Integration
 *  @{
 */


/*!
 * \class IntegPoints
 * \brief Class implemmenting integration points for numerical appoximation of integrals in weak-form models.
 */

class IntegPoints{
public:

    /*!
     * \brief IntegPoints constructor.
     */
    IntegPoints();


    /*!
     * \brief IntegPoints destructor.
     */
    virtual ~IntegPoints();


    /*!
     * \brief Generate the integration points for a given tetrahedral mesh.
     * \param [in] tetramesh The tetrahedral mesh for which integration points will be generated.
     * \param [in] options The integration options for the generation of integration points.
     * \param [in] support_dom The support (influence) domains of the tetrahedral mesh nodes. Used for adaptive integration.
     * \return [void]
     */
    void Generate(const TetraMesh &tetramesh, const IntegOptions &options, const InfSupportDomain &support_dom);

    void Generate2(const TetraMesh &tetramesh, short num_per_elem);

    void LoadFromFile(const std::string & ip_file);


    /*!
     * \brief Get the coordinates of the integration points.
     * \return [std::vector<CLOUDEA::Vec3<double> >] The coordinates of the integration points.
     */
    inline const std::vector<Vec3<double> > & Coordinates() const { return this->coordinates_; }


    /*!
     * \brief Get the weights of the integration points.
     *
     * The weights carry volumetric information associated with the integration points.
     *
     * \return [std::vector<double>] The weights of the integration points.
     */
    inline const std::vector<double> & Weights() const { return this->weights_; }


    /*!
     * \brief Get the number of integration points.
     * \return [int] The number of integration points.
     */
    inline int PointsNum() const {
        if (this->coordinates_.size() != this->weights_.size()) {
            std::string error = "[CLOUDEA ERROR] Cannot retrieve number of integration points. The lists of integration"
                                " points' coordinates and weights do not have the same size.";
            throw std::runtime_error(error.c_str());
        }

        return static_cast<int>(this->coordinates_.size());
    }



protected:


    /*!
     * \brief Generate one integration point per tetrahedron of a tetrahedral mesh.
     * \param [in] tetramesh The tetrahedral mesh.
     * \return [void]
     */
    void GenerateOnePointPerTetra(const TetraMesh &tetramesh);


    /*!
     * \brief Generate four integration points per tetrahedron of a tetrahedral mesh.
     * \param [in] tetramesh The tetrahedral mesh.
     * \return [void]
     */
    void GenerateFourPointsPerTetra(const TetraMesh &tetramesh);


    /*!
     * \overload
     * \brief Generate four integration points for a single tetrahedron.
     * \param [in] tet The nodes of the tetrahedron.
     * \param [out] ip_coords The coordinates of the tetrahedron's integration points.
     * \param [out] ip_weights The weights of the tetrahedron's integration points.
     * \return [void]
     */
    void GenerateFourPointsPerTetra(const std::vector<Node> &tet_nodes, std::vector<Vec3<double> > &ip_coords,
                                    std::vector<double> &ip_weights);


    /*!
     * \brief Generate five integration points per tetrahedron of a tetrahedral mesh.
     * \param [in] tetramesh The tetrahedral mesh.
     * \return [void]
     */
    void GenerateFivePointsPerTetra(const TetraMesh &tetramesh);


    /*!
     * \brief Generate adaptive integration points per tetrahedron of a tetrahedral mesh.
     * \param [in] tetramesh The tetrahedral mesh.
     * \param [in] options The integration options.
     * \param [in] support_dom The support domains of the tetrahedral mesh nodes.
     * \return [void]
     */
    void GenerateAdaptivePointsPerTetra(const TetraMesh &tetramesh, const IntegOptions &options,
                                        const InfSupportDomain &support_dom);

    /*!
     * \brief Adaptive integration with two divisions over a tetrahedron.
     * \param [in] tetramesh The tetrahedral mesh.
     * \param [in] tet_nodes The tetrahedron nodes.
     * \param [in] support_dom The support domains of the tetrahedral mesh nodes.
     * \param [in] adaptive_eps The minimum estimated approximate error. Use for adaptation termination.
     * \param [in] adaptive_level The level of recursion of the adaptive integration.
     * \param [out] tet_ip_coords The coordinates of the generated integration points of the tetrahedron.
     * \param [out] tet_ip_weights The weights of the generated integration points of the tetrahedron.
     * \return [void]
     */
    void AdaptiveTwoTetraDivisions(const TetraMesh &tetramesh, const std::vector<Node> &tet_nodes,
                                   const InfSupportDomain &support_dom, const double &adaptive_eps, int &adaptive_level,
                                   std::vector<Vec3<double> > &tet_ip_coords, std::vector<double> &tet_ip_weights);


    /*!
     * \brief Adaptive integration with four divisions over a tetrahedron.
     * \param [in] tetramesh The tetrahedral mesh.
     * \param [in] tet_nodes The tetrahedron nodes.
     * \param [in] support_dom The support domains of the tetrahedral mesh nodes.
     * \param [in] adaptive_eps The minimum estimated approximate error. Use for adaptation termination.
     * \param [in] adaptive_level The level of recursion of the adaptive integration.
     * \param [out] tet_ip_coords The coordinates of the generated integration points of the tetrahedron.
     * \param [out] tet_ip_weights The weights of the generated integration points of the tetrahedron.
     * \return [void]
     * \todo support_dom should be passed as rvalue reference.
     */
    void AdaptiveFourTetraDivisions(const TetraMesh &tetramesh, const std::vector<Node> &tet_nodes,
                                    const InfSupportDomain &support_dom, const double &adaptive_eps, int &adaptive_level,
                                    std::vector<Vec3<double> > &tet_ip_coords, std::vector<double> &tet_ip_weights);


    /*!
     * \brief Adaptive integration with eight divisions over a tetrahedron.
     * \param [in] tetramesh The tetrahedral mesh.
     * \param [in] tet_nodes The tetrahedron nodes.
     * \param [in] support_dom The support domains of the tetrahedral mesh nodes.
     * \param [in] adaptive_eps The minimum estimated approximate error. Use for adaptation termination.
     * \param [in] adaptive_level The level of recursion of the adaptive integration.
     * \param [out] tet_ip_coords The coordinates of the generated integration points of the tetrahedron.
     * \param [out] tet_ip_weights The weights of the generated integration points of the tetrahedron.
     * \return [void]
     * \todo support_dom should be passed as rvalue reference.
     */
    void AdaptiveEightTetraDivisions(const TetraMesh &tetramesh, const std::vector<Node> &tet_nodes,
                                     const InfSupportDomain &support_dom, const double &adaptive_eps,
                                     std::vector<Vec3<double> > &tet_ip_coords, std::vector<double> &tet_ip_weights);


    /*!
     * \brief Calculate the integration value.
     * \param [in] tetramesh The tetrahedral mesh.
     * \param [in] ip_coords The coordinates of the integration points.
     * \param [in] ip_weights The weights of the integration points.
     * \param [in] neigh_ids The indices of neighbor tetrahedral nodes to the integration points.
     * \param [in] support_dom The support (influence) domains of the tetrahedral mesh nodes. Used for adaptive integration.
     * \return [CLOUDEA::Vec3<double>] The integration value.
     */
    Vec3<double> IntegrationValue(const TetraMesh &tetramesh, const std::vector<Vec3<double> > &ip_coords,
                                  const std::vector<double> &ip_weights, const InfSupportDomain &support_dom);

private:

    std::vector<Vec3<double> > coordinates_;      /*!< The coordinates of the integration points. */

    std::vector<double> weights_;                 /*!< The weight of the integration points. */
};


/*! @} End of Doxygen Groups*/
} //end of namespace CLOUDEA

#endif //CLOUDEA_INTEGRATION_INTEG_POINTS_HPP_
