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
   \file support_domain.hpp
   \brief SupportDomain class header file.
   \author Konstantinos A. Mountris
   \date 15/11/2019
*/

#ifndef CLOUDEA_SUPPORT_DOMAIN_SUPPORT_DOMAIN_HPP_
#define CLOUDEA_SUPPORT_DOMAIN_SUPPORT_DOMAIN_HPP_


#include "CLOUDEA/engine/utilities/logger.hpp"

#include <IMP/Vectors>
#include <IMP/Tesselations>

//#ifdef CLOUDEA_WITH_CGAL
    #include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
    #include <CGAL/Simple_cartesian.h>
    #include <CGAL/Kd_tree.h>
    #include <CGAL/algorithm.h>
    #include <CGAL/Fuzzy_sphere.h>
    #include <CGAL/Search_traits_2.h>
    #include <CGAL/Search_traits_3.h>
    #include <CGAL/Search_traits_adapter.h>
    #include <CGAL/property_map.h>
    #include <CGAL/Orthogonal_k_neighbor_search.h>
    #include <boost/iterator/zip_iterator.hpp>
//#endif

#include <cmath>
#include <cstddef>
#include <string>
#include <vector>
#include <sstream>
#include <fstream>
#include <stdexcept>
#include <exception>
#include <limits>
#include <unordered_set>
#include <algorithm>
#include <utility>

namespace CLOUDEA {

/** \addtogroup Meshfree \{ */


/**
 * \class SupportDomain
 * \author Konstantinos A. Mountris
 * \brief Construction of support domain for field nodes in meshfree-based models.
 * \tparam DIM The spatial dimensions of the support domain.
 */
template <short DIM>
class SupportDomain
{

private:

    std::vector<std::vector<int>> influence_node_ids_;     /**< The indices of the support nodes for each field node */

    std::vector<double> radius_;                           /**< The radius of the support domain for each field node */
    
    std::vector<double> dilate_coeff_;                     /**< The dilatation coefficient of the support domain radius */

    int field_nodes_num_;                                  /**< The number of field nodes for which support domains are constructed */

    int min_influence_nodes_num_;                          /**< The number of support nodes in the smallest support domain */

    int max_influence_nodes_num_;                          /**< The number of support nodes in the largest support domain */


protected:

    /**
     * \brief 
     * \param nodes 
     */
    inline void ExhaustiveInfluenceNodesSearch(const std::vector<IMP::Vec<DIM, double>> &points, const std::vector<IMP::Vec<DIM, double>> &field_nodes);

    //#ifdef CLOUDEA_WITH_CGAL

    /**
     * @brief 
     * 
     * @param nodes 
     */
    inline void FastInfluenceNodesSearch_2d(const std::vector<IMP::Vec<DIM, double>> &points, const std::vector<IMP::Vec<DIM, double>> &field_nodes);
    
    /**
     * @brief 
     * 
     * @param nodes 
     */
    inline void FastInfluenceNodesSearch_3d(const std::vector<IMP::Vec<DIM, double>> &points, const std::vector<IMP::Vec<DIM, double>> &field_nodes);


    /**
     * @brief 
     * @param points 
     * @param field_nodes 
     */
    inline void FastSearchNearestInfluenceNodes_3d(const std::vector<IMP::Vec<DIM, double>> &field_nodes, const IMP::NodeSet &surf_nodeset, int neigh_num);
    //#endif

public:

    /**
     * \brief SupportDomain constructor.
     */
    SupportDomain();


    /**
     * \brief SupportDomain destructor.
     */
    virtual ~SupportDomain();


    /**
     * \brief Set the number of the field nodes for which support domains will be constructed.
     * \param [in] nodes_num The number of field nodes.
     * \return [void] 
     */
    inline void SetFieldNodesNum(int field_nodes_num);
    

    /**
     * \brief Set a global support domain radius for all field nodes.
     * It computes the radius magnitude from the nodes' distance. To be used when the distance between field nodes is constant.
     * \param [in] nodal_dist The distance between field nodes.
     * \return [void]
     */
    inline void SetRadius(IMP::Vec<DIM, double> nodal_dist);


    /**
     * \brief Set a global support domain radius for all field nodes.
     * \param [in] radius The common radius of all the field nodes.
     * \return [void]
     */
    inline void SetRadius(double radius);


    /**
     * \brief Set a dilatation coefficient for the support domain of each field node. May varying for each field node.
     * \param [in] dilate_coeff The dilatation coefficient of the support domain of each field node.
     * \return [void]
     */
    inline void SetDilateCoeff(std::vector<double> dilate_coeff);


    /**
     * \brief Set a constant dilatation coefficient for the support domain of each field node.
     * \param [in] dilate_coeff The dilatation coefficient of the support domain of each field node.
     * \return [void]
     */
    inline void SetDilateCoeff(double dilate_coeff);
    
    
    /**
     * \brief Compute the support radius from a regular grid of field nodes. Same support radius is assigned for all field nodes.
     * \param [in] field_nodes The field nodes forming a regular grid for which the support domain radius will be computed.
     * \return [void]
     */
    inline void ComputeRadiusFromRegularGrid(const std::vector<IMP::Vec<DIM, double>> &field_nodes);


    /**
     * \brief 
     * \tparam CELL_NODES 
     * \param mesh 
     */
    template<short CELL_NODES>
    void ComputeRadiusFromIrregularGrid(const IMP::Grid<DIM, CELL_NODES> &grid);


    /**
     * \brief 
     * \tparam CELL_NODES 
     * \param mesh 
     */
    template<short CELL_NODES>
    void ComputeRadiusFromImmersedGrid(const IMP::Grid<DIM, CELL_NODES> &grid);

    
    /**
     * \brief Identify the support nodes for each of the field nodes.
     * \param [in] nodes The field nodes.
     * \return [void] 
     */
    inline void IdentifyInfluenceNodesInRange(const std::vector<IMP::Vec<DIM, double>> &points, const std::vector<IMP::Vec<DIM, double>> &field_nodes);


    /**
     * @brief 
     * @param points 
     * @param field_nodes 
     */
    inline void IdentifyNearestInfluenceNodes(const std::vector<IMP::Vec<DIM, double>> &field_nodes, const IMP::NodeSet &surf_nodeset, int neigh_num);


    /**
     * @brief 
     * @param points 
     * @param field_nodes 
     */
    inline void IdentifyInfluenceNodesVoronoi(const IMP::Voronoi<DIM> &voronoi);


    /**
     * \brief Print the support nodes indices and the support radius in a file.
     * \param [in] output_file 
     * \return [void]
     */
    inline void PrintInfluenceNodesAndRadius(const std::string &output_file) const;
    
    
    /**
     * \brief 
     * \return int 
     */
    inline int FieldNodesNum() const { return this->field_nodes_num_; }


    /**
     * \brief 
     * \return double 
     */
    inline const std::vector<double> & Radius() const { return this->radius_; }


    /**
     * \brief 
     * \return double 
     */
    inline double Radius(std::size_t id) const { return this->radius_[id]; }


    /**
     * \brief 
     * \return int 
     */
    inline double RadiusAt(std::size_t id) const { return this->radius_.at(id); }


    /**
     * \brief 
     * \return double 
     */
    inline const std::vector<double> & DilateCoeff() const { return this->dilate_coeff_; }


    /**
     * \brief 
     * \return double 
     */
    inline double DilateCoeff(std::size_t id) const { return this->dilate_coeff_[id]; }


    /**
     * \brief 
     * \return double 
     */
    inline double DilateCoeffAt(std::size_t id) const { return this->dilate_coeff_.at(id); }


    /**
     * \brief Get the indices of the nodes in the support domain for all field nodes. 
     * \return [const std::vector<std::vector<int>>&] The indices of the nodes in the support domain for all field nodes.
     */
    inline const std::vector<std::vector<int>> & InfluenceNodeIds() const { return this->influence_node_ids_; }


    /**
     * \brief Get the indices of the nodes in the support domain for a specified field node.
     * Fast access with no range check.
     * \return [const std::vector<int>&] The indices of the nodes in the support domain for a specified field node.
     */
    inline const std::vector<int> & InfluenceNodeIds(std::size_t id) const { return this->influence_node_ids_[id]; }


    /**
     * \brief Get the indices of the nodes in the support domain for a specified field node.
     * Slower access with range check.
     * \return [const std::vector<std::vector<int>>&] The indices of the nodes in the support domain for a specified field node.
     */
    inline const std::vector<int> & InfluenceNodeIdsAt(std::size_t id) const { return this->influence_node_ids_.at(id); }


    /**
     * \brief Get the number of support nodes contained in the smallest support domain.
     * \return [int] The number of support nodes contained in the smallest support domain.
     */
    inline int MinInfluenceNodesNum() const { return this->min_influence_nodes_num_; }


    /**
     * \brief Get the number of support nodes contained in the largest support domain.
     * \return [int] The number of support nodes contained in the largest support domain.
     */
    inline int MaxInfluenceNodesNum() const { return this->max_influence_nodes_num_; }

};


/*! \} End of Doxygen Groups*/
} //end of namespace CLOUDEA


#endif //CLOUDEA_SUPPORT_DOMAIN_SUPPORT_DOMAIN_HPP_

#include "CLOUDEA/engine/support_domain/support_domain.tpp"