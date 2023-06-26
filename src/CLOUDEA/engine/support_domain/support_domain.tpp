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

#ifndef CLOUDEA_SUPPORT_DOMAIN_SUPPORT_DOMAIN_TPP_
#define CLOUDEA_SUPPORT_DOMAIN_SUPPORT_DOMAIN_TPP_


#include "CLOUDEA/engine/support_domain/support_domain.hpp"

namespace CLOUDEA {


template<short DIM>
void SupportDomain<DIM>::ExhaustiveInfluenceNodesSearch(const std::vector<IMP::Vec<DIM, double>> &points, const std::vector<IMP::Vec<DIM, double>> &field_nodes)
{
    // Check nodes number consistency.
    if (this->field_nodes_num_ != static_cast<int>(field_nodes.size())) {
        std::string error_str = "Could not perform exhaustive search for support nodes. " 
                                "The support domain nodes number does not match with the size of the given nodes.";
        
        throw std::invalid_argument(Logger::Error(error_str, __FILE__, __LINE__));
    }

    // Re initialize support nodes indices containers.
    this->influence_node_ids_.clear();
    this->influence_node_ids_.resize(points.size(), std::vector<int>());

    // Reset min and max support domain nodes number.
    this->min_influence_nodes_num_ = std::numeric_limits<int>::max();
    this->max_influence_nodes_num_ = 0;

    // Iterate over the points nodes.
    for (const auto &point : points) {

        auto pid = &point - &points[0];

        // Iterate over the field nodes to find neighbors.
        bool in_support;
        for (const auto &neigh_node : field_nodes) {

            auto neigh_id = &neigh_node - &field_nodes[0];

            // Initialy consider neigh as part of the influence domain.
            in_support = true;

            // Check if neighbour nodes is out of the region defined by the influence domain's dilated radius.
            if (std::sqrt(point.Distance2(neigh_node)) > (this->dilate_coeff_[neigh_id]*this->radius_[neigh_id])) {
                in_support = false;
            }

            // Add the index of the neighbor id in the influence domain of the point.
            if (in_support) {
                this->influence_node_ids_[pid].emplace_back(neigh_id);
            }

        } // End of Iterate over the nodes to find neighbors.

        // Update min - max support domain nodes number.
        if (static_cast<int>(influence_node_ids_[pid].size()) > this->max_influence_nodes_num_) { this->max_influence_nodes_num_ = static_cast<int>(influence_node_ids_[pid].size()); }
        if (static_cast<int>(influence_node_ids_[pid].size()) < this->min_influence_nodes_num_) { this->min_influence_nodes_num_ = static_cast<int>(influence_node_ids_[pid].size()); }
    

    } // End of Iterate over the points.
}


//#ifdef CLOUDEA_WITH_CGAL
template<short DIM>
void SupportDomain<DIM>::FastInfluenceNodesSearch_2d(const std::vector<IMP::Vec<DIM, double>> &points, const std::vector<IMP::Vec<DIM, double>> &field_nodes)
{
    // CGAL types declaration.
    typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
    typedef K::Point_2 Point;
    typedef boost::tuple<Point, int> Point_and_int;
    typedef CGAL::Search_traits_2<K> Traits_base;
    typedef CGAL::Search_traits_adapter<Point_and_int, CGAL::Nth_of_tuple_property_map<0, Point_and_int>, Traits_base> Traits;
    typedef CGAL::Kd_tree<Traits> Tree;
    typedef CGAL::Fuzzy_sphere<Traits> Fuzzy_sphere;

    // Check nodes number consistency.
    if (this->field_nodes_num_ != static_cast<int>(field_nodes.size())) {
        std::string error_str = "Could not perform fast search for support nodes. " 
                                "The support domain nodes number does not match with the size of the given nodes.";
        
        throw std::invalid_argument(Logger::Error(error_str, __FILE__, __LINE__));
    }

    // Re initialize support nodes indices containers.
    this->influence_node_ids_.clear();
    this->influence_node_ids_.resize(points.size(), std::vector<int>());

    // Initialize tree tuple points coordinates.
    std::vector<Point> points_coords;
    points_coords.reserve(points.size());

    // Initialize tree tuple points indices.
    std::vector<int> points_ids;
    points_ids.reserve(points.size());

    auto pid = 0;
    for (auto &point : points) {
        points_coords.emplace_back(Point(point[0], point[1]));
        points_ids.emplace_back(pid++);
    }
    
    // Insert <point,id> tuples in the searching tree.
    Tree tree(boost::make_zip_iterator(boost::make_tuple(points_coords.begin(), points_ids.begin())),
              boost::make_zip_iterator(boost::make_tuple(points_coords.end(), points_ids.end())));    
    
    // Find the influence field nodes for influenced point.
    std::vector<Point_and_int> neighs;
    for (auto &field_node : field_nodes) {
        auto fn_id = &field_node - &field_nodes[0];

        // Set field node as the influence domain's center.
        Point influence_center(field_node[0], field_node[1]);

        // Set the influence domain range.
        Fuzzy_sphere influence_range(influence_center, this->dilate_coeff_[fn_id]*this->radius_[fn_id]);

        // Get the influence domain nodes.
        tree.search(std::back_inserter(neighs), influence_range);

        // Establish the domains as influence or influence domain.
        for (const auto &neigh : neighs) {
                // Store the influencing field node index in the support of the neighbor (influenced) node.
                this->influence_node_ids_[boost::get<1>(neigh)].emplace_back(fn_id);
        }

        // Clear support nodes vector for next iteration.
        neighs.clear();
    }

    // Update min - max support domain nodes number.
    std::size_t min_sd = std::numeric_limits<std::size_t>::max();
    std::size_t max_sd = 0;
    for (const auto &sd_nodes : this->influence_node_ids_) {
        if (sd_nodes.size() > max_sd) { max_sd = sd_nodes.size(); }
        if (sd_nodes.size() < min_sd) { min_sd = sd_nodes.size(); }
    }

    this->min_influence_nodes_num_ = static_cast<int>(min_sd);
    this->max_influence_nodes_num_ = static_cast<int>(max_sd);

}


template<short DIM>
void SupportDomain<DIM>::FastInfluenceNodesSearch_3d(const std::vector<IMP::Vec<DIM, double>> &points, const std::vector<IMP::Vec<DIM, double>> &field_nodes)
{
    // CGAL types declaration.
    typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
    typedef K::Point_3 Point;
    typedef boost::tuple<Point, int> Point_and_int;
    typedef CGAL::Search_traits_3<K> Traits_base;
    typedef CGAL::Search_traits_adapter<Point_and_int, CGAL::Nth_of_tuple_property_map<0, Point_and_int>, Traits_base> Traits;
    typedef CGAL::Kd_tree<Traits> Tree;
    typedef CGAL::Fuzzy_sphere<Traits> Fuzzy_sphere;

    // Check nodes number consistency.
    if (this->field_nodes_num_ != static_cast<int>(field_nodes.size())) {
        std::string error_str = "Could not perform fast search for support nodes. " 
                                "The support domain nodes number does not match with the size of the given nodes.";
        
        throw std::invalid_argument(Logger::Error(error_str, __FILE__, __LINE__));
    }

    // Reset support nodes indices containers.
    this->influence_node_ids_.clear();
    this->influence_node_ids_.resize(points.size(), std::vector<int>());

    // Initialize tree tuple points coordinates.
    std::vector<Point> points_coords;
    points_coords.reserve(points.size());

    // Initialize tree tuple points indices.
    std::vector<int> points_ids;
    points_ids.reserve(points.size());

    for (const auto &point : points) {
        auto pid = &point - &points[0];

        points_coords.emplace_back(Point(point[0], point[1], point[2]));
        points_ids.emplace_back(pid);
    }
    
    // Insert <point,id> tuples in the searching tree.
    Tree tree(boost::make_zip_iterator(boost::make_tuple(points_coords.begin(), points_ids.begin())),
              boost::make_zip_iterator(boost::make_tuple(points_coords.end(), points_ids.end())));   
    
    // Find the influence nodes for each influenced point.
    std::vector<Point_and_int> neighs;
    for (const auto &field_node : field_nodes) {
        auto field_nid = &field_node - &field_nodes[0];

        // Set field node as the influence domain's center.
        Point influence_center(field_node[0], field_node[1], field_node[2]);

        // Set the influence domain range.
        Fuzzy_sphere influence_range(influence_center, this->dilate_coeff_[field_nid]*this->radius_[field_nid]);

        // Get the influence domain nodes.
        tree.search(std::back_inserter(neighs), influence_range);
        
        // Add the node id to the influence domain of the influenced points.
        for (const auto &neigh : neighs) {
            this->influence_node_ids_[boost::get<1>(neigh)].emplace_back(field_nid);
        }
        
        // Clear neighbor nodes vector for next iteration.
        neighs.clear();
    }

    // Update min - max influence domain nodes number.
    std::size_t min_sd = std::numeric_limits<std::size_t>::max();
    std::size_t max_sd = 0;
    for (const auto &sd_nodes : this->influence_node_ids_) {
        if (sd_nodes.size() > max_sd) { max_sd = sd_nodes.size(); }
        if (sd_nodes.size() < min_sd) { min_sd = sd_nodes.size(); }
    }
    this->min_influence_nodes_num_ = static_cast<int>(min_sd);
    this->max_influence_nodes_num_ = static_cast<int>(max_sd);
}


template<short DIM>
void SupportDomain<DIM>::FastSearchNearestInfluenceNodes_3d(const std::vector<IMP::Vec<DIM, double>> &field_nodes, const IMP::NodeSet &surf_nodeset, int neigh_num)
{
    if (neigh_num < 8) {
        std::string error_str = "Could not compute nearest influence nodes. Too few nearest neighbors requested. Provide a number larger than 7.";
        throw std::invalid_argument(Logger::Error(error_str));
    }

    typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
    typedef K::Point_3                                          Point;
    typedef boost::tuple<Point,int>                             Point_and_int;
    typedef CGAL::Search_traits_3<K>                            Traits_base;
    typedef CGAL::Search_traits_adapter<Point_and_int,CGAL::Nth_of_tuple_property_map<0, Point_and_int>, Traits_base> Traits;
    typedef CGAL::Orthogonal_k_neighbor_search<Traits>          K_neighbor_search;
    typedef K_neighbor_search::Tree                             Tree;
    typedef K_neighbor_search::Distance                         Distance;

    // Reset support domain influence nodes indices.
    this->influence_node_ids_.clear();
    this->influence_node_ids_.resize(field_nodes.size(), std::vector<int>());

    // Reset support domain radius.
    this->radius_.clear();
    this->radius_.resize(field_nodes.size(), 0.);

    // Reset support domain dilatation coefficient.
    this->dilate_coeff_.clear();
    this->dilate_coeff_.resize(field_nodes.size(), 0.);

    // Distinguish boundary from internal field nodes.
    std::vector<bool> onsurf_flag(field_nodes.size(), false);
    for (const auto &nid : surf_nodeset.NodeIds()) {
        onsurf_flag[nid] = true;
    }

    // Initialize tree structure of the influence field nodes of the meshfree grid.
    std::vector<Point> points, onsurf_points, inside_points;
    std::vector<int> point_ids, onsurf_point_ids, inside_point_ids;
    for (const auto &node : field_nodes) {
        auto nid = &node - &field_nodes[0];

        // Store in the full neighbor points containers.
        points.emplace_back(Point(node[0], node[1], node[2]));
        point_ids.emplace_back(nid);

        // Store in the partial neighbor points containers.
        if (onsurf_flag[nid]) {
            onsurf_points.emplace_back(Point(node[0], node[1], node[2]));
            onsurf_point_ids.emplace_back(nid);
        } else {
            inside_points.emplace_back(Point(node[0], node[1], node[2]));
            inside_point_ids.emplace_back(nid);
        }
    }

    // Compute partial neighbor nodes numbers.
    int onsurf_neigh_num = neigh_num/3;
    int inside_neigh_num = neigh_num - onsurf_neigh_num;
    
    // Construct neighbors tree structures.
    Tree full_tree(boost::make_zip_iterator(boost::make_tuple(std::begin(points), std::begin(point_ids))),
                   boost::make_zip_iterator(boost::make_tuple(std::end(points), std::end(point_ids))) );

    Tree onsurf_tree(boost::make_zip_iterator(boost::make_tuple(std::begin(onsurf_points), std::begin(onsurf_point_ids))),
                     boost::make_zip_iterator(boost::make_tuple(std::end(onsurf_points), std::end(onsurf_point_ids))) );

    Tree inside_tree(boost::make_zip_iterator(boost::make_tuple(std::begin(inside_points), std::begin(inside_point_ids))),
                     boost::make_zip_iterator(boost::make_tuple(std::end(inside_points), std::end(inside_point_ids))) );

    // Search for the nearest influence nodes to each point.
    Distance tr_dist;
    for (const auto &query : points) {
        auto id = &query - &points[0];

        if (onsurf_flag[id]) {
            // Get the nearest influence node for the query point from the surface.
            K_neighbor_search onsurf_search(onsurf_tree, query, onsurf_neigh_num);

            // Get the nearest influence node for the query point from inside.
            K_neighbor_search inside_search(inside_tree, query, inside_neigh_num);

            for (K_neighbor_search::iterator it = onsurf_search.begin(); it != onsurf_search.end(); it++) {
                this->influence_node_ids_[id].emplace_back(boost::get<1>(it->first));
            }

            for (K_neighbor_search::iterator it = inside_search.begin(); it != inside_search.end(); it++) {
                this->influence_node_ids_[id].emplace_back(boost::get<1>(it->first));
            }

            // Set support domain radius and dilatation coefficient according to nearest neighbors.
            double min_neigh_dist = std::sqrt(tr_dist.transformed_distance((onsurf_search.begin()+1)->second));
            double max_neigh_dist = std::sqrt(tr_dist.transformed_distance((onsurf_search.end()-1)->second));

            if (min_neigh_dist > std::sqrt(tr_dist.transformed_distance((inside_search.begin()+1)->second))) {
                min_neigh_dist = std::sqrt(tr_dist.transformed_distance((inside_search.begin()+1)->second));
            }

            if (max_neigh_dist < std::sqrt(tr_dist.transformed_distance((inside_search.end()-1)->second))) {
                max_neigh_dist = std::sqrt(tr_dist.transformed_distance((inside_search.end()-1)->second));
            }

            this->radius_[id] = min_neigh_dist;
            this->dilate_coeff_[id] = max_neigh_dist/min_neigh_dist;

        } else {
            // Get the nearest influence node for the query point.
            K_neighbor_search full_search(full_tree, query, neigh_num);
            for (K_neighbor_search::iterator it = full_search.begin(); it != full_search.end(); it++) {
                this->influence_node_ids_[id].emplace_back(boost::get<1>(it->first));
            }

            // Set support domain radius and dilatation coefficient according to nearest neighbors.
            double min_neigh_dist = std::sqrt(tr_dist.transformed_distance((full_search.begin()+1)->second));
            double max_neigh_dist = std::sqrt(tr_dist.transformed_distance((full_search.end()-1)->second));

            this->radius_[id] = min_neigh_dist;
            this->dilate_coeff_[id] = max_neigh_dist/min_neigh_dist;
        }
    }

    // Update min - max influence domain nodes number.
    this->min_influence_nodes_num_ = neigh_num;
    this->max_influence_nodes_num_ = neigh_num;

}
//#endif


template<short DIM>
SupportDomain<DIM>::SupportDomain() : influence_node_ids_(), radius_(), dilate_coeff_(), field_nodes_num_(0),
    min_influence_nodes_num_(std::numeric_limits<int>::max()), max_influence_nodes_num_(0)
{}


template<short DIM>
SupportDomain<DIM>::~SupportDomain()
{}


template<short DIM>
void SupportDomain<DIM>::SetFieldNodesNum(int field_nodes_num)
{
    if (field_nodes_num < 0) {
        std::string error_str = "Could not set field nodes number in SupportDomain. A negative number was given.";
        throw std::invalid_argument(Logger::Error(error_str, __FILE__, __LINE__));
    }

    this->field_nodes_num_ = field_nodes_num;
    
}


template<short DIM>
void SupportDomain<DIM>::SetRadius(IMP::Vec<DIM, double> nodal_dist)
{
    if (this->field_nodes_num_ < 1) {
        std::string error_str = "Could not set support domain radius. Set the number of nodes first.";
        throw std::runtime_error(Logger::Error(error_str, __FILE__, __LINE__));
    }

    this->radius_.clear();
    this->radius_.resize(this->field_nodes_num_, std::sqrt(nodal_dist.Length2()));
}


template<short DIM>
void SupportDomain<DIM>::SetRadius(double radius)
{
    if (this->field_nodes_num_ < 1) {
        std::string error_str = "Could not set support domain radius. Set the number of nodes first.";
        throw std::runtime_error(Logger::Error(error_str, __FILE__, __LINE__));
    }

    this->radius_.clear();
    this->radius_.resize(this->field_nodes_num_, radius);
}


template<short DIM>
void SupportDomain<DIM>::SetDilateCoeff(std::vector<double> dilate_coeff)
{
    if (this->field_nodes_num_ < 1) {
        std::string error_str = "Could not set support domain radius. Set the number of nodes first.";
        throw std::runtime_error(Logger::Error(error_str, __FILE__, __LINE__));
    }

    if (this->field_nodes_num_ != static_cast<int>(dilate_coeff.size()) && dilate_coeff.size() != 1) {
        std::string error_str = "Could not set support domain dilatation coefficient. Not provided for all the field nodes.";
        throw std::runtime_error(Logger::Error(error_str, __FILE__, __LINE__));
    }

    if (dilate_coeff.size() == 1) {
        this->dilate_coeff_.clear();
        this->dilate_coeff_.resize(this->field_nodes_num_, dilate_coeff[0]);
    } else {
        this->dilate_coeff_ = dilate_coeff;
    }

}


template<short DIM>
void SupportDomain<DIM>::SetDilateCoeff(double dilate_coeff)
{
    if (this->field_nodes_num_ < 1) {
        std::string error_str = "Could not set support domain dilatation coefficient. Set the number of nodes first.";
        throw std::runtime_error(Logger::Error(error_str, __FILE__, __LINE__));
    }

    if (dilate_coeff < 2*std::numeric_limits<double>::epsilon()) {
        std::string error_str = "Could not set support domain dilatation coefficient. A positive value is required.";
        throw std::runtime_error(Logger::Error(error_str, __FILE__, __LINE__));
    }

    this->dilate_coeff_.clear();
    this->dilate_coeff_.resize(this->field_nodes_num_, dilate_coeff);

}


template<short DIM>
void SupportDomain<DIM>::ComputeRadiusFromRegularGrid(const std::vector<IMP::Vec<DIM, double>> &field_nodes)
{
    // Check nodes number consistency.
    if (this->field_nodes_num_ != static_cast<int>(field_nodes.size())) {
        std::string error_str = "Could not compute support domains' radius from regular grid. " 
                                "The support domain nodes number does not match with the number of the given field nodes.";
        
        throw std::invalid_argument(Logger::Error(error_str, __FILE__, __LINE__));
    }

    // Initialize nodal spacing to maximum possible value.
    double spacing = std::numeric_limits<double>::max();
    
    // Get the minimum spacing between nodes for each spatial dimension.
    double internode_spacing = 0.;
    for (std::size_t i = 1; i != field_nodes.size(); ++i) {

        internode_spacing = std::sqrt(field_nodes[i].Distance2(field_nodes[i-1]));

        if (internode_spacing < spacing && internode_spacing > 2.*std::numeric_limits<double>::epsilon()) { spacing = internode_spacing; }
    }

    // Set the support radius for all field nodes.
    this->SetRadius(spacing);

}


template<short DIM>
template<short CELL_NODES>
void SupportDomain<DIM>::ComputeRadiusFromIrregularGrid(const IMP::Grid<DIM, CELL_NODES> &grid)
{

    // Check nodes number consistency.
    if (this->field_nodes_num_ != grid.NodesNum()) {
        std::string error_str = "Could not compute support domains' radius from irregular grid. " 
                                "The support domain nodes number does not match with the number of the given field nodes.";
        
        throw std::invalid_argument(Logger::Error(error_str, __FILE__, __LINE__));
    }

    // Check that cell is supported.
    if (CELL_NODES > 4) {
        std::string error_str = "Could not compute support domains' radius from irregular grid. Not support ghost cell type.";        
        throw std::invalid_argument(Logger::Error(error_str, __FILE__, __LINE__));
    }


    // Clear the influence radiuses container.
    this->radius_.clear();
    this->radius_.resize(this->field_nodes_num_, 0.);

    // The number of connected nodes to each influence node.
    std::vector<int> connected_nodes_number(this->field_nodes_num_, 0);

    // Extract edges of triangle or quadrilateral cells.
    auto ExtractEdges2D = [] (const IMP::Cell<DIM, CELL_NODES> &cell) { 
        std::vector<IMP::Cell<DIM,2>> edges(CELL_NODES, IMP::Cell<DIM,2>());
        for (int i = 0; i != CELL_NODES-1; ++i) {
            edges[i].SetConnectivity({cell.N(i), cell.N(i+1)});
        }
        edges[CELL_NODES-1].SetConnectivity({cell.N(CELL_NODES-1), cell.N(0)});
        return edges;
    };

    // Extract edges of tetrahedral cell.
    auto ExtractEdgesTet = [] (const IMP::Cell<DIM, CELL_NODES> &cell) { 
        std::vector<IMP::Cell<DIM,2>> edges(6, IMP::Cell<DIM,2>());
        int pos = 0;
        for (int i = 0; i != 3; ++i) {
            for (int j = i+1; j != 4; ++j) {
                edges[pos].SetConnectivity({cell.N(i), cell.N(j)});
                pos++;
            }
        }
        return edges;
    };

    // Iterate over the cells of the grid.
    std::vector<IMP::Cell<DIM,2>> cell_edges;
    for (const auto &cell : grid.GhostCells()){

        // Get the edges of the cell.
        if (CELL_NODES == 2) {
            IMP::Cell<DIM,2> edge;
            edge.SetConnectivity({cell.N(0), cell.N(1)});
            cell_edges.clear();
            cell_edges.emplace_back(edge);
        }
        else if (DIM == 2 && CELL_NODES > 2) {
            cell_edges = ExtractEdges2D(cell);
        } else {
            cell_edges = ExtractEdgesTet(cell);
        }

        // Iterate over the cell's edges
        for (const auto &edge : cell_edges) {
            double dist = std::sqrt( grid.Nodes(edge.N(0)).Distance2(grid.Nodes(edge.N(1))) );

            // Add the distance value to the support of the influence nodes at the corners of the edge.
            this->radius_[edge.N(0)] += dist;
            this->radius_[edge.N(1)] += dist;
            
            // Increase the number of connected nodes to each of the nodes at the corners of the edge.
            connected_nodes_number[edge.N(0)] += 1;
            connected_nodes_number[edge.N(1)] += 1;
        }
    }

    // Normalize the influence radiuses with the number of connected nodes for each influence node.
    for (int id = 0; id != this->field_nodes_num_; ++id) {
        // Normalize the influence radius of the current influence node.
        this->radius_[id] /= connected_nodes_number[id];
    }

}


template<short DIM>
template<short CELL_NODES>
void SupportDomain<DIM>::ComputeRadiusFromImmersedGrid(const IMP::Grid<DIM, CELL_NODES> &grid)
{
    // Check nodes number consistency.
    if (this->field_nodes_num_ != grid.NodesNum()) {
        std::string error_str = "Could not compute support domains' radius from immersed grid. " 
                                "The support domain nodes number does not match with the number of the given field nodes.";
        
        throw std::invalid_argument(Logger::Error(error_str, __FILE__, __LINE__));
    }

    // Clear the influence radiuses container.
    this->radius_.clear();
    this->radius_.resize(this->field_nodes_num_, 0.);

    // The number of connected nodes to each influence node.
    std::vector<int> connected_nodes_number(this->field_nodes_num_, 0);

    // Assign flag=1 to surface nodes.
    std::vector<int> flags(this->field_nodes_num_, 0);
    for (const auto &cell : grid.GhostCells()){
        for (const auto &node_id : cell.Connectivity()) {
            flags[node_id] = 1;
        }
    }

    // Collect interior nodes.
    std::vector<IMP::Vec<DIM, double>> inter_nodes;
    inter_nodes.reserve(this->field_nodes_num_);
    for (const auto &flag : flags) {
        auto id = &flag - &flags[0];

        if (flag == 0) {
            inter_nodes.emplace_back(grid.Nodes(id));
        }
    }
    inter_nodes.shrink_to_fit();

    flags.clear();
    flags.shrink_to_fit();

    // Initialize interior nodal spacing to maximum possible value.
    double spacing = std::numeric_limits<double>::max();
    
    // Get the minimum spacing between interior nodes for each spatial dimension.
    double internode_spacing = 0.;
    for (std::size_t i = 1; i != inter_nodes.size(); ++i) {
        internode_spacing = std::sqrt(inter_nodes[i].Distance2(inter_nodes[i-1]));
        if (internode_spacing < spacing && internode_spacing > 2*std::numeric_limits<double>::epsilon()) { spacing = internode_spacing; }
    }

    // Set the support radius for all field nodes.
    this->SetRadius(spacing);

}


template<short DIM>
void SupportDomain<DIM>::IdentifyInfluenceNodesInRange(const std::vector<IMP::Vec<DIM, double>> &points, const std::vector<IMP::Vec<DIM, double>> &field_nodes)
{
    //#ifdef CLOUDEA_WITH_CGAL
        if (DIM == 1) {
            this->ExhaustiveInfluenceNodesSearch(points, field_nodes);
        }
        else if (DIM == 2) {
            this->FastInfluenceNodesSearch_2d(points, field_nodes);
        }
        else if (DIM == 3) {
            this->FastInfluenceNodesSearch_3d(points, field_nodes);
        } else {
            std::string error_str = "Could not identify support nodes in support domain. Supported domain dimensions [1-3].";
            throw std::runtime_error(Logger::Error(error_str, __FILE__, __LINE__));
        }
    // #else
    //     if (DIM == 2 || DIM == 3) {
    //         std::cout << Logger::Warning("Support nodes identification in support domain runs with exhaustive search. Build CLOUDEA with CGAL for higher efficiency.\n");
    //     }
    //     else if (DIM < 1 && DIM > 3) {
    //         std::string error_str = "Could not identify support nodes in support domain.  Supported domain dimensions [1-3].";
    //         throw std::runtime_error(Logger::Error(error_str, __FILE__, __LINE__));
    //     }

    //     this->ExhaustiveInfluenceNodesSearch(nodes);
    
    // #endif

}


template<short DIM>
void SupportDomain<DIM>::IdentifyNearestInfluenceNodes(const std::vector<IMP::Vec<DIM, double>> &field_nodes, const IMP::NodeSet &surf_nodeset, int neigh_num)
{
        if (DIM == 1) {
            this->ExhaustiveInfluenceNodesSearch(field_nodes, field_nodes);
        }
        else if (DIM == 2) {
            this->FastInfluenceNodesSearch_2d(field_nodes, field_nodes);
        }
        else if (DIM == 3) {
            this->FastSearchNearestInfluenceNodes_3d(field_nodes, surf_nodeset, neigh_num);
        } else {
            std::string error_str = "Could not identify support nodes in support domain. Supported domain dimensions [1-3].";
            throw std::runtime_error(Logger::Error(error_str, __FILE__, __LINE__));
        }

}


template<short DIM>
void SupportDomain<DIM>::IdentifyInfluenceNodesVoronoi(const IMP::Voronoi<DIM> &voronoi)
{
    if (voronoi.FacetsNum() == 0) {
        std::string error_str = "Could not identify influence nodes from voronoi tesselation. Extract the facets of the voronoi tesselation first.";
        throw std::runtime_error(Logger::Error(error_str, __FILE__, __LINE__));
    }

    // Set field nodes number.
    this->field_nodes_num_ = voronoi.NodesNum();

    // Reset support domain influence nodes indices.
    this->influence_node_ids_.clear();
    this->influence_node_ids_.resize(this->field_nodes_num_, std::vector<int>());

    std::vector<std::vector<int>> neigh_nodes(this->field_nodes_num_, std::vector<int>());
    
    int n1 = 0, n2 = 0;
    for (const auto &facet : voronoi.Facets()) {
        if (!facet.IsFree()) {
            // Get node indices sharing the internal facet.
            n1 = facet.ParentCellId();
            n2 = facet.NeighCellId();
            
            // Add each node as neighbor to the other.
            neigh_nodes[n1].emplace_back(n2);
            neigh_nodes[n2].emplace_back(n1);
        }
    }

    int id = 0;
    for (auto &neighs : neigh_nodes) {
        // Sort the neighs indices.
        std::sort(neighs.begin(), neighs.end());

        // Store the influence nodes indices. First influence node is the nth node itself.
        this->influence_node_ids_[id].emplace_back(id);
        for (const auto &n : neighs) { this->influence_node_ids_[id].emplace_back(n); }
        id++;
    }

    // Update min - max influence domain nodes number.
    std::size_t min_sd = std::numeric_limits<std::size_t>::max();
    std::size_t max_sd = 0;
    for (const auto &inf_nodes : this->influence_node_ids_) {
        if (inf_nodes.size() > max_sd)  max_sd = inf_nodes.size(); 
        if (inf_nodes.size() < min_sd)  min_sd = inf_nodes.size(); 
    }
    this->min_influence_nodes_num_ = static_cast<int>(min_sd);
    this->max_influence_nodes_num_ = static_cast<int>(max_sd);

}


template<short DIM>
void SupportDomain<DIM>::PrintInfluenceNodesAndRadius(const std::string &output_file) const {
    
    std::ofstream out(output_file, std::ios::out);

    out << "# File format: Radius n_1 n_2 .... n_n, where 1-n are the n support node indices of the field node.";
    out << "# Data per field node are seperated with a empty line.\n#\n";

    for (std::size_t i = 0; i != this->radius_.size(); ++i) {
        out << this->dilate_coeff_[i] * this->radius_[i] << " ";

        for (const auto &node_id : this->influence_node_ids_[i]) {
            out << node_id << " ";
        }

        out << "\n\n";
    }

    out.close();
}



} // End of namespace CLOUDEA.

#endif //CLOUDEA_SUPPORT_DOMAIN_SUPPORT_DOMAIN_TPP_