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



#ifndef CLOUDEA_APPROXIMANTS_FEM_D2V4_TPP_
#define CLOUDEA_APPROXIMANTS_FEM_D2V4_TPP_


#include "CLOUDEA/engine/approximants/fem_d2v4.hpp"


namespace CLOUDEA {

template<short DIM>
FemD2v4<DIM>::FemD2v4()
{}


template<short DIM>
FemD2v4<DIM>::~FemD2v4()
{}


template<short DIM>
void FemD2v4<DIM>::SetQuadrature(std::initializer_list<short> quad_points_per_axis)
{
    // Check if the desired quadrature points are set for all the axis.
    if (quad_points_per_axis.size() != 2) {
        throw std::invalid_argument(Logger::Error("Could not set quadrature for D3V8 element. The desired quadrature for all the 2 axis should be given."));
    }

    this->Quadrature().SetForQuadrilateralNat(quad_points_per_axis.begin()[0], quad_points_per_axis.begin()[1]);
}


template<short DIM>
void FemD2v4<DIM>::ComputeJacobians(const std::vector<IMP::Vec<DIM, double> > &nodes)
{
    // Check if shape function derivatives in natural coordinates have been computed for all the quadrature points.
    if (static_cast<int>(this->DerivsNatural().size()) != this->Quadrature().PointsNum()) {
        throw std::runtime_error(Logger::Error("Could not compute jacobians for D2V4 element. Compute first the shape function derivatives in the natural coordinates."));
    }

    // Reset jacobians to zero for all quadrature points.
    this->Jacobians().clear();
    this->Jacobians().resize(this->Quadrature().PointsNum(), Eigen::MatrixXd::Zero(2,2));

    // Reset inverse jacobians to zero for all quadrature points.
    this->InvJacobians().clear();
    this->InvJacobians().resize(this->Quadrature().PointsNum(), Eigen::MatrixXd::Zero(2,2));

    // Reset determinants of jacobians to zero for all quadrature points.
    this->DetJacobians().clear();
    this->DetJacobians().resize(this->Quadrature().PointsNum(), 0.);

    // Iterate over the element's quadrature points.
    for (int qid = 0; qid != this->Quadrature().PointsNum(); ++qid) {

        // Iterate over the element's nodes to compute the jacobian for the current quadrature point.
        auto nid = int{0};
        for (const auto &node : nodes) {
            this->Jacobians(qid).coeffRef(0,0) += this->DerivsNatural(qid).coeff(nid,0) * node[0];
            this->Jacobians(qid).coeffRef(0,1) += this->DerivsNatural(qid).coeff(nid,0) * node[1];

            this->Jacobians(qid).coeffRef(1,0) += this->DerivsNatural(qid).coeff(nid,1) * node[0];
            this->Jacobians(qid).coeffRef(1,1) += this->DerivsNatural(qid).coeff(nid,1) * node[1];
            nid++;
        }

        // Compute the inverse of the jacobian for the current quadrature point.
        this->InvJacobians(qid) = this->Jacobians(qid).inverse();

        // Compute the determinant of the jacobian for the current quadrature point.
        this->DetJacobians(qid) = this->Jacobians(qid).determinant();

        // Check for negative jacobian.
        if (this->DetJacobians(qid) < 0.) {
            std::string wrn_msg = "Negative jacobian occurred in bad shaped D2V4 element. Accuracy deterioration in the simulation is expected.";
            std::cout << Logger::Warning(wrn_msg) << "\n";
        }

    } // End of Iterate over the element's quadrature points.

}


template<short DIM>
void FemD2v4<DIM>::ComputeShapeFunctions()
{
    // Check if quadrature points have been set.
    if (this->Quadrature().PointsNum() == 0) {
        throw std::runtime_error(Logger::Error("Could not compute shape functions for D2V4 element. Set first the element's quadrature."));
    }

    // Reset the shape functions to zero for all quadrature points.
    this->ShapeFunctions().clear();
    this->ShapeFunctions().resize(this->Quadrature().PointsNum(), Eigen::VectorXd::Zero(4));

    // Compute shape functions for quadrature point.
    auto qid = int{0};
    for (const auto &quad_point : this->Quadrature().Points()) {
        this->ShapeFunctions(qid++) << 0.25*(1. - quad_point[0])*(1. - quad_point[1]),
                                       0.25*(1. + quad_point[0])*(1. - quad_point[1]),
                                       0.25*(1. + quad_point[0])*(1. + quad_point[1]),
                                       0.25*(1. - quad_point[0])*(1. + quad_point[1]);
    }

}


template<short DIM>
void FemD2v4<DIM>::ComputeDerivsNatural()
{
    // Check if quadrature points have been set.
    if (this->Quadrature().PointsNum() == 0) {
        throw std::runtime_error(Logger::Error("Could not compute derivatives in natural coordinates for D2V4 element. First must set the element's quadrature."));
    }

    // Reset natural derivatives to zero for all quadrature points.
    this->DerivsNatural().clear();
    this->DerivsNatural().resize(this->Quadrature().PointsNum(), Eigen::MatrixXd::Zero(4,2));

    // Compute the natural derivative for each quadrature point.
    auto qid = int{0};
    for (const auto &quad_point : this->Quadrature().Points()) {
        this->DerivsNatural(qid++) << -0.25*(1. - quad_point[1]), -0.25*(1. - quad_point[0]),
                                       0.25*(1. - quad_point[1]), -0.25*(1. + quad_point[0]),
                                       0.25*(1. + quad_point[1]),  0.25*(1. + quad_point[0]),
                                      -0.25*(1. + quad_point[1]),  0.25*(1. - quad_point[0]);

    }
}


template<short DIM>
void FemD2v4<DIM>::ComputeDerivs()
{
    // Check if jacobians have been computed for all the quadrature points.
    if (static_cast<int>(this->Jacobians().size()) != this->Quadrature().PointsNum()) {
        throw std::runtime_error(Logger::Error("Could not compute derivatives in physical coordinates for D2V4 element. First compute the element's jacobians."));
    }

    // Reset physical derivatives to zero for all quadrature points.
    this->Derivs().clear();
    this->Derivs().resize(this->Quadrature().PointsNum(), Eigen::MatrixXd::Zero(4,2));

    // Iterate over the element's quadrature points.
    for (int qid = 0; qid != this->Quadrature().PointsNum(); ++qid) {

        // Compute the first X derivative in physical coordinates.
        this->Derivs(qid).coeffRef(0,0) = this->InvJacobians(qid).coeff(0,0)*this->DerivsNatural(qid).coeff(0,0) +
                                          this->InvJacobians(qid).coeff(0,1)*this->DerivsNatural(qid).coeff(0,1);
        this->Derivs(qid).coeffRef(1,0) = this->InvJacobians(qid).coeff(0,0)*this->DerivsNatural(qid).coeff(1,0) +
                                          this->InvJacobians(qid).coeff(0,1)*this->DerivsNatural(qid).coeff(1,1);
        this->Derivs(qid).coeffRef(2,0) = this->InvJacobians(qid).coeff(0,0)*this->DerivsNatural(qid).coeff(2,0) +
                                          this->InvJacobians(qid).coeff(0,1)*this->DerivsNatural(qid).coeff(2,1);
        this->Derivs(qid).coeffRef(3,0) = this->InvJacobians(qid).coeff(0,0)*this->DerivsNatural(qid).coeff(3,0) +
                                          this->InvJacobians(qid).coeff(0,1)*this->DerivsNatural(qid).coeff(3,1);

        // Compute the first Y derivative in physical coordinates.
        this->Derivs(qid).coeffRef(0,1) = this->InvJacobians(qid).coeff(1,0)*this->DerivsNatural(qid).coeff(0,0) +
                                          this->InvJacobians(qid).coeff(1,1)*this->DerivsNatural(qid).coeff(0,1);
        this->Derivs(qid).coeffRef(1,1) = this->InvJacobians(qid).coeff(1,0)*this->DerivsNatural(qid).coeff(1,0) +
                                          this->InvJacobians(qid).coeff(1,1)*this->DerivsNatural(qid).coeff(1,1);
        this->Derivs(qid).coeffRef(2,1) = this->InvJacobians(qid).coeff(1,0)*this->DerivsNatural(qid).coeff(2,0) +
                                          this->InvJacobians(qid).coeff(1,1)*this->DerivsNatural(qid).coeff(2,1);
        this->Derivs(qid).coeffRef(3,1) = this->InvJacobians(qid).coeff(1,0)*this->DerivsNatural(qid).coeff(3,0) +
                                          this->InvJacobians(qid).coeff(1,1)*this->DerivsNatural(qid).coeff(3,1);
    }

}

} // End of namespace CLOUDEA

#endif //CLOUDEA_APPROXIMANTS_FEM_D2V4_TPP_