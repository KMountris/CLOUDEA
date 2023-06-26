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



#ifndef CLOUDEA_APPROXIMANTS_FEM_D3V4_TPP_
#define CLOUDEA_APPROXIMANTS_FEM_D3V4_TPP_


#include "CLOUDEA/engine/approximants/fem_d3v4.hpp"


namespace CLOUDEA {

template<short DIM>
FemD3v4<DIM>::FemD3v4()
{}


template<short DIM>
FemD3v4<DIM>::~FemD3v4()
{}


template<short DIM>
void FemD3v4<DIM>::SetQuadrature(std::initializer_list<short> quad_points_num)
{
    this->Quadrature().SetForTetrahedronNat(quad_points_num.begin()[0]);
}


template<short DIM>
void FemD3v4<DIM>::ComputeJacobians(const std::vector<IMP::Vec<DIM, double> > &nodes)
{
    // Check if shape function derivatives in natural coordinates have been computed for all the quadrature points.
    if (static_cast<int>(this->DerivsNatural().size()) != this->Quadrature().PointsNum()) {
        throw std::runtime_error(Logger::Error("Could not compute jacobians for D3v4 element. Compute first the shape function derivatives in the natural coordinates."));
    }

    // Reset jacobians to zero for all quadrature points.
    this->Jacobians().clear();
    this->Jacobians().resize(this->Quadrature().PointsNum(), Eigen::MatrixXd::Zero(3,3));

    // Reset inverse jacobians to zero for all quadrature points.
    this->InvJacobians().clear();
    this->InvJacobians().resize(this->Quadrature().PointsNum(), Eigen::MatrixXd::Zero(3,3));

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
            this->Jacobians(qid).coeffRef(0,2) += this->DerivsNatural(qid).coeff(nid,0) * node[2];

            this->Jacobians(qid).coeffRef(1,0) += this->DerivsNatural(qid).coeff(nid,1) * node[0];
            this->Jacobians(qid).coeffRef(1,1) += this->DerivsNatural(qid).coeff(nid,1) * node[1];
            this->Jacobians(qid).coeffRef(1,2) += this->DerivsNatural(qid).coeff(nid,1) * node[2];

            this->Jacobians(qid).coeffRef(2,0) += this->DerivsNatural(qid).coeff(nid,2) * node[0];
            this->Jacobians(qid).coeffRef(2,1) += this->DerivsNatural(qid).coeff(nid,2) * node[1];
            this->Jacobians(qid).coeffRef(2,2) += this->DerivsNatural(qid).coeff(nid,2) * node[2];
            nid++;
        }

        // Compute the inverse of the jacobian for the current quadrature point.
        this->InvJacobians(qid) = this->Jacobians(qid).inverse();

        // Compute the determinant of the jacobian for the current quadrature point.
        this->DetJacobians(qid) = this->Jacobians(qid).determinant();

        // Check for negative jacobian.
        if (this->DetJacobians(qid) < 0.) {
            std::string wrn_msg = "Negative jacobian occurred in bad shaped D3v4 element. Accuracy deterioration in the simulation is expected.";
            std::cout << Logger::Warning(wrn_msg) << "\n";
        }

    } // End of Iterate over the element's quadrature points.

}


template<short DIM>
void FemD3v4<DIM>::ComputeShapeFunctions()
{
    // Check if quadrature points have been set.
    if (this->Quadrature().PointsNum() == 0) {
        throw std::runtime_error(Logger::Error("Could not compute shape functions for D3v4 element. Set first the element's quadrature."));
    }

    // Reset the shape functions to zero for all quadrature points.
    this->ShapeFunctions().clear();
    this->ShapeFunctions().resize(this->Quadrature().PointsNum(), Eigen::VectorXd::Zero(4));

    // Iterate over quadrature points.
    auto qid = int{0};
    for (const auto &quad_point : this->Quadrature().Points()) {
        // Shape function of quadrature point.
        this->ShapeFunctions(qid++) << 1. - quad_point[0] - quad_point[1] - quad_point[2],
                                       quad_point[0],
                                       quad_point[1],
                                       quad_point[2];
    }

}


template<short DIM>
void FemD3v4<DIM>::ComputeDerivsNatural()
{
    // Check if quadrature points have been set.
    if (this->Quadrature().PointsNum() == 0) {
        throw std::runtime_error(Logger::Error("Could not compute derivatives in natural coordinates for D3v4 element. First must set the element's quadrature."));
    }

    // Reset natural derivatives to zero for all quadrature points.
    this->DerivsNatural().clear();
    this->DerivsNatural().resize(this->Quadrature().PointsNum(), Eigen::MatrixXd::Zero(4,3));

    // Iterate over the element's quadrature points.
    for (int qid = 0; qid != this->Quadrature().PointsNum(); ++qid) {
        this->DerivsNatural(qid) << -1., -1., -1.,
                                     1.,  0.,  0.,
                                     0.,  1.,  0.,
                                     0.,  0.,  1.;
    } // End of Iterate over the element's quadrature points.
}


template<short DIM>
void FemD3v4<DIM>::ComputeDerivs()
{
    // Check if jacobians have been computed for all the quadrature points.
    if (static_cast<int>(this->Jacobians().size()) != this->Quadrature().PointsNum()) {
        throw std::runtime_error(Logger::Error("Could not compute derivatives in physical coordinates for D3v4 element. First compute the element's jacobians."));
    }

    // Reset physical derivatives to zero for all quadrature points.
    this->Derivs().clear();
    this->Derivs().resize(this->Quadrature().PointsNum(), Eigen::MatrixXd::Zero(4,3));

    // Iterate over the element's quadrature points.
    for (int qid = 0; qid != this->Quadrature().PointsNum(); ++qid) {

        // Compute the first X derivative in physical coordinates.
        this->Derivs(qid).coeffRef(0,0) = this->InvJacobians(qid).coeff(0,0)*this->DerivsNatural(qid).coeff(0,0) +
                                          this->InvJacobians(qid).coeff(0,1)*this->DerivsNatural(qid).coeff(0,1) +
                                          this->InvJacobians(qid).coeff(0,2)*this->DerivsNatural(qid).coeff(0,2);
        this->Derivs(qid).coeffRef(1,0) = this->InvJacobians(qid).coeff(0,0)*this->DerivsNatural(qid).coeff(1,0) +
                                          this->InvJacobians(qid).coeff(0,1)*this->DerivsNatural(qid).coeff(1,1) +
                                          this->InvJacobians(qid).coeff(0,2)*this->DerivsNatural(qid).coeff(1,2);
        this->Derivs(qid).coeffRef(2,0) = this->InvJacobians(qid).coeff(0,0)*this->DerivsNatural(qid).coeff(2,0) +
                                          this->InvJacobians(qid).coeff(0,1)*this->DerivsNatural(qid).coeff(2,1) +
                                          this->InvJacobians(qid).coeff(0,2)*this->DerivsNatural(qid).coeff(2,2);
        this->Derivs(qid).coeffRef(3,0) = this->InvJacobians(qid).coeff(0,0)*this->DerivsNatural(qid).coeff(3,0) +
                                          this->InvJacobians(qid).coeff(0,1)*this->DerivsNatural(qid).coeff(3,1) +
                                          this->InvJacobians(qid).coeff(0,2)*this->DerivsNatural(qid).coeff(3,2);

        // Compute the first Y derivative in physical coordinates.
        this->Derivs(qid).coeffRef(0,1) = this->InvJacobians(qid).coeff(1,0)*this->DerivsNatural(qid).coeff(0,0) +
                                          this->InvJacobians(qid).coeff(1,1)*this->DerivsNatural(qid).coeff(0,1) +
                                          this->InvJacobians(qid).coeff(1,2)*this->DerivsNatural(qid).coeff(0,2);
        this->Derivs(qid).coeffRef(1,1) = this->InvJacobians(qid).coeff(1,0)*this->DerivsNatural(qid).coeff(1,0) +
                                          this->InvJacobians(qid).coeff(1,1)*this->DerivsNatural(qid).coeff(1,1) +
                                          this->InvJacobians(qid).coeff(1,2)*this->DerivsNatural(qid).coeff(1,2);
        this->Derivs(qid).coeffRef(2,1) = this->InvJacobians(qid).coeff(1,0)*this->DerivsNatural(qid).coeff(2,0) +
                                          this->InvJacobians(qid).coeff(1,1)*this->DerivsNatural(qid).coeff(2,1) +
                                          this->InvJacobians(qid).coeff(1,2)*this->DerivsNatural(qid).coeff(2,2);
        this->Derivs(qid).coeffRef(3,1) = this->InvJacobians(qid).coeff(1,0)*this->DerivsNatural(qid).coeff(3,0) +
                                          this->InvJacobians(qid).coeff(1,1)*this->DerivsNatural(qid).coeff(3,1) +
                                          this->InvJacobians(qid).coeff(1,2)*this->DerivsNatural(qid).coeff(3,2);

        // Compute the first Z derivative in physical coordinates.
        this->Derivs(qid).coeffRef(0,2) = this->InvJacobians(qid).coeff(2,0)*this->DerivsNatural(qid).coeff(0,0) +
                                          this->InvJacobians(qid).coeff(2,1)*this->DerivsNatural(qid).coeff(0,1) +
                                          this->InvJacobians(qid).coeff(2,2)*this->DerivsNatural(qid).coeff(0,2);
        this->Derivs(qid).coeffRef(1,2) = this->InvJacobians(qid).coeff(2,0)*this->DerivsNatural(qid).coeff(1,0) +
                                          this->InvJacobians(qid).coeff(2,1)*this->DerivsNatural(qid).coeff(1,1) +
                                          this->InvJacobians(qid).coeff(2,2)*this->DerivsNatural(qid).coeff(1,2);
        this->Derivs(qid).coeffRef(2,2) = this->InvJacobians(qid).coeff(2,0)*this->DerivsNatural(qid).coeff(2,0) +
                                          this->InvJacobians(qid).coeff(2,1)*this->DerivsNatural(qid).coeff(2,1) +
                                          this->InvJacobians(qid).coeff(2,2)*this->DerivsNatural(qid).coeff(2,2);
        this->Derivs(qid).coeffRef(3,2) = this->InvJacobians(qid).coeff(2,0)*this->DerivsNatural(qid).coeff(3,0) +
                                          this->InvJacobians(qid).coeff(2,1)*this->DerivsNatural(qid).coeff(3,1) +
                                          this->InvJacobians(qid).coeff(2,2)*this->DerivsNatural(qid).coeff(3,2);
    }
}


} // End of namespace CLOUDEA


#endif //CLOUDEA_APPROXIMANTS_FEM_D3V4_TPP_