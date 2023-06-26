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
   \file fem_mats.hpp
   \brief FemMats class header file.
   \author Konstantinos A. Mountris
   \date 28/02/2022
*/

#ifndef CLOUDEA_APPROXIMANTS_FEM_MATS_TPP_
#define CLOUDEA_APPROXIMANTS_FEM_MATS_TPP_


#include "CLOUDEA/engine/approximants/fem_mats.hpp"

namespace CLOUDEA {


template<short DIM, short CELL_VERTS>
FemMats<DIM,CELL_VERTS>::FemMats() : fbar_patches_(), shape_funcs_(), grads_(), jacobian_dets_(), gauss_weights_(), gauss_per_cell_(0)
{}


template<short DIM, short CELL_VERTS>
FemMats<DIM,CELL_VERTS>::~FemMats()
{}


template<short DIM, short CELL_VERTS>
void FemMats<DIM,CELL_VERTS>::ComputeFbarPatches(IMP::Mesh<DIM,CELL_VERTS> mesh)
{
    this->fbar_patches_.clear();
    this->fbar_patches_.resize(mesh.CellsNum());

    auto patch_size = 0;
    if (DIM == 2 && CELL_VERTS == 3) {
        patch_size = 2;
    } else if (DIM == 3 && CELL_VERTS == 4) {
        patch_size = 8;
    }

    auto in_patch_flag = std::vector<int>(mesh.CellsNum(), 0);
    auto attached_cell_ids = std::vector<int>{};
    for (auto i = int{0}; i != mesh.NodesNum(); ++i) {
        attached_cell_ids = mesh.AttachedCellIdsToNode(i);

        // Iterate attached cell ids to populate the patches.
        for (const auto &cid : attached_cell_ids) {
            for (const auto &neigh : attached_cell_ids) {
                if (this->fbar_patches_[cid].size() != patch_size) {
                    this->fbar_patches_[cid].emplace_back(neigh);
                }
            }
            in_patch_flag[cid] = 1;
        }
    }

    for (const auto &patch : this->fbar_patches_) {
        std::cout << patch.size() << std::endl;
    }


}


template<short DIM, short CELL_VERTS>
void FemMats<DIM,CELL_VERTS>::ComputeMatrices(const IMP::Mesh<DIM,CELL_VERTS> &mesh)
{
    // Create a finite element of the given type.
    auto elem = CLOUDEA::FemFactory<DIM>::Create(mesh.CellsShape());

    // Set the quadrature of the finite element.
    switch (mesh.CellsShape()) {
        case IMP::CellShape::tri :
            elem->SetQuadrature({3});
            this->gauss_per_cell_ = 3;
        break;
        case IMP::CellShape::quad :
            elem->SetQuadrature({2,2});
            this->gauss_per_cell_ = 4;
        break;
        case IMP::CellShape::tet :
            elem->SetQuadrature({1});
            this->gauss_per_cell_ = 1;
        break;
        case IMP::CellShape::hex :
            elem->SetQuadrature({2,2,2});
            this->gauss_per_cell_ = 8;
        break;
        default:
            auto err_msg = "Could not compute FEM Matrices. It is required a mesh with cells of type: [tri | quad | tet | hex].";
            throw std::invalid_argument(Logger::Error(err_msg));
        break;
    }

    // Compute natural derivatives and shape functions of the element.
    elem->ComputeDerivsNatural();
    elem->ComputeShapeFunctions();

    // Initialize FEM matrices containers.
    auto n = mesh.CellsNum()*elem->Quadrature().PointsNum();
    this->shape_funcs_.clear();  this->shape_funcs_.reserve(n);
    this->grads_.clear();      this->grads_.reserve(n);
    this->jacobian_dets_.clear();  this->jacobian_dets_.reserve(n);
    this->gauss_weights_.clear();  this->gauss_weights_.reserve(n);

    // Compute shape functions and derivatives.
    auto nodes = std::vector<IMP::Vec<DIM,double>>{CELL_VERTS};
    for (int cid = 0; cid != mesh.CellsNum(); ++cid) {

        // Get node coordinates of the current cell.
        for (int i = 0; i != CELL_VERTS; ++i) {
            for (short d = 0; d != DIM; ++d) { nodes[i][d] = mesh.Nodes(mesh.Cells(cid).N(i))[d]; }
        }

        // Compute the jacobian and physical derivatives for the current cell.
        elem->ComputeJacobians(nodes);
        elem->ComputeDerivs();

        // Store in FEM matrices containers.
        this->shape_funcs_.insert(this->shape_funcs_.end(), elem->ShapeFunctions().begin(), elem->ShapeFunctions().end());
        this->grads_.insert(this->grads_.end(), elem->Derivs().begin(), elem->Derivs().end());
        this->jacobian_dets_.insert(this->jacobian_dets_.end(), elem->DetJacobians().begin(), elem->DetJacobians().end());
        this->gauss_weights_.insert(this->gauss_weights_.end(), elem->Quadrature().Weights().begin(), elem->Quadrature().Weights().end());
    } // End of Iterate over the tissue mesh cells.
}


} // End of namespace CLOUDEA

#endif // CLOUDEA_APPROXIMANTS_FEM_MATS_TPP_