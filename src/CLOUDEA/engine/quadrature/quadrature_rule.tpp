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

#ifndef CLOUDEA_QUADRATURE_QUADRATURE_RULE_TPP_
#define CLOUDEA_QUADRATURE_QUADRATURE_RULE_TPP_

#include "CLOUDEA/engine/quadrature/quadrature_rule.hpp"


namespace CLOUDEA {

template<short DIM>
QuadratureRule<DIM>::QuadratureRule() : points_(), weights_(), order_(0)
{}


template<short DIM>
QuadratureRule<DIM>::~QuadratureRule()
{}


template<short DIM>
void QuadratureRule<DIM>::SetForLineNat(short quad_points_num)
{
    // Check if dimensions are at least one.
    if (DIM < 1) {
        throw std::invalid_argument(Logger::Error("Could not generate quadrature for line in natural coordinate system. Quadrature dimensions (DIM) must be at least 1."));
    }

    // Reset the quadrature points and weights containers to zero.
    this->points_.clear();  this->points_.resize(quad_points_num, IMP::Vec<DIM, double>());
    this->weights_.clear(); this->weights_.resize(quad_points_num, 0.);

    // Assign the requested corresponding quadrature points and weights.
    switch (quad_points_num) {
        case 1:
            this->order_ = 1;

            this->points_[0].Set({0.});

            this->weights_[0] = 2.;
        break;

        case 2:
            this->order_ = 2;

            this->points_[0].Set({-0.577350269189625764509149});
            this->points_[1].Set({0.577350269189625764509149});

            this->weights_[0] = 1.;
            this->weights_[1] = 1.;
        break;

        case 3:
            this->order_ = 3;

            this->points_[0].Set({-0.774596669241483377035835});
            this->points_[1].Set({0.});
            this->points_[2].Set({0.774596669241483377035835});

            this->weights_[0] = 0.555555555555555555555556;
            this->weights_[1] = 0.888888888888888888888889;
            this->weights_[2] = 0.555555555555555555555556;
        break;

        case 4:
            this->order_ = 4;

            this->points_[0].Set({-0.861136311594052575223946});
            this->points_[1].Set({-0.339981043584856264802666});
            this->points_[2].Set({0.339981043584856264802666});
            this->points_[3].Set({0.861136311594052575223946});

            this->weights_[0] = 0.34785484513745385737306;
            this->weights_[1] = 0.65214515486254614262694;
            this->weights_[2] = 0.65214515486254614262694;
            this->weights_[3] = 0.34785484513745385737306;
        break;

        case 5:
            this->order_ = 5;

            this->points_[0].Set({-0.906179845938663992797627});
            this->points_[1].Set({-0.538469310105683091036314});
            this->points_[2].Set({0.});
            this->points_[3].Set({0.538469310105683091036314});
            this->points_[4].Set({0.906179845938663992797627});

            this->weights_[0] = 0.23692688505618908751426;
            this->weights_[1] = 0.47862867049936646804129;
            this->weights_[2] = 0.56888888888888888888889;
            this->weights_[3] = 0.47862867049936646804129;
            this->weights_[4] = 0.23692688505618908751426;
        break;

        case 6:
            this->order_ = 6;

            this->points_[0].Set({-0.932469514203152027812302});
            this->points_[1].Set({-0.6612093864662645136614});
            this->points_[2].Set({-0.238619186083196908630502});
            this->points_[3].Set({0.238619186083196908630502});
            this->points_[4].Set({0.6612093864662645136614});
            this->points_[5].Set({0.932469514203152027812302});

            this->weights_[0] = 0.1713244923791703450403;
            this->weights_[1] = 0.36076157304813860756983;
            this->weights_[2] = 0.46791393457269104738987;
            this->weights_[3] = 0.46791393457269104738987;
            this->weights_[4] = 0.36076157304813860756983;
            this->weights_[5] = 0.1713244923791703450403;
        break;

        case 7:
            this->order_ = 7;

            this->points_[0].Set({-0.94910791234275852452619});
            this->points_[1].Set({-0.741531185599394439863865});
            this->points_[2].Set({-0.405845151377397166906607});
            this->points_[3].Set({0.});
            this->points_[4].Set({0.405845151377397166906607});
            this->points_[5].Set({0.741531185599394439863865});
            this->points_[6].Set({0.94910791234275852452619});

            this->weights_[0] = 0.12948496616886969327061;
            this->weights_[1] = 0.27970539148927666790147;
            this->weights_[2] = 0.38183005050511894495037;
            this->weights_[3] = 0.4179591836734693877551;
            this->weights_[4] = 0.38183005050511894495037;
            this->weights_[5] = 0.27970539148927666790147;
            this->weights_[6] = 0.12948496616886969327061;
        break;

        case 8:
            this->order_ = 8;

            this->points_[0].Set({-0.960289856497536231683561});
            this->points_[1].Set({-0.796666477413626739591554});
            this->points_[2].Set({-0.525532409916328985817739});
            this->points_[3].Set({-0.183434642495649804939476});
            this->points_[4].Set({0.183434642495649804939476});
            this->points_[5].Set({0.525532409916328985817739});
            this->points_[6].Set({0.796666477413626739591554});
            this->points_[7].Set({0.960289856497536231683561});

            this->weights_[0] = 0.10122853629037625915253;
            this->weights_[1] = 0.22238103445337447054436;
            this->weights_[2] = 0.31370664587788728733796;
            this->weights_[3] = 0.36268378337836198296515;
            this->weights_[4] = 0.36268378337836198296515;
            this->weights_[5] = 0.31370664587788728733796;
            this->weights_[6] = 0.22238103445337447054436;
            this->weights_[7] = 0.10122853629037625915253;
        break;

        case 9:
            this->order_ = 9;

            this->points_[0].Set({-0.968160239507626089835576});
            this->points_[1].Set({-0.83603110732663579429943});
            this->points_[2].Set({-0.613371432700590397308702});
            this->points_[3].Set({-0.324253423403808929038538});
            this->points_[4].Set({0.});
            this->points_[5].Set({0.324253423403808929038538});
            this->points_[6].Set({0.613371432700590397308702});
            this->points_[7].Set({0.83603110732663579429943});
            this->points_[8].Set({0.968160239507626089835576});

            this->weights_[0] = 0.08127438836157441197189;
            this->weights_[1] = 0.18064816069485740405847;
            this->weights_[2] = 0.26061069640293546231874;
            this->weights_[3] = 0.31234707704000284006863;
            this->weights_[4] = 0.33023935500125976316453;
            this->weights_[5] = 0.31234707704000284006863;
            this->weights_[6] = 0.26061069640293546231874;
            this->weights_[7] = 0.18064816069485740405847;
            this->weights_[8] = 0.08127438836157441197189;
        break;

        case 10:
            this->order_ = 10;

            this->points_[0].Set({-0.973906528517171720077964});
            this->points_[1].Set({-0.865063366688984510732097});
            this->points_[2].Set({-0.679409568299024406234327});
            this->points_[3].Set({-0.433395394129247190799266});
            this->points_[4].Set({-0.148874338981631210884826});
            this->points_[5].Set({0.148874338981631210884826});
            this->points_[6].Set({0.433395394129247190799266});
            this->points_[7].Set({0.865063366688984510732097});
            this->points_[8].Set({0.679409568299024406234327});
            this->points_[9].Set({0.973906528517171720077964});

            this->weights_[0] = 0.06667134430868813759357;
            this->weights_[1] = 0.14945134915058059314578;
            this->weights_[2] = 0.21908636251598204399554;
            this->weights_[3] = 0.26926671930999635509123;
            this->weights_[4] = 0.29552422471475287017389;
            this->weights_[5] = 0.29552422471475287017389;
            this->weights_[6] = 0.26926671930999635509123;
            this->weights_[7] = 0.14945134915058059314578;
            this->weights_[8] = 0.21908636251598204399554;
            this->weights_[9] = 0.06667134430868813759357;
        break;

        case 11:
            this->order_ = 11;

            this->points_[0].Set({-0.987228658146056992803938});
            this->points_[1].Set({-0.887062599768095299075158});
            this->points_[2].Set({-0.730152005574049324093416});
            this->points_[3].Set({-0.519096129206811815925726});
            this->points_[4].Set({-0.269543155952344972331532});
            this->points_[5].Set({0.});
            this->points_[6].Set({0.269543155952344972331532});
            this->points_[7].Set({0.519096129206811815925726});
            this->points_[8].Set({0.730152005574049324093416});
            this->points_[9].Set({0.887062599768095299075158});
            this->points_[10].Set({0.987228658146056992803938});

            this->weights_[0] = 0.05566856711617366648275;
            this->weights_[1] = 0.12558036946490462463469;
            this->weights_[2] = 0.1862902109277342514261;
            this->weights_[3] = 0.23319376459199047991852;
            this->weights_[4] = 0.26280454451024666218069;
            this->weights_[5] = 0.27292508677790063071448;
            this->weights_[6] = 0.26280454451024666218069;
            this->weights_[7] = 0.23319376459199047991852;
            this->weights_[8] = 0.1862902109277342514261;
            this->weights_[9] = 0.12558036946490462463469;
            this->weights_[10] = 0.05566856711617366648275;
        break;

        case 12:
            this->order_ = 12;

            this->points_[0].Set({-0.981560634246719250690549});
            this->points_[1].Set({-0.904117256370474856678466});
            this->points_[2].Set({-0.769002674194304687036894});
            this->points_[3].Set({-0.587317954286617447296702});
            this->points_[4].Set({-0.367831498998180193752692});
            this->points_[5].Set({-0.125233408511468915472441});
            this->points_[6].Set({0.125233408511468915472441});
            this->points_[7].Set({0.367831498998180193752692});
            this->points_[8].Set({0.587317954286617447296702});
            this->points_[9].Set({0.769002674194304687036894});
            this->points_[10].Set({0.904117256370474856678466});
            this->points_[11].Set({0.981560634246719250690549});

            this->weights_[0] = 0.04717533638651182719462;
            this->weights_[1] = 0.10693932599531843096025;
            this->weights_[2] = 0.16007832854334622633465;
            this->weights_[3] = 0.20316742672306592174906;
            this->weights_[4] = 0.23349253653835480876085;
            this->weights_[5] = 0.24914704581340278500056;
            this->weights_[6] = 0.24914704581340278500056;
            this->weights_[7] = 0.23349253653835480876085;
            this->weights_[8] = 0.20316742672306592174906;
            this->weights_[9] = 0.16007832854334622633465;
            this->weights_[10] = 0.10693932599531843096025;
            this->weights_[11] = 0.04717533638651182719462;
        break;

        default:
            std::string error_msg = "Could not generate quadrature for line. Unsupported quadrature points number requested. Supported: 1 up to 12";
            throw std::invalid_argument(Logger::Error(error_msg));
        break;
    }

}


template<short DIM>
void QuadratureRule<DIM>::SetForQuadrilateralNat(short quad_points_num_x, short quad_points_num_y)
{
    // Check if dimensions are at least two.
    if (DIM < 2) {
        throw std::invalid_argument(Logger::Error("Could not generate quadrature for quadrilateral in natural coordinate system. Quadrature dimensions (DIM) must be at least 2."));
    }

    // Reset the quadrature points and weights containers to zero.
    this->points_.clear();  this->points_.resize(quad_points_num_x*quad_points_num_y, IMP::Vec<DIM, double>());
    this->weights_.clear(); this->weights_.resize(quad_points_num_x*quad_points_num_y, 0.);

    // Set as order of the element the minimum per axis.
    this->order_ = std::min(quad_points_num_x, quad_points_num_y);

    // Get the quadrature over X axis.
    QuadratureRule<1> quadrature_x;
    quadrature_x.SetForLineNat(quad_points_num_x);

    // Get the quadrature over Y axis.
    QuadratureRule<1> quadrature_y;
    quadrature_y.SetForLineNat(quad_points_num_y);

    // Set the quadrature for the quadrilateral element.
    for (int i = 0; i != quad_points_num_x; ++i) {
        for (int j = 0; j != quad_points_num_y; ++j) {
            // Get the index of the current quadrature point.
            int id = i*quad_points_num_y + j;

            // Set the coordianates of the quadrature point with index id.
            this->points_[id].Set({quadrature_x.Points()[i][0], quadrature_y.Points()[j][0]});

            // Set the weight of the quadrature point with index id.
            this->weights_[id] = quadrature_x.Weights()[i] * quadrature_y.Weights()[j];
        }
    }

}


template<short DIM>
void QuadratureRule<DIM>::SetForTriangleNat(short quad_points_num)
{
    // Check if dimensions are at least two.
    if (DIM < 2) {
        throw std::invalid_argument(Logger::Error("Could not generate quadrature for triangle in natural coordinate system. Quadrature dimensions (DIM) must be at least 2."));
    }

    // Reset the quadrature points and weights containers to zero.
    this->points_.clear();  this->points_.resize(quad_points_num, IMP::Vec<DIM, double>());
    this->weights_.clear(); this->weights_.resize(quad_points_num, 0.);

    // Assign the requested corresponding quadrature points and weights.
    switch (quad_points_num) {
        case 1:
            this->order_ = 1;

            this->points_[0].Set({0.333333333333333333, 0.333333333333333333});

            this->weights_[0] = 0.5;
        break;

        case 3:
            this->order_ = 2;

            this->points_[0].Set({0.666666666666666667, 0.166666666666666667});
            this->points_[1].Set({0.166666666666666667, 0.666666666666666667});
            this->points_[2].Set({0.166666666666666667, 0.166666666666666667});

            this->weights_[0] = 0.166666666666666666;
            this->weights_[1] = 0.166666666666666667;
            this->weights_[2] = 0.166666666666666667;
        break;

        case 4:
            this->order_ = 3;
            this->points_[0].Set({0.333333333333333333, 0.333333333333333333});
            this->points_[1].Set({0.6, 0.2});
            this->points_[2].Set({0.2, 0.6});
            this->points_[3].Set({0.2, 0.2});

            this->weights_[0] = -0.28125;
            this->weights_[1] = 0.260416666666666667;
            this->weights_[2] = 0.260416666666666667;
            this->weights_[3] = 0.260416666666666666;
        break;

        case 6:
            this->order_ = 4;

            this->points_[0].Set({0.10810301816807, 0.4459484909159650});
            this->points_[1].Set({0.445948490915965, 0.10810301816807});
            this->points_[2].Set({0.445948490915965, 0.4459484909159650});
            this->points_[3].Set({0.816847572980459, 0.091576213509771});
            this->points_[4].Set({0.091576213509771, 0.816847572980459});
            this->points_[5].Set({0.091576213509771, 0.091576213509771});

            this->weights_[0] = 0.1116907948390055;
            this->weights_[1] = 0.1116907948390055;
            this->weights_[2] = 0.1116907948390055;
            this->weights_[3] = 0.054975871827661;
            this->weights_[4] = 0.054975871827661;
            this->weights_[5] = 0.054975871827661;

        break;

        case 7:
            this->order_ = 5;

            this->points_[0].Set({0.333333333333333, 0.333333333333333});
            this->points_[1].Set({0.05971587178977, 0.470142064105115});
            this->points_[2].Set({0.470142064105115, 0.05971587178977});
            this->points_[3].Set({0.470142064105115, 0.470142064105115});
            this->points_[4].Set({0.797426985353087, 0.101286507323456});
            this->points_[5].Set({0.101286507323456, 0.797426985353087});
            this->points_[6].Set({0.101286507323456, 0.101286507323456});

            this->weights_[0] = 0.1125;
            this->weights_[1] = 0.066197076394253;
            this->weights_[2] = 0.066197076394253;
            this->weights_[3] = 0.066197076394253;
            this->weights_[4] = 0.0629695902724135;
            this->weights_[5] = 0.0629695902724135;
            this->weights_[6] = 0.0629695902724135;
        break;

        case 12:
            this->order_ = 6;

            this->points_[0].Set({0.501426509658179, 0.24928674517091});
            this->points_[1].Set({0.24928674517091, 0.501426509658179});
            this->points_[2].Set({0.24928674517091, 0.24928674517091});
            this->points_[3].Set({0.873821971016996, 0.063089014491502});
            this->points_[4].Set({0.063089014491502, 0.873821971016996});
            this->points_[5].Set({0.063089014491502, 0.063089014491502});
            this->points_[6].Set({0.053145049844817, 0.310352451033784});
            this->points_[7].Set({0.636502499121399, 0.053145049844817});
            this->points_[8].Set({0.310352451033784, 0.636502499121399});
            this->points_[9].Set({0.053145049844817, 0.636502499121399});
            this->points_[10].Set({0.636502499121399, 0.310352451033784});
            this->points_[11].Set({0.310352451033784, 0.053145049844817});

            this->weights_[0] = 0.0583931378631895;
            this->weights_[1] = 0.0583931378631895;
            this->weights_[2] = 0.0583931378631895;
            this->weights_[3] = 0.0254224531851035;
            this->weights_[4] = 0.0254224531851035;
            this->weights_[5] = 0.0254224531851035;
            this->weights_[6] = 0.041425537809187;
            this->weights_[7] = 0.041425537809187;
            this->weights_[8] = 0.041425537809187;
            this->weights_[9] = 0.041425537809187;
            this->weights_[10] = 0.041425537809187;
            this->weights_[11] = 0.041425537809187;
        break;

        default:
            std::string error_msg = "Could not generate quadrature for triangle. Unsupported quadrature points number requested. "
                                     "Supported: 1 | 3 | 4 | 6 | 7 | 12 | 13 | 16 | 19 | 25 | 27 | 33 | 37 | 42 | 48 | 52 | 61";
            throw std::invalid_argument(Logger::Error(error_msg));
        break;
    }
}


template<short DIM>
void QuadratureRule<DIM>::SetForTetrahedronNat(short quad_points_num)
{
    // Check if dimensions are at least three.
    if (DIM < 3) {
        throw std::invalid_argument(Logger::Error("Could not generate quadrature for tetrahedron in natural coordinate system. Quadrature dimensions (DIM) must be at least 3."));
    }

    // Reset the quadrature points and weights containers to zero.
    this->points_.clear();  this->points_.resize(quad_points_num, IMP::Vec<DIM, double>());
    this->weights_.clear(); this->weights_.resize(quad_points_num, 0.);

    // Assign the requested corresponding quadrature points and weights.
    switch (quad_points_num) {
        case 1:
            this->order_ = 0;

            this->points_[0].Set({0.25, 0.25, 0.25});

            this->weights_[0] = 0.166666666666666666;
        break;

        case 4:
            this->order_ = 1;

            this->points_[0].Set({0.585410196624968515, 0.138196601125010504, 0.138196601125010504});
            this->points_[1].Set({0.138196601125010504, 0.585410196624968515, 0.138196601125010504});
            this->points_[2].Set({0.138196601125010504, 0.138196601125010504, 0.585410196624968515});
            this->points_[3].Set({0.138196601125010504, 0.138196601125010504, 0.138196601125010504});

            this->weights_[0] = 0.0416666666666666667;
            this->weights_[1] = 0.0416666666666666667;
            this->weights_[2] = 0.0416666666666666667;
            this->weights_[3] = 0.0416666666666666667;
        break;

        case 5:
            this->order_ = 2;

            this->points_[0].Set({0.25, 0.25, 0.25});
            this->points_[1].Set({0.5, 0.166666666666666667, 0.166666666666666667});
            this->points_[2].Set({0.166666666666666667, 0.5, 0.166666666666666667});
            this->points_[3].Set({0.166666666666666667, 0.166666666666666667, 0.5});
            this->points_[4].Set({0.166666666666666667, 0.166666666666666667, 0.166666666666666667});

            this->weights_[0] = -0.133333333333333333;
            this->weights_[1] = 0.0749999999999999999;
            this->weights_[2] = 0.0749999999999999999;
            this->weights_[3] = 0.0749999999999999999;
            this->weights_[4] = 0.0749999999999999999;
        break;

        case 10:
            this->order_ = 3;

            this->points_[0].Set({0.568430584196844446, 0.143856471934385194, 0.143856471934385194});
            this->points_[1].Set({0.143856471934385194, 0.568430584196844446, 0.143856471934385194});
            this->points_[2].Set({0.143856471934385194, 0.143856471934385194, 0.568430584196844446});
            this->points_[3].Set({0.143856471934385194, 0.143856471934385194, 0.143856471934385194});
            this->points_[4].Set({0.5, 0.5, 0.});
            this->points_[5].Set({0.5, 0., 0.5});
            this->points_[6].Set({0.5, 0., 0.});
            this->points_[7].Set({0., 0.5, 0.5});
            this->points_[8].Set({0., 0.5, 0.});
            this->points_[9].Set({0., 0., 0.5});

            this->weights_[0] = 0.0362941783134008988;
            this->weights_[1] = 0.0362941783134008988;
            this->weights_[2] = 0.0362941783134008988;
            this->weights_[3] = 0.0362941783134008988;
            this->weights_[4] = 0.00358165890217718337;
            this->weights_[5] = 0.00358165890217718337;
            this->weights_[6] = 0.00358165890217718337;
            this->weights_[7] = 0.00358165890217718337;
            this->weights_[8] = 0.00358165890217718337;
            this->weights_[9] = 0.00358165890217718337;
        break;

        default:
            std::string error_msg = "Could not generate quadrature for tetrahedron. Unsupported quadrature points number requested. "
                                     "Supported: 1 | 4 | 5 | 10 | 11 | 15 | 24 | 31 | 45";
            throw std::invalid_argument(Logger::Error(error_msg));
        break;
    }

}


template<short DIM>
void QuadratureRule<DIM>::SetForHexahedronNat(short quad_points_num_x, short quad_points_num_y, short quad_points_num_z)
{
    // Check if dimensions are at least two.
    if (DIM < 3) {
        throw std::invalid_argument(Logger::Error("Could not generate quadrature for hexahedron in natural coordinate system. Quadrature dimensions (DIM) must be at least 2."));
    }

    // Reset the quadrature points and weights containers to zero.
    this->points_.clear();  this->points_.resize(quad_points_num_x*quad_points_num_y*quad_points_num_z, IMP::Vec<DIM, double>());
    this->weights_.clear(); this->weights_.resize(quad_points_num_x*quad_points_num_y*quad_points_num_z, 0.);

    // Set as order of the element the minimum per axis.
    this->order_ = std::min({quad_points_num_x, quad_points_num_y, quad_points_num_z});

    // Get the quadrature over X axis.
    QuadratureRule<1> quadrature_x;
    quadrature_x.SetForLineNat(quad_points_num_x);

    // Get the quadrature over Y axis.
    QuadratureRule<1> quadrature_y;
    quadrature_y.SetForLineNat(quad_points_num_y);

    // Get the quadrature over Z axis.
    QuadratureRule<1> quadrature_z;
    quadrature_z.SetForLineNat(quad_points_num_z);

    // Set the quadrature for the hexahedral element.
    for (int i = 0; i != quad_points_num_x; ++i) {
        for (int j = 0; j != quad_points_num_y; ++j) {
            for (int k = 0; k != quad_points_num_z; ++k) {
                // Get the index of the current quadrature point.
                int id = i*quad_points_num_y + j*quad_points_num_z + k;

                // Set the coordianates of the quadrature point with index id.
                this->points_[id].Set({quadrature_x.Points()[i][0], quadrature_y.Points()[j][0], quadrature_z.Points()[i][0]});

                // Set the weight of the quadrature point with index id.
                this->weights_[id] = quadrature_x.Weights()[i] * quadrature_y.Weights()[j] * quadrature_z.Weights()[k];
            }
        }
    }

}


} // End of namespace CLOUDEA


#endif //CLOUDEA_QUADRATURE_QUADRATURE_RULE_TPP_