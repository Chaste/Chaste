/*

Copyright (C) University of Oxford, 2005-2011

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

*/

#ifndef HEATEQUATIONWITHSOURCEFORCOUPLEDODESYSTEM_HPP_
#define HEATEQUATIONWITHSOURCEFORCOUPLEDODESYSTEM_HPP_

#include "AbstractLinearParabolicPdeSystemForCoupledOdeSystem.hpp"
#include "ChastePoint.hpp"
#include "Element.hpp"

/**
 * A PDE system defining the heat equation
 *
 * u_t = del^2 u + v,
 *
 * that couples with the single ODE
 *
 * dv/dt = a*v, v(0) = 1,
 *
 * which is defined in the separate class OdeSystemForCoupledHeatEquationWithSource.
 */
template<unsigned SPACE_DIM>
class HeatEquationWithSourceForCoupledOdeSystem : public AbstractLinearParabolicPdeSystemForCoupledOdeSystem<SPACE_DIM, SPACE_DIM, 1>
{
public:

    HeatEquationWithSourceForCoupledOdeSystem()
        : AbstractLinearParabolicPdeSystemForCoupledOdeSystem<SPACE_DIM, SPACE_DIM, 1>()
    {
    }

    double ComputeDuDtCoefficientFunction(const ChastePoint<SPACE_DIM>& rX, unsigned index)
    {
        return 1.0;
    }

    double ComputeSourceTerm(const ChastePoint<SPACE_DIM>& rX, c_vector<double,1>& rU, std::vector<double>& rOdeSolution, unsigned index)
    {
        return rOdeSolution[0];
    }

    c_matrix<double, SPACE_DIM, SPACE_DIM> ComputeDiffusionTerm(const ChastePoint<SPACE_DIM>& rX, unsigned index, Element<SPACE_DIM,SPACE_DIM>* pElement=NULL)
    {
        return identity_matrix<double>(SPACE_DIM);
    }
};

#endif /*HEATEQUATIONWITHSOURCEFORCOUPLEDODESYSTEM_HPP_*/
