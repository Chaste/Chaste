/*

Copyright (C) University of Oxford, 2005-2012

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

#ifdef CHASTE_CVODE
#ifndef _CVODEFIRSTORDER_HPP_
#define _CVODEFIRSTORDER_HPP_

#include "AbstractCvodeSystem.hpp"
#include "OdeSystemInformation.hpp"

/**
 * dy/dt = y,
 * so with y(0) = 1.0,
 * y = exp(t)
 * is the exact solution.
 */
class CvodeFirstOrder : public AbstractCvodeSystem
{
public:
    CvodeFirstOrder() : AbstractCvodeSystem(1) // 1 here is the number of variables
    {
        mpSystemInfo = OdeSystemInformation<CvodeFirstOrder>::Instance();
        Init();
    }

    void EvaluateYDerivatives(double time, const N_Vector rY, N_Vector rDY)
    {
        NV_Ith_S(rDY, 0) = NV_Ith_S(rY, 0);
    }
};

template<>
void OdeSystemInformation<CvodeFirstOrder>::Initialise()
{
    this->mVariableNames.push_back("Variable_1");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(1.0);

    this->mInitialised = true;
}

#endif //_CVODEFIRSTORDER_HPP_
#endif // CHASTE_CVODE

