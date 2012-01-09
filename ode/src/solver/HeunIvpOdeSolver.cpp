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

//Heun's Method
//Megan Lewis, University of Saskatchewan
//January 2009
#include "HeunIvpOdeSolver.hpp"

void HeunIvpOdeSolver::CalculateNextYValue(AbstractOdeSystem* pAbstractOdeSystem,
                                                  double timeStep,
                                                  double time,
                                                  std::vector<double>& rCurrentYValues,
                                                  std::vector<double>& rNextYValues)
{
    /*
     * Apply Heun 2nd order method for each time step in AbstractOneStepIvpSolver.
     * Calculates a vector containing the next Y value from the current one for each
     * equation in the system.
     */

    const unsigned num_equations = pAbstractOdeSystem->GetNumberOfStateVariables();

    std::vector<double> k1(num_equations);
    std::vector<double> k2(num_equations);
    std::vector<double>& dy = rNextYValues; // re-use memory

    // Work out k1
    pAbstractOdeSystem->EvaluateYDerivatives(time, rCurrentYValues, k1);

    // Add current y values to k1 values
    for (unsigned i=0; i<num_equations; i++)
    {
        dy[i] = timeStep*k1[i] + rCurrentYValues[i];
    }
    //Work out k2
    pAbstractOdeSystem->EvaluateYDerivatives(time+timeStep, dy, k2);

    // New solution
    for (unsigned i=0; i<num_equations; i++)
    {
        rNextYValues[i] = rCurrentYValues[i] + timeStep*0.5*(k1[i] + k2[i]);
    }
}


// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(HeunIvpOdeSolver)
