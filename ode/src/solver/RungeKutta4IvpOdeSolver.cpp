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

#include "RungeKutta4IvpOdeSolver.hpp"

void RungeKutta4IvpOdeSolver::CalculateNextYValue(AbstractOdeSystem* pAbstractOdeSystem,
                                                  double timeStep,
                                                  double time,
                                                  std::vector<double>& rCurrentYValues,
                                                  std::vector<double>& rNextYValues)
{
    /*
     * Apply Runge-Kutta 4th order method for each timestep in AbstractOneStepIvpSolver.
     * Calculates a vector containing the next Y value from the current one for each
     * equation in the system.
     */

    const unsigned num_equations = pAbstractOdeSystem->GetNumberOfStateVariables();

    if (num_equations != k1.size())
    {
        k1.resize(num_equations);
        k2.resize(num_equations);
        k3.resize(num_equations);
        k4.resize(num_equations);
        yki.resize(num_equations);
    }

    std::vector<double>& dy = rNextYValues; // re-use memory

    pAbstractOdeSystem->EvaluateYDerivatives(time, rCurrentYValues, dy);
    for (unsigned i=0; i<num_equations; i++)
    {
        k1[i] = timeStep*dy[i];
        yki[i] = rCurrentYValues[i] + 0.5*k1[i];
    }

    pAbstractOdeSystem->EvaluateYDerivatives(time+0.5*timeStep, yki, dy);
    for (unsigned i=0; i<num_equations; i++)
    {
        k2[i] = timeStep*dy[i];
        yki[i] = rCurrentYValues[i] + 0.5*k2[i];
    }

    pAbstractOdeSystem->EvaluateYDerivatives(time+0.5*timeStep, yki, dy);
    for (unsigned i=0; i<num_equations; i++)
    {
        k3[i] = timeStep*dy[i];
        yki[i] = rCurrentYValues[i] + k3[i];
    }

    pAbstractOdeSystem->EvaluateYDerivatives(time+timeStep, yki, dy);
    for (unsigned i=0; i<num_equations; i++)
    {
        k4[i] = timeStep*dy[i];
        rNextYValues[i] = rCurrentYValues[i] + (k1[i]+2*k2[i]+2*k3[i]+k4[i])/6.0;
    }
}


// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(RungeKutta4IvpOdeSolver)
