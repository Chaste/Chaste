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


#ifndef _ABSTRACTODESYSTEMWITHANALYTICJACOBIAN_HPP_
#define _ABSTRACTODESYSTEMWITHANALYTICJACOBIAN_HPP_

#include "Exception.hpp"
#include "AbstractOdeSystem.hpp"

/**
 * Abstract analytic Jacobian class.
 *
 * Represents an ODE system with an analytic Jacobian available,
 * which can be computed using the method AnalyticJacobian.
 */
class AbstractOdeSystemWithAnalyticJacobian : public AbstractOdeSystem
{
public:

    /**
     * Constructor.
     *
     * @param numberOfStateVariables  the number of state variables in the ODE system (defaults to 0)
     */
    AbstractOdeSystemWithAnalyticJacobian(unsigned numberOfStateVariables = 0)
        : AbstractOdeSystem(numberOfStateVariables)
    {
        mUseAnalyticJacobian = true;
    }

    /**
     * Compute the analytic Jacobian matrix of the ODE system.
     *
     * @param rSolutionGuess  the current guess at the solution for this time step
     * @param jacobian  will be filled in with the Jacobian matrix entries
     * @param time  the current simulation time
     * @param timeStep  the time step in use by the integrator at present
     */
    virtual void AnalyticJacobian(const std::vector<double>& rSolutionGuess,
                                  double** jacobian,
                                  double time,
                                  double timeStep) = 0;

};

#endif //_ABSTRACTODESYSTEMWITHANALYTICJACOBIAN_HPP_
