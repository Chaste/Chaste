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

#ifndef _ABSTRACTODESYSTEMFORCOUPLEDPDESYSTEM_HPP_
#define _ABSTRACTODESYSTEMFORCOUPLEDPDESYSTEM_HPP_

#include "AbstractOdeSystem.hpp"
#include "OdeSystemInformation.hpp"
#include "Exception.hpp"

/**
 * Represents an ODE system of the form
 *
 * d/dt (v_j) = g_j(x, u_1, ..., u_p, v_1, ..., v_q),  j=1,...,q,
 *
 * where the variables u_1, ..., u_p are assumed to satisfy a coupled PDE system of the form
 *
 * d/dt (u_i) = div (D(x) grad (u_i)) + f_i (x, u_1, ..., u_p, v_1, ..., v_q),  i=1,...,p.
 *
 * Such systems may be solved using LinearParabolicPdeSystemWithCoupledOdeSystemSolver.
 */
class AbstractOdeSystemForCoupledPdeSystem : public AbstractOdeSystem
{
protected:

    /**
     * Current solution to the PDE problem.
     */
    std::vector<double> mPdeSolution;

    /**
     * The size of the PDE solution at a point in space.
     */
    unsigned mPdeSolutionSize;

public:

    /**
     * Constructor.
     *
     * @param numberOfStateVariables  the number of state variables in the ODE system (defaults to 0)
     * @param pdeSolutionSize  the size of the PDE solution at a point in space (defaults to 0)
     */
    AbstractOdeSystemForCoupledPdeSystem(unsigned numberOfStateVariables = 0,
                                         unsigned pdeSolutionSize = 0)
        : AbstractOdeSystem(numberOfStateVariables),
          mPdeSolutionSize(pdeSolutionSize)
    {
        mPdeSolution.clear();
    }

    /**
     * Get #mPdeSolution.
     */
    std::vector<double>& rGetPdeSolution()
    {
        return mPdeSolution;
    }

    /**
     * Set #mPdeSolution.
     *
     * @param pdeSolution the PDE solution at a point in space
     */
    void SetPdeSolution(std::vector<double> pdeSolution)
    {
        if (pdeSolution.size() != mPdeSolutionSize)
        {
            EXCEPTION("The supplied vector is not the correct size.");
        }
        mPdeSolution = pdeSolution;
    }

    /**
     * Get #mPdeSolutionSize.
     */
    unsigned GetPdeSolutionSize()
    {
        return mPdeSolutionSize;
    }
};

#endif //_ABSTRACTODESYSTEMFORCOUPLEDPDESYSTEM_HPP_
