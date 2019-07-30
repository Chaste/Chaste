/*

Copyright (c) 2005-2019, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

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
     * @return #mPdeSolution.
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
     * @return #mPdeSolutionSize.
     */
    unsigned GetPdeSolutionSize()
    {
        return mPdeSolutionSize;
    }
};

#endif //_ABSTRACTODESYSTEMFORCOUPLEDPDESYSTEM_HPP_
