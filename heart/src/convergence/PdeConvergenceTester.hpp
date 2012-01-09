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


#ifndef PDECONVERGENCETESTER_HPP_
#define PDECONVERGENCETESTER_HPP_

#include "AbstractConvergenceTester.hpp"
/**
 * Drop the PDE time-step until a convergence criterion is met
 */
template<class CELL, class CARDIAC_PROBLEM, unsigned DIM, unsigned PROBLEM_DIM>
class PdeConvergenceTester : public AbstractConvergenceTester<CELL, CARDIAC_PROBLEM, DIM, PROBLEM_DIM>
{
public:
    /** The Pde time-step is setup.  Note that the Ode time-step is
     * set low, so that it won't need to be altered again.  We are assuming
     * that the convergence criterion will be met before PdeTimeStep<OdeTimeStep.
     */
    void SetInitialConvergenceParameters()
    {
        this->PdeTimeStep = 0.04;
        this->OdeTimeStep = 0.0025;
        //For extreme convergence tests this->OdeTimeStep = 7.8125e-5;
    }
    /**
     * Each new run has the #PdeTimeStep halved.
     */
    void UpdateConvergenceParameters()
    {
        this->PdeTimeStep *= 0.5;
    }
    /**
     * @return true to give up convergence when the PdeTimeStep>=OdeTimeStep requirement is violated.
     * This gives us the option to run with various (reasonable) OdeTimeStep parameters.
     */
    bool GiveUpConvergence()
    {
        return (this->PdeTimeStep<this->OdeTimeStep);
    }
    /**
     * @return the #PdeTimeStep as abcissa
     */
    double Abscissa()
    {
        return this->PdeTimeStep;
    }
};
#endif /*PDECONVERGENCETESTER_HPP_*/
