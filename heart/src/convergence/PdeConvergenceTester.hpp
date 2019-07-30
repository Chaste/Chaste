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
