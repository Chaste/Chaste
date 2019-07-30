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

#ifndef ABSTRACTTIMEADAPTIVITYCONTROLLER_HPP_
#define ABSTRACTTIMEADAPTIVITYCONTROLLER_HPP_

#include "PetscVecTools.hpp"

/**
 *  Abstract class for defining rules for adaptive time-stepping.
 *
 */
class AbstractTimeAdaptivityController
{
private:

    /** Minimum timestep to be used */
    double mMinimumTimeStep;

    /** Maximum timestep to be used */
    double mMaximumTimeStep;

    /**
     * @return the timestep based on the
     * current solution and current time. Doesn't need to be checked to be
     * between the minimum and maximum timesteps.
     *
     * @param currentTime current time
     * @param currentSolution current solution
     */
    virtual double ComputeTimeStep(double currentTime, Vec currentSolution)=0;

public:

    /**
     * Constructor.
     *
     * @param minimumTimeStep minimum timestep to be used
     * @param maximumTimeStep maximum timestep to be used
     */
    AbstractTimeAdaptivityController(double minimumTimeStep, double maximumTimeStep)
        : mMinimumTimeStep(minimumTimeStep),
          mMaximumTimeStep(maximumTimeStep)
    {
        assert(minimumTimeStep>0.0);
        assert(maximumTimeStep>0.0);
        assert(minimumTimeStep < maximumTimeStep);
    }

    /** Destructor. */
    virtual ~AbstractTimeAdaptivityController()
    {
    }

    /**
     * @return the actual timestep to be used.
     *
     * @param currentTime current time
     * @param currentSolution current solution
     */
    double GetNextTimeStep(double currentTime, Vec currentSolution)
    {
        double dt = ComputeTimeStep(currentTime, currentSolution);
        if (dt < mMinimumTimeStep)
        {
            dt = mMinimumTimeStep;
        }
        if (dt > mMaximumTimeStep)
        {
            dt = mMaximumTimeStep;
        }
        return dt;
    }
};

#endif /*ABSTRACTTIMEADAPTIVITYCONTROLLER_HPP_*/
