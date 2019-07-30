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

#ifndef PDESIMULATIONTIME_HPP_
#define PDESIMULATIONTIME_HPP_

/**
 * \todo Integrate with other time stepping classes or remove?
 * A small convenience class providing a consistent global time to the
 * PDE solver classes.
 *
 * This isn't technically a singleton, as it's implemented with static
 * data and methods.
 */
class PdeSimulationTime
{
public:

    /**
     * Set the current time.
     *
     * @param time  the current time
     */
    static void SetTime(double time);

    /** @return the current time. */
    static double GetTime();

    /**
     * Set the current PDE timestep.
     *
     * The method checks that next_time ~= mTime + timestep
     * @param timestep  the current timestep
     * @param next_time  the next time (as given by the PDE time stepper).
     */
    static void SetPdeTimeStepAndNextTime(double timestep, double next_time);

    /** @return the current PDE timestep. */
    static double GetPdeTimeStep();

    /** @return the next time (time after time-step has been made).*/
    static double GetNextTime();

    /** @return 1/dt. */
    static double GetPdeTimeStepInverse();

private:

    /** The current time. */
    static double mTime;

    /** The timestep used in the PDE solve. */
    static double mPdeTimeStep;

    /** The next time (from the original PDE time-stepper).
     * mNextTime ~= mTime + mPdeTimeStep.
     * Note that this is stored explicitly because if we do the addition
     * then the answer will be off by mTime*DBL_EPSILON.
     */
    static double mNextTime;

    /** 1/dt. */
    static double mPdeTimeStepInverse;
};

#endif /*PDESIMULATIONTIME_HPP_*/
