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


#ifndef ABSTRACTALGEBRAICCONTRACTIONMODEL_HPP_
#define ABSTRACTALGEBRAICCONTRACTIONMODEL_HPP_


#include "AbstractContractionModel.hpp"

/**
 *  Contraction models that give the active tension as an algebraic function
 *  of [Ca] or voltage, stretch, stretch rate and time; with no ODEs to
 *  integrate.
 */
class AbstractAlgebraicContractionModel : public AbstractContractionModel
{
protected:
    /** The time at the next timestep. Set in RunDoNotUpdate() */
    double mTime;

public:
    /** Constructor does nothing */
    AbstractAlgebraicContractionModel()
     : AbstractContractionModel()
    {
        mTime = 0.0;
    }

    /** No ODE to run, so this does nothing except save the time (using the
     *  time at the next timestep)
     *  @param startTime start time
     *  @param endTime end time
     *  @param timestep timestep for integrating ODEs if there were any
     */
    void RunDoNotUpdate(double startTime, double endTime, double timestep)
    {
        mTime = endTime;
    }

    /** No ODE to run, so this does nothing except save the time (using the
     *  time at the next timestep)
     *  @param startTime start time
     *  @param endTime end time
     *  @param timestep timestep for integrating ODEs if there were any
     */
    void RunAndUpdate(double startTime, double endTime, double timestep)
    {
        mTime = endTime;
    }

    /**
     *  @return same as GetActiveTension() for algebraic models (uses which stretch and
     *  and stretch rate has been passed in).
     */
    double GetNextActiveTension()
    {
        return GetActiveTension();
    }

    /** No ODE so does nothing. */
    void UpdateStateVariables()
    {
    }
};

#endif /*ABSTRACTALGEBRAICCONTRACTIONMODEL_HPP_*/
