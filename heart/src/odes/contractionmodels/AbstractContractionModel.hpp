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

#ifndef ABSTRACTCONTRACTIONMODEL_HPP_
#define ABSTRACTCONTRACTIONMODEL_HPP_

#include <cassert>
#include "ContractionModelName.hpp"

/**
 *  Struct storing the input parameters that might be used by a contraction model (excl stretch and stretch-rate,
 *  as these may be set several times using the current deformation guess by the implicit assembler, and time).
 */
typedef struct ContractionModelInputParameters_
{
    double voltage;                           /**< Input voltage (mV)*/
    double intracellularCalciumConcentration; /**< Input calcium concentration (mMol) */
} ContractionModelInputParameters;


/**
 *  General interface for contraction models (models on the cell level determining
 *  the active tension (actually a stress) induced in a cell in response to
 *  electrical activity or mechanical deformation).
 */
class AbstractContractionModel
{
public:
    /**
     *  Constructor does nothing.
     */
    AbstractContractionModel()
    {
    }

    /**
     * Abstract destructor does nothing.
     */
    virtual ~AbstractContractionModel()
    {
    }

    /**
     *  @return whether the model depend on the stretch. (Pure, to be implemented in the concrete class).
     */
    virtual bool IsStretchDependent()=0;

    /**
     *  @return whether the model depend on the stretch-rate. (Pure, to be implemented in the concrete class).
     */
    virtual bool IsStretchRateDependent()=0;

    /**
     *  Set any input parameters (excl stretch and stretch rate). (Pure, to be implemented in the concrete class).
     *
     *  @param rInputParameters  contains various parameters: voltage, intracellular calcium concentration
     */
    virtual void SetInputParameters(ContractionModelInputParameters& rInputParameters)=0;

    /**
     *  Set the stretch and stretch rate. (Pure, to be implemented in the concrete class).
     *
     *  @param stretch  fibre stretch (dimensionless)
     *  @param stretchRate  fibre stretch rate (1/ms)
     */
    virtual void SetStretchAndStretchRate(double stretch, double stretchRate)=0;

    /** Safe setting of stretch-only, for stretch-rate independent models ONLY.
     *  @param stretch Stretch in fibre direction
     */
    void SetStretch(double stretch)
    {
        assert(!IsStretchRateDependent());
        SetStretchAndStretchRate(stretch, 0.0);
    }

    /**
     *  Run the contraction (ie if an ODE system) between the given times. This should NOT update any
     *  state variables. Call UpdateStateVariables() afterwards to update. For use in the implicit
     *  electromechanics algorithm
     *
     *  @param startTime start time
     *  @param endTime end time
     *  @param timeStep timestep to use in ODE solving.
     */
    virtual void RunDoNotUpdate(double startTime, double endTime, double timeStep)=0;

    /**
     *  Run the contraction (ie if an ODE system) between the given times, and update the state variables.
     *  For use in the explicit electromechanics algorithm
     *
     *  @param startTime start time
     *  @param endTime end time
     *  @param timeStep timestep to use in ODE solving.
     */
    virtual void RunAndUpdate(double startTime, double endTime, double timeStep)=0;


    /**
     *  After calling RunDoNotUpdate, which ran but should not have updated the state variables,
     *  the state variables are updated if this is called.
     */
    virtual void UpdateStateVariables()=0;

    /**
     *  @return the current active tension (note, actually a stress). (Pure, to be implemented in the concrete class).
     *
     *  DO NOT call inbetween RunDoNotUpdate() and UpdateStateVariables() as the old state variables but
     *  the next time would then be used in calculating Ta. Instead, use GetNextActiveTension(), or
     *  call UpdateStateVariables() and then this.
     */
    virtual double GetActiveTension()=0;

    /**
     *  @return the current active tension (note, actually a stress), using the current temporary
     *  state variables (ie those obtained after RunDoNotUpdate has been called).
     */
    virtual double GetNextActiveTension()=0;
};


#endif /*ABSTRACTCONTRACTIONMODEL_HPP_*/
