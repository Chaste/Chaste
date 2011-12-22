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
     *  Does the model depend on the stretch. (Pure, to be implemented in the concrete class).
     */
    virtual bool IsStretchDependent()=0;

    /**
     *  Does the model depend on the stretch-rate. (Pure, to be implemented in the concrete class).
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
     *  Get the current active tension (note, actually a stress). (Pure, to be implemented in the concrete class).
     *
     *  DO NOT call inbetween RunDoNotUpdate() and UpdateStateVariables() as the old state variables but
     *  the next time would then be used in calculating Ta. Instead, use GetNextActiveTension(), or
     *  call UpdateStateVariables() and then this.
     */
    virtual double GetActiveTension()=0;

    /**
     *  Get the current active tension (note, actually a stress), using the current temporary
     *  state variables (ie those obtained after RunDoNotUpdate has been called).
     */
    virtual double GetNextActiveTension()=0;
};


#endif /*ABSTRACTCONTRACTIONMODEL_HPP_*/
