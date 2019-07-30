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

#ifndef ABSTRACTACINARUNIT_HPP_
#define ABSTRACTACINARUNIT_HPP_


/**
 * This is the base class for ode-based acinar unit models
 */
class AbstractAcinarUnit
{

public:

    /**
     * Virtual destructor
     */
    virtual ~AbstractAcinarUnit() {};

    /**
     * Set the timestep to use for simulating this acinar unit's.
     *
     * @param dt  the timestep
     */
    virtual void SetTimestep(double dt) = 0;

    /**
     * Simulate this acinar unit's behaviour between the time interval [tStart, tEnd],
     * with timestep from SetTimestep, updating the internal state variable values.
     *
     * @param tStart  beginning of the time interval to simulate
     * @param tEnd  end of the time interval to simulate
     */
    virtual void SolveAndUpdateState(double tStart, double tEnd) = 0;

    /**
     * Simulate this acinar's behaviour between the time interval [tStart, tEnd],
     * with timestep #mDt. The flow into the acinus, and so its volume, are kept fixed.
     *
     * @param tStart The starting time
     * @param tEnd The ending time
     */
    virtual void ComputeExceptFlow(double tStart, double tEnd) = 0;

    /**
     * Update the flow to the acinus across a time interval.
     * Note that this will update the acinus' volume and other
     * state variables.
     *
     * @param tStart The starting time
     * @param tEnd The ending time
     */
    virtual void UpdateFlow(double tStart, double tEnd) = 0;

        /** Set the air flow
     * @param flow  new value
     */
    virtual void SetFlow(double flow) = 0;

    /**
     * @return the current value of the airflow.
     */
    virtual double GetFlow() = 0;

    /** Set the air pressure
     * @param pressure new value
     */
    virtual void SetAirwayPressure(double pressure) = 0;

    /** Return the airway pressure
     * @return Airway pressure
     */
    virtual double GetAirwayPressure() = 0;

    /** Set the pleural pressure
     * @param pressure new value
     */
    virtual void SetPleuralPressure(double pressure) = 0;

    /** Set the resistance of the terminal bronchiole entering the acinar unit
     * @param raw new value
     */
    virtual void SetTerminalBronchioleResistance(double raw) = 0;

    /**
     * @return the current stretch ratio of the acinus
     */
    virtual double GetStretchRatio() = 0;

    /**
     * @param lambda The new stretch ratio
     */
    virtual void SetStretchRatio(double lambda) = 0;

    /**
     * @return the current volume of the acinus
     */
    virtual double GetVolume() = 0;
};

#endif /*ABSTRACTACINARUNIT_HPP_*/
