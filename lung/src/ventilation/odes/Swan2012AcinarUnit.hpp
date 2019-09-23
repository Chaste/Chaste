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

#ifndef SWAN2012ACINARUNIT_HPP_
#define SWAN2012ACINARUNIT_HPP_

#include "AbstractAcinarUnit.hpp"

/**
 * Implementation of the acinar dynamics model presented in Swan et. al. 2012
 *
 * "A computational model of the topographic distribution of ventilation in
 *  healthy human lungs"
 *
 *  Note, that this model is implemented to match the equations underlying
 *  the model presented in Swan 2012, but does not seek to implement
 *  the numerical coupling method used. Instead, the method described in
 *  Coleman 1977 is applied.
 */
class Swan2012AcinarUnit : public AbstractAcinarUnit
{
friend class TestAcinarUnitModels;

public:
    /** Create a new acinar unit.
     */
    Swan2012AcinarUnit();

    /** Virtual destructor */
    virtual ~Swan2012AcinarUnit();

    /**
     * Set the timestep to use for simulating this acinus.
     *
     * @param dt  the timestep
     */
    void SetTimestep(double dt);

    /**
     * Simulate this acinar's behaviour between the time interval [tStart, tEnd],
     * with timestemp #mDt, updating the internal state variable values.
     *
     * @param tStart  beginning of the time interval to simulate
     * @param tEnd  end of the time interval to simulate
     */
    virtual void SolveAndUpdateState(double tStart, double tEnd);

    /**
     * Simulate this acinar's behaviour between the time interval [tStart, tEnd],
     * with timestep #mDt. The flow into the acinus, and so its volume, are kept fixed.
     *
     * @param tStart The starting time
     * @param tEnd The ending time
     */
    virtual void ComputeExceptFlow(double tStart, double tEnd);

    /**
     * Update the flow to the acinus across a time interval.
     * Note that this will update the acinus' volume and other
     * state variables.
     *
     * @param tStart The starting time
     * @param tEnd The ending time
     */
    virtual void UpdateFlow(double tStart, double tEnd);


    /** Set the air flow
     * @param flow  new value
     */
    void SetFlow(double flow);

    /**
     * @return the current value of the airflow, as given
     * in our state variable vector.
     */
    double GetFlow();

    /** Set the air pressure
     * @param pressure new value
     */
    void SetAirwayPressure(double pressure);

    /** Return the airway pressure
     * @return Airway pressure
     */
    double GetAirwayPressure();

    /** Set the pleural pressure
     * @param pressure new value
     */
    void SetPleuralPressure(double pressure);

    /** Set the resistance of the terminal bronchiole entering the acinar unit
     * @param raw new value
     */
    void SetTerminalBronchioleResistance(double raw);

    /**
     * @return the current stretch ratio of the acinus
     */
    double GetStretchRatio();

    /**
     * @param lambda The new stretch ratio
     */
    void SetStretchRatio(double lambda);

    /**
     * @return the current volume of the acinus
     */
    double GetVolume();

    /**
     * @param v0 The undeformed volume
     */
    void SetUndeformedVolume(double v0);

    /**
     * Calculates the static recoil pressure for the current stretch
     */
    double CalculateStaticRecoilPressure(double lambda);


private:
    /** The flow into the acinar unit */
    double mQ;

    /** The stretch ratio of the acinar unit */
    double mLambda;

    /** Timestep size */
    double mDt;

    /** The current air pressure in the acinar duct */
    double mPaw;

    /** The current pleural pressure (Pa) */
    double mPpl;

    /** The resistance of the terminal bronchiole entering the acinus */
    double mRaw;

    /** Fung material law constant */
    double mA;

    /** Fung material law constant */
    double mB;

    /** Fung material law constant (Pa) */
    double mXi;

    /** The initial volume of the acinus */
    double mV0;

    /**
     * @return dV / dlambda
     */
    double CalculateDerivativeVolumeByStrain();


    /**
     * @return d Pe /dlambda
     */
    double CalculateDerivativeStaticRecoilPressureByStrain();

    /**
     * @return gamma
     */
    double CalculateGamma(double lambda);
};


#endif /*SWAN2012ACINARUNIT*/
