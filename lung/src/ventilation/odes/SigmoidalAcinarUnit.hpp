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

#ifndef SIGMOIDALACINARUNIT_HPP_
#define SIGMOIDALACINARUNIT_HPP_

#include "AbstractAcinarUnit.hpp"

/**
 * Implementation of an acinar balloon model with a sigmoidal compliance curves.
 *
 * Obtained by rearranging & differentiating equation 1 given in Venegas et al.
 * "A comprehensive equation for the pulmonary pressure-volume curve" JAP 1998
 */
class SigmoidalAcinarUnit : public AbstractAcinarUnit
{
friend class TestAcinarUnitModels;

public:
    /** Create a new acinar unit.
     */
    SigmoidalAcinarUnit();

    /** Virtual destructor */
    virtual ~SigmoidalAcinarUnit();

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
     */
    virtual void ComputeExceptFlow(double tStart, double tEnd);

    /**
     * Update the flow to the acinus across a time interval.
     * Note that this will update the acinus' volume and other
     * state variables.
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

    /** Get the pressure at the airway
     * @return The pressure at the airway
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
     * @param a The minimum volume of the acinus
     */
    void SetA(double a);

    /**
     * @param b The maximum volume of the acinus
     */
    void SetB(double b);

    /**
     * @param c The pressure for which the acinus has maximal compliance
     */
    void SetC(double c);

    /**
     * @param d The pressure range in which most of the volume change takes place
     */
    void SetD(double d);


private:
    /** The flow into the acinar unit */
    double mQ;

    /** Timestep size */
    double mDt;

    /** The current air pressure in the acinar duct */
    double mPaw;

    /** The current pleural pressure (Pa) */
    double mPpl;

    /** The resistance of the terminal bronchiole entering the acinus */
    double mRaw;

    /** The volume of the acinus */
    double mV;

    /** numerical parameter */
    double mA;

    /** numerical parameter */
    double mB;

    /** numerical parameter */
    double mC;

    /** numerical parameter */
    double mD;

    /** @return the compliance at the current volume. */
    double CalculateCurrentCompliance();

    /** @return the derivative of the volume/compliance wrt to volume at the current volume. */
    double CalculateCurrentDerivativeVolumeOverCompliance();
};


#endif /*SIGMOIDALACINARUNIT*/
