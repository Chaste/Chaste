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

#ifndef LAMBERTAIRWAYWALL_HPP_
#define LAMBERTAIRWAYWALL_HPP_

#include "AbstractAirwayWall.hpp"
#include "MathsCustomFunctions.hpp"

/**
 * Implements a simple dynamic airway wall model based on the equations in Lambert et al.
 * "A Computational Model of Expiratory Flow" Journal of Applied Physiology 1982
 */
class LambertAirwayWall : public AbstractAirwayWall
{
    friend class TestAirwayWallModels;

public:

    /**
     * Constructor
     */
    LambertAirwayWall();

    /**
     * Virtual destructor
     */
    virtual ~LambertAirwayWall();

    /**
     * Set the timestep to use for simulating the airway wall.
     *
     * @param dt  the timestep
     */
    virtual void SetTimestep(double dt);

    /**
     * Simulate this airway wall's behaviour between the time interval [tStart, tEnd],
     * with timestep from SetTimestep, updating the internal state variable values.
     *
     * @param tStart  beginning of the time interval to simulate
     * @param tEnd  end of the time interval to simulate
     */
    virtual void SolveAndUpdateState(double tStart, double tEnd);

    /**
     * @return the current value of the lumen radius.
     */
    virtual double GetLumenRadius();

    /** Set the air pressure
     * @param pressure new value
     */
    virtual void SetAirwayPressure(double pressure);

    /** Set the pleural pressure
     * @param pressure new value
     */
    virtual void SetPleuralPressure(double pressure);

    /**
     * @param n1 The parameter n1
     */
    void SetN1(double n1);

    /**
     * @param n2 The parameter n2
     */
    void SetN2(double n2);

    /**
     * @param p1 The parameter P1
     *
     * Parameter p2 will be automatically set based on p1.
     */
    void SetP1(double p1);

    /**
     * @param rmax The parameter alphai
     */
    void SetRiMax(double riMax);

    /**
     * @param R1 The parameter R1
     */
    void SetRi(double Ri);


    /** The airway radius at zero transpulmonary pressure */
    double mRi;

    /** The maximum airway radius */
    double mRiMax;

    /** The parameter P1 */
    double mP1;

    /** The parameter P2 */
    double mP2;

    /** The parameter n1 */
    double mN1;

    /** The parameter n2 */
    double mN2;


private:

    /** The pressure inside the airway */
    double mAirwayPressure;

    /** The pleural pressure */
    double mPleuralPressure;

    /** The current airway radius */
    double mDeformedAirwayRadius;
};


#endif /*LAMBERTAIRWAYWALL_HPP_*/
