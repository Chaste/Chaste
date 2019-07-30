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

#ifndef HIORNSAIRWAYWALL_HPP_
#define HIORNSAIRWAYWALL_HPP_

#include "AbstractAirwayWall.hpp"
#include "MathsCustomFunctions.hpp"

/**
 * Implements a nonlienar airway wall model based on the strain energy function in Hiorns et al. 2014 "Nonlinear Compliance
 * Modulates Dynamic Bronchoconstriction in a Multiscale Airway Model" -- see equation S18 in the paper
 */
class HiornsAirwayWall : public AbstractAirwayWall
{
    friend class TestAirwayWallModels;

public:

    /**
     * Constructor
     */
    HiornsAirwayWall();

    /**
     * Virtual destructor
     */
    virtual ~HiornsAirwayWall();

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
     * @param RIn The parameter RIn
     */
    void SetRIn(double);

     /**
     * @param ROut The parameter ROut
     */
    void SetROut(double);

     /**
     * @param mu The parameter mu
     */
    void Setmu(double);

     /**
     * @param phi1 The parameter phi1
     */
    void Setphi1(double);

     /**
     * @param phi2 The parameter phi2
     */
    void Setphi2(double);

     /**
     * @param C1 The parameter C1
     */
    void SetC1(double);

     /**
     * @param C2 The parameter C2
     */
    void SetC2(double);

     /**
     * @param A The parameter A
     */
    void SetA(double);

    /**
    * Works out the difference between the pressure and the pressure needed to produce a particular radius
    */
    double CalculatePressureRadiusResidual(double radius);

private:
    /** The airway radius at zero transpulmonary pressure */
    double mRi;

    /** The maximum airway radius */
    double mRiMax;

    /** The pressure inside the airway */
    double mAirwayPressure;

    /** The pleural pressure */
    double mPleuralPressure;

    /** The current airway radius */
    double mDeformedAirwayRadius;

    /** The outward pointing pressure applied to the airway wall*/
    double mTargetPressure;

    /** The airway inner radius pre deformation*/
    double mRIn;

    /** The airway outer radius pre deformation*/
    double mROut;

    /** The parameter mu from Hiorns et al. 2014 Biophys J*/
    double mmu;

    /** The parameter phi1 from Hiorns et al. 2014 Biophys J - Angle of collagen to hoop direction */
    double mphi1;

    /** The parameter phi2 from Hiorns et al. 2014 Biophys J - Angle of ASM to hoop direction */
    double mphi2;

    /** The parameter C1 from Hiorns et al. 2014 Biophys J*/
    double mC1;

    /** The parameter C2 from Hiorns et al. 2014 Biophys J*/
    double mC2;

    /** The parameter A from Hiorns et al. 2014 Biophys J*/
    double mA;
};

#endif /*HIORNSAIRWAYWALL_HPP_*/
