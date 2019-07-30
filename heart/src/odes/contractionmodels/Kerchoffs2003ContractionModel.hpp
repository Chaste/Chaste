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

#ifndef KERCHOFFS2003CONTRACTIONMODEL_HPP_
#define KERCHOFFS2003CONTRACTIONMODEL_HPP_

#include "AbstractOdeBasedContractionModel.hpp"
#include "OdeSystemInformation.hpp"
#include <math.h>

/**
 *  Implementation of the ODE-based, stretch-dependent, stretch-rate-INdependent, contraction
 *  model detailed in the appendix of Kerchoffs 2003 "Intra- and interventricular asynchrony of
 *  electromechanics in the ventricularly paced heart".
 */
class Kerchoffs2003ContractionModel : public AbstractOdeBasedContractionModel
{
friend class TestContractionModels;

private:
    static const double a6;  /**< See reference. 2.0 um^{-1} */
    static const double a7;  /**< See reference. 1.5 um */
    static const double T0;  /**< See reference. 180 kPa */
    static const double Ea;  /**< See reference. 20 um^{-1} */
    static const double v0;  /**< See reference. 0.0075 um/ms */
    static const double ls0; /**< See reference. 1.9 um */
    //static const double tr;  /*< See reference. 75 ms.  Is now a parameter*/
    //static const double td;  /*< See reference. 75 ms.  Is now a parameter */
    //static const double b;   /*< See reference. 150 ms/um.  Is now a parameter */
    static const double ld;  /**< See reference. -0.4 um */

    /**Voltage threshold above which the cell is activated (mV) - note hysteresis*/
    static const double mActivationVoltage; //= 0;
    /**Voltage threshold below which the cell is deactivated (mV)*/
    static const double mDeactivationVoltage;// = -70;

    /** Length of the sarcomere in um. Variable "ls" in reference. Fibre-stretch is ls/ls0. */
    double mSarcomereLength;
    /** Time (ms) of electrical activation (= time the voltage at this cell reached mActivationVoltage) */
    double mActivationTime;
    /** Whether the cell is activated - whether the voltage has gone above mActivationTime
     *  and not gone below without going below mDeactivationVoltage, and
     *  that the cell has stopped producing force. */
    bool mIsActivated;
    /** Whether the cell is electrically unactivated yet (whether the voltage has gone
     *  below mDeactivationVoltage */
    bool mElectricallyUnactivated;

    /**
     *  @return the active tension as a function of length of contractile element. This is private. The public
     *  GetActiveTension() calls this using the value of lc in the state variable
     *  @param lengthOfContractileElement length of contractile element (the state variable in this model).
     */
    double GetActiveTension(double lengthOfContractileElement);


public:
    /** Constructor */
    Kerchoffs2003ContractionModel();


    /**
     *  The derivative function of the one state variable: "lc" in reference, the length of the contractile element
     *  @param time time
     *  @param rY 1D vector containing lc
     *  @param rDY 1D vector in which dlc/dt is set
     */
    void EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY);

    /**
     *  Set the input parameters. The calcium concentration is not used. If the voltage is such that activation
     *  has occured (gone above mActivationVoltage), the state is set to active and the time saved as the activation time.
     *  @param rInputParameters reference to the input parameters
     */
    void SetInputParameters(ContractionModelInputParameters& rInputParameters);

    /**
     *  Take the stretch and compute the sarcomere length (stretch rate is not used).
     *  @param stretch stretch
     *  @param stretchRate stretch rate
     */
    void SetStretchAndStretchRate(double stretch, double stretchRate);

    /**
     *  @return the active tension (note: actually a stress), ie kPa
     */
    double GetActiveTension();

    /**
     *  @return whether model is stretch-dependent
     */
    bool IsStretchDependent()
    {
        return true;
    }

    /**
     *  @return whether model is stretch-rate-independent
     */
    bool IsStretchRateDependent()
    {
        return false;
    }

    /**
     *  @return the active tension corresponding to the temporary stored state variables
     *  produced by calling RunDoNotUpdate (and before calling UpdateStateVariables())
     */
    double GetNextActiveTension();
};


#endif /*KERCHOFFS2003CONTRACTIONMODEL_HPP_*/
