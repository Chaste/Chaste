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

#ifndef NASH2004CONTRACTIONMODEL_HPP_
#define NASH2004CONTRACTIONMODEL_HPP_

#include "AbstractOdeBasedContractionModel.hpp"
#include "OdeSystemInformation.hpp"
#include "Exception.hpp"

/**
 *  Nash2004 contraction model - Nash & Panfilov "Electromechanical mode of excitable tissue to
 *  study reentrant cardiac arrthrymias", Progress in Biophysics and Molecular Biology, 2004.
 *
 *  A simple, stretch- and stretch-rate-independent contraction model, just dependent on the voltage.
 *  If v is the non-dimensionalised voltage, Ta is given by
 *  dTa/dt = eps(v) (kTa*v - Ta);
 *  where eps(v) = e0 if v < 0.05, = 10*e0 v >= 0.05.
 *
 *  We use the non-dimensionalisation: V in [-85,40] ---> v in [0,1], ie v=(V+85)/125
 *
 *  Not sure what the appropriate value of e0 is: as the paper uses non-dimensionalised time
 *  the e0 above corresponds to "eps0/t0", where eps0 is value used in the paper (1.0) and
 *  t0 is the characteristic time. The paper suggests using t0=25.9 (?), which gives Ta growing too
 *  quickly at the beginning (as rapidly as the voltage). t0 = 100 coded at the moment..
 *
 */
class Nash2004ContractionModel : public AbstractOdeBasedContractionModel
{
    /** Stiffness parameter. See reference. kPa */
    static const double kTa;

    /** Other parameter. See above and reference. */
    static const double e0ByT0;

    /** Non-dimensionalised voltage. See above. */
    double mScaledVoltage;

public:
    /** Constructor */
    Nash2004ContractionModel() : AbstractOdeBasedContractionModel(1)
    {
        this->mpSystemInfo = OdeSystemInformation<Nash2004ContractionModel>::Instance();

        mScaledVoltage = 0.0;
        this->mStateVariables.push_back(0.0);
    }

    /**
     *  Calculate the derivative of the Ta (the state variable)
     *  @param time time
     *  @param rY 1D vector containing Ta
     *  @param rDY 1D vector in which dTa/dt is set
     */
    void EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY)
    {
        double epsilon = (mScaledVoltage < 0.05 ? e0ByT0 : 10*e0ByT0);
        rDY[0] = epsilon * (kTa*mScaledVoltage - rY[0]);
    }

    /**
     *  Set the input parameters. Only the voltage is used
     *  @param rInputParameters reference to the input parameters
     */
    void SetInputParameters(ContractionModelInputParameters& rInputParameters)
    {
        assert(rInputParameters.voltage != DOUBLE_UNSET);
        mScaledVoltage = (rInputParameters.voltage+85)/125;
    }

    /**
     *  Neither stretch nor stretch rate are used so this method does nothing
     *  @param stretch stretch
     *  @param stretchRate stretch rate
     */
    void SetStretchAndStretchRate(double stretch, double stretchRate)
    {
    }

    /**
     *  @return the current active tension
     */
    double GetActiveTension()
    {
        return rGetStateVariables()[0];
    }

    /**
     *  @return the active tension corresponding to the temporary stored state variables
     *  produced by calling RunDoNotUpdate (and before calling UpdateStateVariables())
     */
    double GetNextActiveTension()
    {
        return mTemporaryStateVariables[0];
    }


    /**
     *  @return whether model is stretch-independent
     */
    bool IsStretchDependent()
    {
        return false;
    }

    /**
     *  @return whether model is stretch-rate-independent
     */
    bool IsStretchRateDependent()
    {
        return false;
    }
};

#endif /*NASH2004CONTRACTIONMODEL_HPP_*/
