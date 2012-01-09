/*

Copyright (C) University of Oxford, 2005-2012

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
     *  Get the current active tension
     */
    double GetActiveTension()
    {
        return rGetStateVariables()[0];
    }

    /**
     *  Get the active tension corresponding to the temporary stored state variables
     *  produced by calling RunDoNotUpdate (and before calling UpdateStateVariables())
     */
    double GetNextActiveTension()
    {
        return mTemporaryStateVariables[0];
    }


    /**
     *  This model is stretch-independent
     */
    bool IsStretchDependent()
    {
        return false;
    }

    /**
     *  This model is stretch-rate-independent
     */
    bool IsStretchRateDependent()
    {
        return false;
    }
};



#endif /*NASH2004CONTRACTIONMODEL_HPP_*/
