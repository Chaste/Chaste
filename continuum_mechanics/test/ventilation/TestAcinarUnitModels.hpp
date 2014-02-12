/*

Copyright (c) 2005-2014, University of Oxford.
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

#ifndef _TESTACINARUNITMODELS_HPP_
#define _TESTACINARUNITMODELS_HPP_

#include <cxxtest/TestSuite.h>

#include "Swan2012AcinarUnit.hpp"
#include "TimeStepper.hpp"

//#include "PetscSetupAndFinalize.hpp"

class TestAcinarUnitModels: public CxxTest::TestSuite
{
public:

    void TestSwan2012AcinarUnitCalculateMethods() throw(Exception)
    {
        Swan2012AcinarUnit acinus;

        acinus.SetAirwayPressure(0.0);
        acinus.SetAirwayPressure(0.0);
        acinus.SetPleuralPressure(-490);
        acinus.SetPleuralPressure(-490);
        acinus.SetFlow(0.0);
        acinus.SetTerminalBronchioleResistance(1.0); //Determine appropriate values

        double V0 = 4.0;

        acinus.SetStretchRatio(1.0);
        acinus.SetUndeformedVolume(V0);
        TS_ASSERT_DELTA(acinus.GetVolume(), V0, 1e-6);

        double lambda = 1.4;
        acinus.SetStretchRatio(lambda);
        TS_ASSERT_DELTA(acinus.GetVolume(), V0*lambda*lambda*lambda, 1e-6);
        TS_ASSERT_DELTA(acinus.CalculateDerivativeVolumeByStrain(), 3*V0*lambda*lambda, 1e-6);

        TS_ASSERT_DELTA(acinus.GetStretchRatio(), lambda, 1e-6);

        double gamma = (3.0/4.0)*(3*0.433 + 0.611)*(lambda*lambda - 1)*(lambda*lambda - 1);
        TS_ASSERT_DELTA(acinus.CalculateGamma(), gamma, 1e-6);

        double xi = 2500;
        double dPedLambda = (3.0*xi/2.0)*(3*0.433 + 0.611)*(3*0.433 + 0.611)*(lambda*lambda - 1)*(lambda*lambda - 1)*std::exp(gamma) +
                             (xi/2.0)*(3*0.433 + 0.611)*(lambda*lambda + 1)*std::exp(gamma)/(lambda*lambda);

        TS_ASSERT_DELTA(acinus.CalculateDerivativeStaticRecoilPressureByStrain(), dPedLambda, 1e-6);
        TS_ASSERT_DELTA(acinus.CalculateAcinarTissueCompliance(), 3*V0*lambda*lambda/dPedLambda, 1e-6);

    }

    void TestSwan2012AcinarUnitInspiration() throw(Exception)
    {
        Swan2012AcinarUnit acinus;

        acinus.SetAirwayPressure(0.0);
        acinus.SetAirwayPressure(0.0);
        acinus.SetPleuralPressure(-490);
        acinus.SetPleuralPressure(-490);
        acinus.SetFlow(0.0);
        acinus.SetStretchRatio(1.0);
        acinus.SetUndeformedVolume(1.0);
        acinus.SetTerminalBronchioleResistance(0.0005); //Determine appropriate values

        acinus.SolveAndUpdateState(0.0, 0.2);
        TS_ASSERT_DELTA(acinus.GetFlow(), 0.0, 1e-6); //With no pressure change we expect no flow
        TS_ASSERT_DELTA(acinus.GetVolume(), 1.0, 1e-6); //With no pressure change we expect no volume change

        //Sinussoidal inspiration followed by fixed pleural pressure
        TimeStepper time_stepper(0.0, 2.0, 0.01);
        double old_compliance = -1;
        double old_flow = DBL_MAX;
        while (!time_stepper.IsTimeAtEnd())
        {
            if(time_stepper.GetNextTime() <= 1.0) //breath in
            {
                double pleural_pressure = -490 - 300*(1 + sin((M_PI/2)*(time_stepper.GetNextTime() - 1)));
                acinus.SetPleuralPressure(pleural_pressure);
                acinus.SolveAndUpdateState(time_stepper.GetTime(), time_stepper.GetNextTime());

                TS_ASSERT_LESS_THAN(old_compliance, acinus.CalculateAcinarTissueCompliance()); //Check compliance monotonicity
                old_compliance = acinus.CalculateAcinarTissueCompliance();
            }
            else //constant pleural pressure
            {
                acinus.SolveAndUpdateState(time_stepper.GetTime(), time_stepper.GetNextTime());

                if(time_stepper.GetNextTime() > 1.1) //Check flow decreases monotonically after the initial inspiration
                {
                    TS_ASSERT_LESS_THAN(acinus.GetFlow(), old_flow);
                    old_flow = acinus.GetFlow();
                }
            }

            time_stepper.AdvanceOneTimeStep();
        }
    }
};
#endif /*_TESTACINARUNITMODELS_HPP_*/

