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

#ifndef _TESTACINARUNITMODELS_HPP_
#define _TESTACINARUNITMODELS_HPP_

#include <cxxtest/TestSuite.h>

#include "SimpleBalloonAcinarUnit.hpp"
#include "SimpleBalloonExplicitAcinarUnit.hpp"
#include "SigmoidalAcinarUnit.hpp"
#include "Swan2012AcinarUnit.hpp"
#include "TimeStepper.hpp"
#include "MathsCustomFunctions.hpp"
#include "OutputFileHandler.hpp"
#include "AbstractOdeSystem.hpp"
#include "OdeSystemInformation.hpp"
#include "BackwardEulerIvpOdeSolver.hpp"

#include <boost/math/tools/roots.hpp>
#include <boost/bind.hpp>
#include <iomanip>

//#include "PetscSetupAndFinalize.hpp"

//Class to benchmark sigmoidal solution against
class MySigmoidalOde : public AbstractOdeSystem
{
public:
    MySigmoidalOde(double raw, double a, double b, double c, double d) : AbstractOdeSystem(1)
    {
        mpSystemInfo = OdeSystemInformation<MySigmoidalOde>::Instance();
        mRaw = raw;
        mA = a;
        mB = b;
        mC = c;
        mD = d;
    }

    void EvaluateYDerivatives(double time, const std::vector<double>& rY,
                              std::vector<double>& rDY)
    {
        assert(rY[0] > mA);

        //double pleural_pressure = -2400*(sin((M_PI)*time));
        double pleural_pressure = -750 - 250*sin(2*M_PI*(time - 0.25));

        rDY[0] = 1/mRaw*(-(mD*log(mB/(rY[0] - mA) - 1)) + mC + pleural_pressure);
    }

    double mRaw, mA, mB, mC, mD;
};

template<>
void OdeSystemInformation<MySigmoidalOde>::Initialise()
{
    this->mVariableNames.push_back("V");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(2.3918/1e3);

    this->mInitialised = true;
}

//Class to benchmark swan solution against
class MySwanOde : public AbstractOdeSystem
{
public:
    MySwanOde(double raw, double a, double b, double xi) : AbstractOdeSystem(1)
    {
        mpSystemInfo = OdeSystemInformation<MySwanOde>::Instance();
        mRaw = raw;
        mA = a;
        mB = b;
        mXi = xi;
    }

    void EvaluateYDerivatives(double time, const std::vector<double>& rY,
                              std::vector<double>& rDY)
    {
        double pleural_pressure = -750 - 250*sin(2*M_PI*(time - 0.25));

        double V0 = 1/1e3;
        double lambda = std::pow(rY[0]/V0, 1.0/3.0);

        double gamma = (3.0/4.0)*(3*mA + mB)*(lambda*lambda - 1)*(lambda*lambda - 1);
        double Pe = mXi*std::exp(gamma)/(2.0*lambda)*(3*mA + mB)*(lambda*lambda -1);

        rDY[0] = 1/mRaw*(-Pe - pleural_pressure);
    }

    double mRaw, mA, mB, mXi;
};

template<>
void OdeSystemInformation<MySwanOde>::Initialise()
{
    this->mVariableNames.push_back("V");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(2.3918/1e3);

    this->mInitialised = true;
}


class TestAcinarUnitModels: public CxxTest::TestSuite
{
public:

    void TestSimpleBalloonAcinarUnitInspiration()
    {
        double viscosity = 1.92e-5;               //Pa s
        double terminal_airway_radius = 0.05;   //m
        double terminal_airway_length  = 0.02;   //m
        double terminal_airway_resistance = 8*viscosity*terminal_airway_length/(M_PI*SmallPow(terminal_airway_radius, 4));

        SimpleBalloonAcinarUnit acinus;

        // Coverage
        TS_ASSERT_DELTA(acinus.GetStretchRatio(), 0.0, 1e-6);
        acinus.SetStretchRatio(15.3);
        TS_ASSERT_DELTA(acinus.GetStretchRatio(), 0.0, 1e-6);

        acinus.SetAirwayPressure(0.0);
        acinus.SetPleuralPressure(0.0);
        acinus.SetFlow(0.0);
        acinus.SetUndeformedVolume(0.0);
        double compliance = 0.1/98.0665/1e3;  //in m^3 / pa. Converted from 0.1 L/cmH2O per lung.

        acinus.SetCompliance(compliance);
        acinus.SetTerminalBronchioleResistance(terminal_airway_resistance);

        acinus.SolveAndUpdateState(0.0, 0.2);
        TS_ASSERT_DELTA(acinus.GetFlow(), 0.0, 1e-6);   //With no pressure change we expect no flow
        TS_ASSERT_DELTA(acinus.GetVolume(), 0.0, 1e-2); //With no pressure change we expect no volume change

        TimeStepper time_stepper(0.0, 1.0, 0.0001);
        acinus.SetAirwayPressure(0.0);
        double pleural_pressure = 0.0;

        double ode_volume = 0.0;
        double flow_integral = 0.0;

        while (!time_stepper.IsTimeAtEnd())
        {
            pleural_pressure = - 2400*(sin((M_PI)*(time_stepper.GetNextTime())));

            //Solve the acinar problem coupled to a single bronchiole
            acinus.SetPleuralPressure(pleural_pressure);

            //acinus.SolveAndUpdateState(time_stepper.GetTime(), time_stepper.GetNextTime());
            acinus.ComputeExceptFlow(time_stepper.GetTime(), time_stepper.GetNextTime());

            double airway_pressure = acinus.GetAirwayPressure();
            double flow = -airway_pressure/terminal_airway_resistance;

            flow_integral += (time_stepper.GetNextTime() - time_stepper.GetTime())*flow;

            acinus.SetFlow(flow);
            acinus.UpdateFlow(time_stepper.GetTime(), time_stepper.GetNextTime());

            //Solve the corresponding ODE problem using backward Euler for testing
            // dv/dt = -1/R*(V/C - (Paw - Ppl))
            // Discretise using backward euler and rearrange to obtain the below
            double dt = time_stepper.GetNextTimeStep();
            ode_volume = (ode_volume - dt*pleural_pressure/terminal_airway_resistance)/(1 + dt/(terminal_airway_resistance*compliance));

            TS_ASSERT_DELTA(acinus.GetVolume(), ode_volume, 1e-8);

            time_stepper.AdvanceOneTimeStep();
        }

        TS_ASSERT_DELTA(ode_volume, -compliance*pleural_pressure, 1e-8);
        TS_ASSERT_DELTA(acinus.GetVolume(), -compliance*pleural_pressure, 1e-8);
        TS_ASSERT_DELTA(flow_integral, -compliance*pleural_pressure, 1e-8);

        // Coverage
        TS_ASSERT_THROWS_NOTHING(acinus.SetTimestep(0.01));
    }

    void TestSimpleBalloonExplicitAcinarUnitInspiration()
    {
       double viscosity = 1.92e-5;               //Pa s
       double terminal_airway_radius = 0.05;   //m
       double terminal_airway_length  = 0.02;   //m
       double terminal_airway_resistance = 8*viscosity*terminal_airway_length/(M_PI*SmallPow(terminal_airway_radius, 4));

       SimpleBalloonExplicitAcinarUnit acinus;

       // Coverage
       TS_ASSERT_DELTA(acinus.GetStretchRatio(), 0.0, 1e-6);
       acinus.SetStretchRatio(15.3);
       TS_ASSERT_DELTA(acinus.GetStretchRatio(), 0.0, 1e-6);

       acinus.SetAirwayPressure(0.0);
       acinus.SetPleuralPressure(0.0);
       acinus.SetFlow(0.0);
       acinus.SetUndeformedVolume(0.0);
       double compliance = 0.1/98.0665/1e3;  //in m^3 / pa. Converted from 0.1 L/cmH2O per lung.

       acinus.SetCompliance(compliance);
       acinus.SetTerminalBronchioleResistance(terminal_airway_resistance);

       acinus.SolveAndUpdateState(0.0, 0.2);
       TS_ASSERT_DELTA(acinus.GetFlow(), 0.0, 1e-6);   //With no pressure change we expect no flow
       TS_ASSERT_DELTA(acinus.GetVolume(), 0.0, 1e-2); //With no pressure change we expect no volume change

       //Uncomment below to find time step bound
       //std::cout << 2*terminal_airway_radius*compliance << std::endl; abort();

       TimeStepper time_stepper(0.0, 0.01, 0.0000001); //Only solve for a very short time due to dt restriction, this test is mostly for coverage
       acinus.SetAirwayPressure(0.0);
       double pleural_pressure = 0.0;

       double ode_volume = 0.0;
       double flow_integral = 0.0;

       while (!time_stepper.IsTimeAtEnd())
       {
           pleural_pressure = - 2400*(sin((M_PI)*(time_stepper.GetNextTime())));

           //Solve the acinar problem coupled to a single bronchiole
           acinus.SetPleuralPressure(pleural_pressure);
           acinus.ComputeExceptFlow(time_stepper.GetTime(), time_stepper.GetNextTime());

           double airway_pressure = acinus.GetAirwayPressure();
           double flow = -airway_pressure/terminal_airway_resistance;

           flow_integral += (time_stepper.GetNextTime() - time_stepper.GetTime())*flow;

           acinus.SetFlow(flow);
           acinus.UpdateFlow(time_stepper.GetTime(), time_stepper.GetNextTime());

           //Solve the corresponding ODE problem using backward Euler for testing
           // dv/dt = -1/R*(V/C - (Paw - Ppl))
           // Discretise using backward euler and rearrange to obtain the below
           double dt = time_stepper.GetNextTimeStep();
           ode_volume = (ode_volume - dt*pleural_pressure/terminal_airway_resistance)/(1 + dt/(terminal_airway_resistance*compliance));

           TS_ASSERT_DELTA(acinus.GetVolume(), ode_volume, 1e-8);

           time_stepper.AdvanceOneTimeStep();
       }

       // Coverage
       TS_ASSERT_THROWS_NOTHING(acinus.SetTimestep(0.01));
   }

    void TestSigmoidalAcinarUnitInspiration()
    {
        double viscosity = 1.92e-5;               //Pa s
        double terminal_airway_radius = 0.05;   //m
        double terminal_airway_length  = 0.02;   //m
        double terminal_airway_resistance = 8*viscosity*terminal_airway_length/(M_PI*SmallPow(terminal_airway_radius, 4));

        double a = 2/1e3; //m^3 (RV 2L)
        double b = (6 - 2)/1e3; //m^3 (TLC 6L, RV 2L)
        double c = 667; //Pa (6.8 cmH2O)
        double d = 300; //Pa (3.8 cmH2O)

        SigmoidalAcinarUnit acinus;

        // Coverage
        TS_ASSERT_DELTA(acinus.GetStretchRatio(), 0.0, 1e-6);
        acinus.SetStretchRatio(15.3);
        TS_ASSERT_DELTA(acinus.GetStretchRatio(), 0.0, 1e-6);

        acinus.SetA(a);
        acinus.SetB(b);
        acinus.SetC(c);
        acinus.SetD(d);

        acinus.SetAirwayPressure(0.0);
        acinus.SetPleuralPressure(-1.0);
        acinus.SetFlow(0.0);
        acinus.SetUndeformedVolume(3.4573/1e3); //Nearly completely deflated acinus. Model is invalid if it becomes fully deflated.
        acinus.SetTerminalBronchioleResistance(terminal_airway_resistance);

        acinus.SolveAndUpdateState(0.0, 0.2);
        TS_ASSERT_DELTA(acinus.GetFlow(), 0.0, 1e-6);   //With no pressure change we expect no flow
        TS_ASSERT_DELTA(acinus.GetVolume(), 0.0, 1e-2); //With no pressure change we expect no volume change

        // Setup corresponding ODE for testing
        MySigmoidalOde my_ode(terminal_airway_resistance,a,b,c,d);
        BackwardEulerIvpOdeSolver euler_solver(1);
        std::vector<double> initial_condition;
        initial_condition.push_back(3.4573/1e3);
        OdeSolution solutions = euler_solver.Solve(&my_ode, initial_condition, 0, 2, 0.001, 0.001);

        TimeStepper time_stepper(0.0, 2.0, 0.001);

        unsigned i = 0;
        while (!time_stepper.IsTimeAtEnd())
        {
            double pleural_pressure = -750 - 250*sin(2*M_PI*(time_stepper.GetNextTime() - 0.25));

            TS_ASSERT_DELTA(acinus.GetVolume(), solutions.rGetSolutions()[i][0], 1e-5);
            ++i;

            acinus.SetPleuralPressure(pleural_pressure);
            acinus.ComputeExceptFlow(time_stepper.GetTime(), time_stepper.GetNextTime());

            double airway_pressure = acinus.GetAirwayPressure();
            double flow = -airway_pressure/terminal_airway_resistance;

            acinus.SetFlow(flow);
            acinus.UpdateFlow(time_stepper.GetTime(), time_stepper.GetNextTime());

            time_stepper.AdvanceOneTimeStep();
        }

        // Coverage
        TS_ASSERT_THROWS_NOTHING(acinus.SetTimestep(0.01));
    }

    void TestSwan2012AcinarUnit()
    {
        Swan2012AcinarUnit acinus;

        // Coverage
        TS_ASSERT_DELTA(acinus.GetFlow(), 0.0, 1e-6);
        TS_ASSERT_DELTA(acinus.GetStretchRatio(), 1.26, 1e-6);

        // Test against corresponding ODE
        double viscosity = 1.92e-5;               //Pa s
        double terminal_airway_radius = 0.005;   //m
        double terminal_airway_length  = 0.02;   //m
        double terminal_airway_resistance = 8*viscosity*terminal_airway_length/(M_PI*SmallPow(terminal_airway_radius, 4));

        acinus.SetFlow(0.0);
        acinus.SetAirwayPressure(0.0);
        acinus.SetPleuralPressure(0.0);
        acinus.SetUndeformedVolume(1/1e3);
        acinus.SetStretchRatio(std::pow(2, 1.0/3.0));
        acinus.SetTerminalBronchioleResistance(terminal_airway_resistance);

        MySwanOde my_ode(terminal_airway_resistance, 0.433, -0.611, 2500);
        BackwardEulerIvpOdeSolver euler_solver(1);
        std::vector<double> initial_condition;
        initial_condition.push_back(2/1e3);
        OdeSolution solutions = euler_solver.Solve(&my_ode, initial_condition, 0, 2, 0.001, 0.001);

        TimeStepper time_stepper(0.0, 2.0, 0.001);

        unsigned i = 0;
        while (!time_stepper.IsTimeAtEnd())
        {
            double pleural_pressure = -750 - 250*sin(2*M_PI*(time_stepper.GetNextTime() - 0.25));

            TS_ASSERT_DELTA(acinus.GetVolume(), solutions.rGetSolutions()[i][0], 1e-5);
            ++i;

            acinus.SetPleuralPressure(pleural_pressure);
            acinus.ComputeExceptFlow(time_stepper.GetTime(), time_stepper.GetNextTime());

            double airway_pressure = acinus.GetAirwayPressure();
            double flow = -airway_pressure/terminal_airway_resistance;

            acinus.SetFlow(flow);
            acinus.UpdateFlow(time_stepper.GetTime(), time_stepper.GetNextTime());

            time_stepper.AdvanceOneTimeStep();
        }

        // Coverage
        TS_ASSERT_THROWS_NOTHING(acinus.SetTimestep(0.01));
        TS_ASSERT_THROWS_NOTHING(acinus.SolveAndUpdateState(0.0, 1.0));
    }
};
#endif /*_TESTACINARUNITMODELS_HPP_*/

