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


#ifndef _TESTIONICMODELS_HPP_
#define _TESTIONICMODELS_HPP_

#include <cxxtest/TestSuite.h>
#include <cmath>
#include <iostream>
#include <vector>
#include <string>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "RunAndCheckIonicModels.hpp"
#include "Exception.hpp"
#include "HeartConfig.hpp"

#include "SimpleStimulus.hpp"
#include "RegularStimulus.hpp"
#include "MultiStimulus.hpp"
#include "ZeroStimulus.hpp"

#include "EulerIvpOdeSolver.hpp"
#include "BackwardEulerIvpOdeSolver.hpp"
#include "RungeKutta2IvpOdeSolver.hpp"
#include "RungeKutta4IvpOdeSolver.hpp"
#include "CellProperties.hpp"

#include "HodgkinHuxley1952.hpp"
#include "HodgkinHuxley1952BackwardEuler.hpp"
#include "FitzHughNagumo1961OdeSystem.hpp"
#include "LuoRudy1991.hpp"
#include "LuoRudy1991BackwardEuler.hpp"

#include "FoxModel2002.hpp"
#include "FoxModel2002BackwardEuler.hpp"

#include "FaberRudy2000.hpp"
#include "FaberRudy2000Opt.hpp"

#include "NobleVargheseKohlNoble1998a.hpp"
#include "NobleVargheseKohlNoble1998WithSac.hpp"
#include "NobleVargheseKohlNoble1998aOpt.hpp"
#include "NobleVargheseKohlNoble1998aBackwardEuler.hpp"
#include "Mahajan2008.hpp"
#include "Mahajan2008BackwardEuler.hpp"
#include "TenTusscher2006Epi.hpp"
#include "TenTusscher2006EpiBackwardEuler.hpp"
#include "DiFrancescoNoble1985.hpp"
#include "Maleckar2008.hpp"

#include "ArchiveLocationInfo.hpp"

//This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

// Note: RunOdeSolverWithIonicModel(), CheckCellModelResults(), CompareCellModelResults()
// are defined in RunAndCheckIonicModels.hpp

class TestIonicModels : public CxxTest::TestSuite
{
public:
    void TestSolveForNoble98WithSimpleStimulus(void)
    {
        clock_t ck_start, ck_end;

        // Set stimulus
        double magnitude_stimulus = -3/0.095;  // uA/cm2
        double duration_stimulus = 3;  // ms
        double start_stimulus = 10.0;   // ms
        boost::shared_ptr<SimpleStimulus> p_stimulus(new SimpleStimulus(
                magnitude_stimulus,
                duration_stimulus,
                start_stimulus));
        boost::shared_ptr<EulerIvpOdeSolver> p_solver; // No solver set yet
        double time_step = 0.01;

        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(time_step, time_step, time_step);

        CellNobleVargheseKohlNoble1998aFromCellML n98_ode_system(p_solver, p_stimulus);

        // Set solver after model creation, for coverage
        TS_ASSERT(!p_solver);
        TS_ASSERT(!n98_ode_system.GetSolver());
        TS_ASSERT(n98_ode_system.GetSolver() == p_solver);
        p_solver.reset(new EulerIvpOdeSolver);
        n98_ode_system.SetSolver(p_solver);
        TS_ASSERT(n98_ode_system.GetSolver());
        TS_ASSERT(n98_ode_system.GetSolver() == p_solver);

        // Solve and write to file
        ck_start = clock();
        RunOdeSolverWithIonicModel(&n98_ode_system,
                                   150.0,
                                   "N98RegResult",
                                   100,
                                   true,
                                   true);
        ck_end = clock();
        double forward = (double)(ck_end - ck_start)/CLOCKS_PER_SEC;
        std::cout << "\n\tForward: " << forward << std::endl;

        CheckCellModelResults("N98RegResult");
        TS_ASSERT_DELTA( n98_ode_system.GetIIonic(), 0.2462, 1e-3);

        std::string error_should_be = "Non fast-slow cell model being used in a fast-slow problem.";
        TS_ASSERT_THROWS_THIS(n98_ode_system.SetState(ALL_VARS), error_should_be);
        TS_ASSERT_THROWS_THIS(n98_ode_system.GetNumSlowValues(), error_should_be);
        TS_ASSERT_THROWS_THIS(n98_ode_system.IsFastOnly(), error_should_be);
        std::vector<double> slows;
        slows.push_back(100.0);
        TS_ASSERT_THROWS_THIS(n98_ode_system.AdjustOutOfRangeSlowValues(slows), error_should_be);
        TS_ASSERT_THROWS_THIS(n98_ode_system.GetSlowValues(slows), error_should_be);
        TS_ASSERT_THROWS_THIS(n98_ode_system.SetSlowValues(slows), error_should_be);

        TS_ASSERT(n98_ode_system.HasCellMLDefaultStimulus());
        boost::shared_ptr<RegularStimulus> p_stim = n98_ode_system.UseCellMLDefaultStimulus();
        TS_ASSERT_DELTA(p_stim->GetMagnitude(), -31.5789, 1e-4);
        TS_ASSERT_DELTA(p_stim->GetPeriod(), 1000, 1e-7);
        TS_ASSERT_DELTA(p_stim->GetStartTime(), 100, 1e-7);
        TS_ASSERT_DELTA(p_stim->GetDuration(), 3, 1e-7);

        // Coverage of some helper methods
        double new_value = -128.73;
        n98_ode_system.SetStateVariable("membrane_voltage", new_value);
        TS_ASSERT_DELTA(n98_ode_system.GetStateVariable("membrane_voltage"),new_value,1e-9);
    }

    void TestSolveForNoble98WithSacWithSimpleStimulus(void)
    {
        clock_t ck_start, ck_end;

        // Set stimulus
        double magnitude_stimulus = -3;  // uA/cm2
        double duration_stimulus = 3;  // ms
        double start_stimulus = 10.0;   // ms
        boost::shared_ptr<SimpleStimulus> p_stimulus(new SimpleStimulus(
                magnitude_stimulus,
                duration_stimulus,
                start_stimulus));
        boost::shared_ptr<EulerIvpOdeSolver> p_solver(new EulerIvpOdeSolver);
        double time_step = 0.01;

        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(time_step, time_step, time_step);

        CML_noble_varghese_kohl_noble_1998_basic_with_sac   n98_with_sac(p_solver, p_stimulus);

        // some models have this implemented so they can be used in mechanics simulations
        TS_ASSERT_DELTA(n98_with_sac.GetIntracellularCalciumConcentration(), 1.4e-5, 2e-6);

        // Solve and write to file
        ck_start = clock();
        RunOdeSolverWithIonicModel(&n98_with_sac,
                                   150.0,
                                   "N98SacResult");
        ck_end = clock();
        double forward = (double)(ck_end - ck_start)/CLOCKS_PER_SEC;
        std::cout << "\n\tForward: " << forward << std::endl;

        // the 'good' data the result is compared against here was just copied from the standard N98 run good data, as they should be identical
        CheckCellModelResults("N98SacResult");


        // get a new ODE system (which still has its state variables set to the initial conditions),
        // and check GetIonic agrees with standard noble98
        CML_noble_varghese_kohl_noble_1998_basic_with_sac   another_n98_with_sac(p_solver, p_stimulus);
        CellNobleVargheseKohlNoble1998aFromCellML   n98_ode_system(p_solver, p_stimulus);
        TS_ASSERT_DELTA( another_n98_with_sac.GetIIonic(), n98_ode_system.GetIIonic(), 1e-3);

        another_n98_with_sac.SetStretch(0.9);
        TS_ASSERT_DELTA( another_n98_with_sac.GetIIonic(), n98_ode_system.GetIIonic(), 1e-3);

        // add stretch, and now should be different
        another_n98_with_sac.SetStretch(1.1);
        TS_ASSERT_DELTA( another_n98_with_sac.GetIIonic(), -20.3267, 1e-3);

        // coverage
        TS_ASSERT_DELTA(another_n98_with_sac.GetStretch(),1.1,1e-9);
    }


    void TestSolveForNoble98WithSacStretchActivated(void)
    {
        clock_t ck_start, ck_end;

        // Set stimulus
        boost::shared_ptr<ZeroStimulus> p_stimulus(new ZeroStimulus);
        boost::shared_ptr<EulerIvpOdeSolver> p_solver(new EulerIvpOdeSolver);
        double time_step = 0.01;

        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(time_step, time_step, time_step);

        CML_noble_varghese_kohl_noble_1998_basic_with_sac   n98_with_sac(p_solver, p_stimulus);

        n98_with_sac.SetStretch(1.1);

        // Solve and write to file
        ck_start = clock();
        RunOdeSolverWithIonicModel(&n98_with_sac,
                                   150.0,
                                   "N98Sac_StretchActivatedResult");
        ck_end = clock();
        double forward = (double)(ck_end - ck_start)/CLOCKS_PER_SEC;
        std::cout << "\n\tForward: " << forward << std::endl;

        CheckCellModelResults("N98Sac_StretchActivatedResult");
    }

    void TestSolveForOptimisedNoble98WithSimpleStimulus(void)
    {
        clock_t ck_start, ck_end;

        // Set stimulus
        double magnitude_stimulus = -3/0.095;  // uA/cm2
        double duration_stimulus = 3;  // ms
        double start_stimulus = 10.0;   // ms
        boost::shared_ptr<SimpleStimulus> p_stimulus(new SimpleStimulus(magnitude_stimulus,
                                                                        duration_stimulus,
                                                                        start_stimulus));
        boost::shared_ptr<EulerIvpOdeSolver> p_solver(new EulerIvpOdeSolver);
        double time_step = 0.01;

        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(time_step, time_step, time_step);

        //Check Optimised
        CellNobleVargheseKohlNoble1998aFromCellMLOpt n98_ode_system(p_solver, p_stimulus);

        // Solve and write to file
        ck_start = clock();
        RunOdeSolverWithIonicModel(&n98_ode_system,
                                   150.0,
                                   "N98RegResultOpt");
        ck_end = clock();
        double forward = (double)(ck_end - ck_start)/CLOCKS_PER_SEC;
        std::cout << "\n\tForward: " << forward << std::endl;

        CheckCellModelResults("N98RegResultOpt", "N98RegResult");
        TS_ASSERT_DELTA( n98_ode_system.GetIIonic(), 0.2462, 1e-3);

        //Stress the lookup table with a silly voltage
        n98_ode_system.rGetStateVariables()[0] = 100.0;
        TS_ASSERT_EQUALS(n98_ode_system.GetVoltage(), 100.0);
        TS_ASSERT_THROWS_EQUALS( n98_ode_system.GetIIonic(), const Exception &err,
                err.GetShortMessage().find("membrane_voltage outside lookup table range",0), 0u);
        n98_ode_system.rGetStateVariables()[0] = 101.0;
        TS_ASSERT_THROWS_EQUALS( n98_ode_system.GetIIonic(), const Exception &err,
                err.GetShortMessage().find("membrane_voltage outside lookup table range",0), 0u);
        n98_ode_system.rGetStateVariables()[0] = 99.0;
        TS_ASSERT_THROWS_NOTHING( n98_ode_system.GetIIonic());
        n98_ode_system.rGetStateVariables()[0] = -100.1;
        TS_ASSERT_THROWS_EQUALS( n98_ode_system.GetIIonic(), const Exception &err,
                err.GetShortMessage().find("membrane_voltage outside lookup table range",0), 0u);
        n98_ode_system.rGetStateVariables()[0] = -100.0;
        TS_ASSERT_THROWS_NOTHING( n98_ode_system.GetIIonic());

        n98_ode_system.rGetStateVariables()[0] = -100.1;

        TS_ASSERT_THROWS_EQUALS( RunOdeSolverWithIonicModel(&n98_ode_system, 150.0, "DoNotRun"),
                const Exception &err, err.GetShortMessage().find("membrane_voltage outside lookup table range",0), 0u);
    }


    void TestSolverForHH52WithSimpleStimulus(void)
    {
        clock_t ck_start, ck_end;

        // Set stimulus
        double magnitude_stimulus = -20.0;  // uA/cm2
        double duration_stimulus = 0.5;  // ms
        double start_stimulus = 10.0;   // ms
        boost::shared_ptr<SimpleStimulus> p_stimulus(new SimpleStimulus(magnitude_stimulus,
                                                                        duration_stimulus,
                                                                        start_stimulus));
        boost::shared_ptr<EulerIvpOdeSolver> p_solver(new EulerIvpOdeSolver);
        CellHodgkinHuxley1952FromCellML hh52_ode_system(p_solver, p_stimulus);

        // Solve and write to file
        const double end_time = 150.0;
        ck_start = clock();
        RunOdeSolverWithIonicModel(&hh52_ode_system, end_time, "HH52RegResult");
        ck_end = clock();
        double forward = (double)(ck_end - ck_start)/CLOCKS_PER_SEC;
        std::cout << "\n\tForward: " << forward << std::endl;

        CheckCellModelResults("HH52RegResult");

        // test GetIionic: (the GetIionic method was first manually tested
        // by changing the EvaluateYDerivatives() code to call it, this verified
        // that GetIionic has no errors, therefore we can test here against
        // a hardcoded result
        RunOdeSolverWithIonicModel(&hh52_ode_system, 15.0, "HhGetIIonic");
        TS_ASSERT_DELTA(hh52_ode_system.GetIIonic(), 40.6341, 1e-3);

        // For coverage of the case where a cell model has no non-linear ODEs, we also test the backward Euler
        // version of this model.
        CellHodgkinHuxley1952FromCellMLBackwardEuler hh52_be(p_solver, p_stimulus);
        RunOdeSolverWithIonicModel(&hh52_be, end_time, "HH52BackwardEuler", 1, true /* check ComputeExceptVoltage too */);
        // Compare end result against using SolveAndUpdateState
        std::vector<double> state_variables_copy = hh52_be.GetStdVecStateVariables();
        hh52_be.ResetToInitialConditions();
        hh52_be.SolveAndUpdateState(0.0, end_time);
        for (unsigned i=0; i<hh52_be.GetNumberOfStateVariables(); i++)
        {
            TS_ASSERT_DELTA(hh52_be.GetStateVariable(i), state_variables_copy[i], 1e-6);
        }
    }


    void TestSolverForFHN61WithSimpleStimulus(void)
    {
        clock_t ck_start, ck_end;

        // Set stimulus
        double magnitude_stimulus = -80.0;   // dimensionless
        double duration_stimulus = 0.5;  // ms
        double start_stimulus = 0.0;   // ms
        boost::shared_ptr<SimpleStimulus> p_stimulus(new SimpleStimulus(magnitude_stimulus,
                                                                        duration_stimulus,
                                                                        start_stimulus));
        boost::shared_ptr<EulerIvpOdeSolver> p_solver(new EulerIvpOdeSolver);
        FitzHughNagumo1961OdeSystem fhn61_ode_system(p_solver, p_stimulus);

        // fhn has no [Ca_i]
        TS_ASSERT_THROWS_THIS(fhn61_ode_system.GetIntracellularCalciumConcentration(),
                "AbstractCardiacCellInterface::GetIntracellularCalciumConcentration() called. "
                "Either model has no [Ca_i] or method has not been implemented yet");

        // Solve and write to file
        ck_start = clock();
        RunOdeSolverWithIonicModel(&fhn61_ode_system,
                                   500.0,
                                   "FHN61RegResult");
        ck_end = clock();
        double forward = (double)(ck_end - ck_start)/CLOCKS_PER_SEC;
        std::cout << "\n\tForward: " << forward << std::endl;

        CheckCellModelResults("FHN61RegResult");

        // test GetIionic ('fake' ionic current) (the GetIionic method was first
        // manually tested by changing the EvaluateYDerivatives() code to call it,
        // this verified that GetIionic has no errors, therefore we can test here
        // against a hardcoded result
        TS_ASSERT_DELTA( fhn61_ode_system.GetIIonic(), -0.0058, 1e-3);

        // some coverage
        boost::shared_ptr<SimpleStimulus> p_another_stimulus(new SimpleStimulus(-200, 1.0, 0.0));
        boost::shared_ptr<SimpleStimulus> p_intra_stimulus(new SimpleStimulus(-100, 1.0, 0.0));
        FitzHughNagumo1961OdeSystem another_fhn61_ode_system(p_solver, p_stimulus);

        another_fhn61_ode_system.SetStimulusFunction(p_another_stimulus);
        TS_ASSERT_DELTA(another_fhn61_ode_system.GetStimulus(0.5), -200, 1e-12);
        TS_ASSERT_DELTA(another_fhn61_ode_system.GetIntracellularStimulus(0.5), -200, 1e-12);

        another_fhn61_ode_system.SetIntracellularStimulusFunction(p_intra_stimulus);
        TS_ASSERT_DELTA(another_fhn61_ode_system.GetStimulus(0.5), -100, 1e-12);
        TS_ASSERT_DELTA(another_fhn61_ode_system.GetIntracellularStimulus(0.5), -100, 1e-12);

        TS_ASSERT_EQUALS(another_fhn61_ode_system.HasCellMLDefaultStimulus(),false);
        TS_ASSERT_THROWS_THIS(another_fhn61_ode_system.UseCellMLDefaultStimulus(),"This class has no default stimulus from CellML metadata.");

    }


    void TestSolverForLR91WithDelayedSimpleStimulus(void)
    {
        clock_t ck_start, ck_end;

        // Set stimulus
        double magnitude = -25.5;
        double duration  = 2.0 ;  // ms
        double when = 50.0; // ms

        boost::shared_ptr<SimpleStimulus> p_stimulus(new SimpleStimulus(magnitude, duration, when));
        boost::shared_ptr<EulerIvpOdeSolver> p_solver(new EulerIvpOdeSolver);

        std::cout << "\n";
        std::cout << "p_solver.use_count() = " << p_solver.use_count() << std::endl;
        std::cout << "p_stimulus.use_count() = " << p_stimulus.use_count() << std::endl;
        std::cout << "\n";

        double end_time = 1000.0; //One second in milliseconds

        CellLuoRudy1991FromCellML lr91_ode_system(p_solver, p_stimulus);
        TS_ASSERT_EQUALS(lr91_ode_system.GetVoltageIndex(), 0u); // For coverage

        std::cout << "\n";
        std::cout << "p_solver.use_count() = " << p_solver.use_count() << std::endl;
        std::cout << "p_stimulus.use_count() = " << p_stimulus.use_count() << std::endl;
        std::cout << "\n";

        // Solve and write to file
        ck_start = clock();
        RunOdeSolverWithIonicModel(&lr91_ode_system,
                                   end_time,
                                   "Lr91DelayedStim");
        ck_end = clock();
        double forward = (double)(ck_end - ck_start)/CLOCKS_PER_SEC;
        std::cout << "\n\tForward: " << forward << std::endl;

        std::cout << "\n";
        std::cout << "p_solver.use_count() = " << p_solver.use_count() << std::endl;
        std::cout << "p_stimulus.use_count() = " << p_stimulus.use_count() << std::endl;
        std::cout << "\n";

        CheckCellModelResults("Lr91DelayedStim");

        // test GetIionic: (the GetIionic method was first manually tested
        // by changing the EvaluateYDerivatives() code to call it, this verified
        // that GetIionic has no errors, therefore we can test here against
        // a hardcoded result
        RunOdeSolverWithIonicModel(&lr91_ode_system,
                                   60.0,
                                   "Lr91GetIIonic");
        TS_ASSERT_DELTA( lr91_ode_system.GetIIonic(), 1.9411, 1e-3);

        // For coverage
        lr91_ode_system.ResetToInitialConditions();
        std::vector<double> inits = lr91_ode_system.GetInitialConditions();
        for (unsigned i=0; i<inits.size(); i++)
        {
            TS_ASSERT_EQUALS(lr91_ode_system.rGetStateVariables()[i], inits[i]);
        }

        TS_ASSERT_DELTA(lr91_ode_system.GetParameter(1u),23,1e-12);
    }

    void TestSolverForLR91WithRegularStimulus(void)
    {
        // Set stimulus
        double magnitude = -25.5;
        double duration  = 2.0 ;  // ms
        double start = 50.0; // ms
        double period = 500; // ms
        boost::shared_ptr<RegularStimulus> p_stimulus(new RegularStimulus(magnitude,
                                                                          duration,
                                                                          period,
                                                                          start));

        double end_time = 1000.0; //One second in milliseconds

        boost::shared_ptr<EulerIvpOdeSolver> p_solver(new EulerIvpOdeSolver);
        CellLuoRudy1991FromCellML lr91_ode_system(p_solver, p_stimulus);

        // some models have this implemented so they can be used in mechanics simulations
        TS_ASSERT_DELTA(lr91_ode_system.GetIntracellularCalciumConcentration(), 0.0002, 1e-5)

        // Solve and write to file
        RunOdeSolverWithIonicModel(&lr91_ode_system,
                                   end_time,
                                   "Lr91RegularStim");

        CheckCellModelResults("Lr91RegularStim");

        // Test SolveAndUpdateState
        double v = lr91_ode_system.GetVoltage();
        lr91_ode_system.ResetToInitialConditions();
        lr91_ode_system.SolveAndUpdateState(0.0, end_time);
        TS_ASSERT_DELTA(lr91_ode_system.GetVoltage(), v, 1e-10);
    }

    void TestBackwardEulerLr91WithDelayedSimpleStimulus(void)
    {
        clock_t ck_start, ck_end;
        // Set stimulus
        double magnitude = -25.5;
        double duration  = 2.0 ;  // ms
        double when = 50.0; // ms

        boost::shared_ptr<SimpleStimulus> p_stimulus(new SimpleStimulus(magnitude, duration, when));
        boost::shared_ptr<EulerIvpOdeSolver> p_solver(new EulerIvpOdeSolver);
        double end_time = 1000.0; //One second in milliseconds

        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.01, 0.01, 0.01);

        // Solve using backward euler (and lookup tables)
        CellLuoRudy1991FromCellMLBackwardEuler lr91_backward_euler(p_solver, p_stimulus);

        // Some models have this implemented so they can be used in mechanics simulations
        TS_ASSERT_DELTA(lr91_backward_euler.GetIntracellularCalciumConcentration(), 0.0002, 1e-5);

        ck_start = clock();
        RunOdeSolverWithIonicModel(&lr91_backward_euler,
                                   end_time,
                                   "Lr91BackwardEuler",
                                   100,
                                   true,
                                   true);
        ck_end = clock();
        double backward1 = (double)(ck_end - ck_start)/CLOCKS_PER_SEC;

        // Solve using forward Euler
        CellLuoRudy1991FromCellML lr91_ode_system(p_solver, p_stimulus);
        ck_start = clock();
        RunOdeSolverWithIonicModel(&lr91_ode_system,
                                   end_time,
                                   "Lr91DelayedStim");
        ck_end = clock();
        double forward = (double)(ck_end - ck_start)/CLOCKS_PER_SEC;

        // Compare results
        CompareCellModelResults("Lr91DelayedStim", "Lr91BackwardEuler", 0.01);

        // Test SolveAndUpdateState
        double v = lr91_backward_euler.GetVoltage();
        lr91_backward_euler.ResetToInitialConditions();
        lr91_backward_euler.SolveAndUpdateState(0.0, end_time);
        TS_ASSERT_DELTA(lr91_backward_euler.GetVoltage(), v, 1e-10);

        // Try with larger timestep and coarser tolerance.
        // We can't use a larger time step than 0.01 for forward Euler - the gating
        // variables go out of range.
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.1, 0.1, 0.1);
        CellLuoRudy1991FromCellMLBackwardEuler lr91_backward_euler2(p_solver, p_stimulus);
        ck_start = clock();
        RunOdeSolverWithIonicModel(&lr91_backward_euler2,
                                   end_time,
                                   "Lr91BackwardEuler2");

        ck_end = clock();
        double backward2 = (double)(ck_end - ck_start)/CLOCKS_PER_SEC;
        CompareCellModelResults("Lr91DelayedStim", "Lr91BackwardEuler2", 0.25);

        std::cout << "Run times:\n\tForward: " << forward << "\n\tBackward: "
                  << backward1 << "\n\tBackward (long dt): " << backward2 << std::endl;


        // cover and check GetIIonic() match for normal and backward euler lr91
        // calc IIonic using initial conditions
        CellLuoRudy1991FromCellML lr91(p_solver, p_stimulus);
        CellLuoRudy1991FromCellMLBackwardEuler backward_lr91(p_solver, p_stimulus);
        TS_ASSERT_DELTA(lr91.GetIIonic(), backward_lr91.GetIIonic(), 1e-3);

        // Reset for next test
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.01, 0.01, 0.01);
    }

    void TestSolverForFR2000WithDelayedSimpleStimulus(void)
    {
        clock_t ck_start, ck_end;

        // Set stimulus
        double magnitude = -25.5;
        double duration  = 2.0;  // ms
        double when = 10.0; // ms
        boost::shared_ptr<SimpleStimulus> p_stimulus(new SimpleStimulus(magnitude, duration, when));

        double end_time = 1000.0; //ms
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.007, 0.007, 0.007);

        boost::shared_ptr<EulerIvpOdeSolver> p_solver(new EulerIvpOdeSolver);

        CellFaberRudy2000FromCellMLOpt fr2000_ode_system_opt(p_solver, p_stimulus);
        CellFaberRudy2000FromCellML fr2000_ode_system(p_solver, p_stimulus);

        // Solve and write to file
        ck_start = clock();
        RunOdeSolverWithIonicModel(&fr2000_ode_system,
                                   end_time,
                                   "FR2000DelayedStim",
                                   500, false);
        ck_end = clock();
        double forward = (double)(ck_end - ck_start)/CLOCKS_PER_SEC;

        ck_start = clock();
        RunOdeSolverWithIonicModel(&fr2000_ode_system_opt,
                                   end_time,
                                   "FR2000DelayedStimOpt",
                                   500, false);
        ck_end = clock();
        double opt = (double)(ck_end - ck_start)/CLOCKS_PER_SEC;

        std::cout << "\n\tForward: " << forward
                  << "\n\tOptimised: " << opt << std::endl;

        CheckCellModelResults("FR2000DelayedStim");
        CompareCellModelResults("FR2000DelayedStim", "FR2000DelayedStimOpt", 2e-3);

        // test GetIionic: the GetIionic method was first manually tested
        // by changing the EvaluateYDerivatives() code to call it, this verified
        // that GetIionic has no errors, therefore we can test here against
        // a hardcoded result
        TS_ASSERT_DELTA(fr2000_ode_system.GetIIonic(), 0.0002, 1e-4);
        TS_ASSERT_DELTA(fr2000_ode_system_opt.GetIIonic(), 0.0002, 1e-4);

        // Check that ComputeExceptVoltage does the correct thing (doesn't change the voltage)
        double voltage = fr2000_ode_system.GetVoltage();
        fr2000_ode_system.ComputeExceptVoltage(end_time, end_time+0.001);
        TS_ASSERT_DELTA(fr2000_ode_system.GetVoltage(), voltage, 1e-10);
        voltage = fr2000_ode_system_opt.GetVoltage();
        fr2000_ode_system_opt.ComputeExceptVoltage(end_time, end_time+0.001);
        TS_ASSERT_DELTA(fr2000_ode_system_opt.GetVoltage(), voltage, 1e-10);

    }


    void TestSolverForFR2000WithVariablePotassiumCurrents(void)
    {
        // Set stimulus
        double magnitude = -25.5;
        double duration  = 2.0;  // ms
        double when = 0.0; // ms
        boost::shared_ptr<SimpleStimulus> p_stimulus(new SimpleStimulus(magnitude, duration, when));

        double end_time = 1000.0; //ms
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.007, 0.007, 0.007);

        boost::shared_ptr<EulerIvpOdeSolver> p_solver(new EulerIvpOdeSolver);
        CellFaberRudy2000FromCellML fr2000_ode_system_endo(p_solver, p_stimulus);
        fr2000_ode_system_endo.SetParameter("ScaleFactorGks",0.462);
        fr2000_ode_system_endo.SetParameter("ScaleFactorIto",0.0);
        fr2000_ode_system_endo.SetParameter("ScaleFactorGkr",1.0);
        // Solve and write to file
        RunOdeSolverWithIonicModel(&fr2000_ode_system_endo,
                                   end_time,
                                   "FR2000Endo",
                                   500, false);

        CheckCellModelResults("FR2000Endo");

        CellFaberRudy2000FromCellML fr2000_ode_system_mid(p_solver, p_stimulus);
        fr2000_ode_system_mid.SetParameter("ScaleFactorGks",1.154);
        fr2000_ode_system_mid.SetParameter("ScaleFactorIto",0.85);
        fr2000_ode_system_mid.SetParameter("ScaleFactorGkr",1.0);

        // Solve and write to file
        RunOdeSolverWithIonicModel(&fr2000_ode_system_mid,
                                   end_time,
                                   "FR2000Mid",
                                   500, false);

        CheckCellModelResults("FR2000Mid");

        CellFaberRudy2000FromCellML fr2000_ode_system_epi(p_solver, p_stimulus);
        fr2000_ode_system_epi.SetParameter("ScaleFactorGks",1.154);
        fr2000_ode_system_epi.SetParameter("ScaleFactorIto",1.0);
        fr2000_ode_system_epi.SetParameter("ScaleFactorGkr",1.0);

        // Solve and write to file
        RunOdeSolverWithIonicModel(&fr2000_ode_system_epi,
                                   end_time,
                                   "FR2000Epi",
                                   500, false);

        CheckCellModelResults("FR2000Epi");
    }


    void TestSolverForFox2002WithRegularStimulus(void)
    {
        clock_t ck_start, ck_end;

        // Set stimulus
        double magnitude = -80.0;
        double duration  = 1.0 ;  // ms
        double start = 50.0; // ms
        double period = 500; // ms
        boost::shared_ptr<RegularStimulus> p_stimulus(new RegularStimulus(magnitude, duration, period, start));

        double end_time = 200.0;  // milliseconds

        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.002, 0.002, 0.002); // 0.005 leads to NaNs.

        boost::shared_ptr<EulerIvpOdeSolver> p_solver(new EulerIvpOdeSolver);
        CellFoxModel2002FromCellML fox_ode_system(p_solver, p_stimulus);

        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.01, 0.01, 0.01);
        CellFoxModel2002FromCellMLBackwardEuler backward_system(p_solver, p_stimulus); // solver ignored

        // Mainly for coverage, and to test consistency of GetIIonic
        TS_ASSERT_DELTA(fox_ode_system.GetIIonic(), backward_system.GetIIonic(), 1e-4);

        // Solve and write to file
        ck_start = clock();
        RunOdeSolverWithIonicModel(&fox_ode_system,
                                   end_time,
                                   "FoxRegularStim",
                                   500);
        ck_end = clock();
        double forward = (double)(ck_end - ck_start)/CLOCKS_PER_SEC;

        CheckCellModelResults("FoxRegularStim");

        // Solve using Backward Euler
        ck_start = clock();
        RunOdeSolverWithIonicModel(&backward_system,
                                   end_time,
                                   "BackwardFoxRegularStim",
                                   100);
        ck_end = clock();
        double backward = (double)(ck_end - ck_start)/CLOCKS_PER_SEC;

        CompareCellModelResults("FoxRegularStim", "BackwardFoxRegularStim", 0.2);

        std::cout << "Run times:\n\tForward: " << forward
                  << "\n\tBackward: " << backward
                  << std::endl;

    }

    void TestSolveForBackwardNoble98WithSimpleStimulus(void)
    {
        clock_t ck_start, ck_end;

        // Set stimulus
        double magnitude_stimulus = -3/0.095;  // uA/cm2
        double duration_stimulus = 3;  // ms
        double start_stimulus = 10.0;   // ms
        boost::shared_ptr<SimpleStimulus> p_stimulus(new SimpleStimulus(magnitude_stimulus,
                                                                        duration_stimulus,
                                                                        start_stimulus));

        // Just adding to check that multi-stim works properly with a cell model.
        boost::shared_ptr<MultiStimulus> p_multi_stim(new MultiStimulus);
        p_multi_stim->AddStimulus(p_stimulus);

        double time_step = 0.2;

        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(time_step, time_step, time_step);

        boost::shared_ptr<AbstractIvpOdeSolver> p_no_solver;
        CellNobleVargheseKohlNoble1998aFromCellMLBackwardEuler n98_backward_system(p_no_solver, p_multi_stim);

        // Solve and write to file
        ck_start = clock();
        RunOdeSolverWithIonicModel(&n98_backward_system,
                                   150.0,
                                   "N98BackwardResult",
                                   1);
        ck_end = clock();
        double backward = (double)(ck_end - ck_start)/CLOCKS_PER_SEC;
        std::cout << "\n\tBackward: " << backward << std::endl;

        CheckCellModelResults("N98BackwardResult");
        //TS_ASSERT_DELTA( n98_backward_system.GetIIonic(), 0.023, 1e-3);
        //This cell now returns a current density
        TS_ASSERT_DELTA( n98_backward_system.GetIIonic(), 0.2462, 1e-3);

        ///\todo compare with the forward results?
    }

    void TestSolveForTT06WithSimpleStimulus(void)
    {

        // This is a shortened test. Longer tests correctly produced AP
        // Full testing for AP in the nightly build
        double simulation_end=40;//end time, in milliseconds for this model

        // Set the stimulus, the following values are appropriate for single cell simulations of this model.
        double magnitude = -38.0;   // pA/pF
        double duration = 1.0;  // ms
        double start = 5;   // ms
        boost::shared_ptr<SimpleStimulus> p_stimulus(new SimpleStimulus(magnitude,
                                                                        duration,
                                                                        start));

        boost::shared_ptr<EulerIvpOdeSolver> p_solver(new EulerIvpOdeSolver); //define the solver
        HeartConfig::Instance()->SetOdeTimeStep(0.001);// with Forward Euler, this must be as small as 0.001.
        CellTenTusscher2006EpiFromCellML TT_model(p_solver, p_stimulus);

        //Default values for the scale factors, other values tested in the nightly build
        TT_model.SetParameter("ScaleFactorIto", 1.0);
        TT_model.SetParameter("ScaleFactorGkr", 1.0);
        TT_model.SetParameter("ScaleFactorGks", 1.0);

        // Solve and write to file
        RunOdeSolverWithIonicModel(&TT_model,
                                   simulation_end,
                                   "TenTusscher",
                                   1000,
                                   true);
        //Check against validated data
        //These data are considered valid after (visually) checking against output from  CellML code of the model for an epicardial cell
        // and also numerically compared against pycml automatically generated code.
        CheckCellModelResults("TenTusscher");

        //Test the GetIIonic method against one hardcoded value.
        TS_ASSERT_DELTA( TT_model.GetIIonic(), -0.1843, 1e-3);

        //Test the GetIIonic method against one hardcoded value for initial values of voltage
        //(mainly for coverage of different if conditions in sodium channel gates for different voltages)
        CellTenTusscher2006EpiFromCellML TT_model_initial(p_solver, p_stimulus);
        TS_ASSERT_DELTA(TT_model_initial.GetIIonic(), 0.0012 , 1e-3);
    }

    void TestBackwardEulerTenTusscher06(void)
    {
        clock_t ck_start, ck_end;
        // Set stimulus
        double magnitude = -1800;
        double duration  = 0.05 ;  // ms
        double when = 5.0; // ms

        boost::shared_ptr<SimpleStimulus> p_stimulus(new SimpleStimulus(magnitude, duration, when));
        double end_time = 50.0;

        HeartConfig::Instance()->SetOdeTimeStep(0.001);

        // Define solver passed in both constructor but used only by forward Euler
        boost::shared_ptr<EulerIvpOdeSolver> p_solver(new EulerIvpOdeSolver);

        // Solve using backward euler
        CellTenTusscher2006EpiFromCellMLBackwardEuler tt06_backward_euler(p_solver, p_stimulus);

        ck_start = clock();
        RunOdeSolverWithIonicModel(&tt06_backward_euler,
                                   end_time,
                                   "TenTusscherBackwardEuler");
        ck_end = clock();
        double backward = (double)(ck_end - ck_start)/CLOCKS_PER_SEC;

        TS_ASSERT_DELTA( tt06_backward_euler.GetIIonic(), -0.0413, 1e-3);

        // Solve using forward euler
        HeartConfig::Instance()->SetOdeTimeStep(0.001);// with Forward Euler, this must be as small as 0.001.

        CellTenTusscher2006EpiFromCellML tt06_ode_system(p_solver, p_stimulus);
        ck_start = clock();
        RunOdeSolverWithIonicModel(&tt06_ode_system,
                                   end_time,
                                   "TenTusscherForward");
        ck_end = clock();
        double forward = (double)(ck_end - ck_start)/CLOCKS_PER_SEC;

        // Compare results
        CompareCellModelResults("TenTusscherForward", "TenTusscherBackwardEuler", 0.03, true);

        std::cout << "Run times:\n\tForward: " << forward << "\n\tBackward: "
          << backward << std::endl;
    }

    void TestDifrancescoNoble1985(void)
    {
        // Set stimulus (no stimulus in this case because this cell is self excitatory)
        boost::shared_ptr<ZeroStimulus> p_stimulus(new ZeroStimulus);

        boost::shared_ptr<EulerIvpOdeSolver> p_solver(new EulerIvpOdeSolver); //define the solver
        HeartConfig::Instance()->SetOdeTimeStep(0.01);
        CellDiFrancescoNoble1985FromCellML purkinje_ode_system(p_solver, p_stimulus);

        // Solve and write to file
        RunOdeSolverWithIonicModel(&purkinje_ode_system,
                                   1800,/*end time, in milliseconds for this model*/
                                   "DiFrancescoNoble",
                                   1000);

        // Check against validated data
        // (the valid data have been checked against CellML code of the model known to be valid).
        CheckCellModelResults("DiFrancescoNoble");

         //Test the GetIIonic method against one hardcoded value.
        TS_ASSERT_DELTA(purkinje_ode_system.GetIIonic(), -0.0141, 1e-3);

        // This model has no stimulus metadata so this method should not exist:
        TS_ASSERT_EQUALS(purkinje_ode_system.HasCellMLDefaultStimulus(),false);
        TS_ASSERT_THROWS_THIS(purkinje_ode_system.UseCellMLDefaultStimulus(),"This class has no default stimulus from CellML metadata.");
    }

    void TestMahajan2008(void)
    {
        // Set stimulus
        double magnitude_stimulus = -1800;
        boost::shared_ptr<RegularStimulus> p_stimulus(new RegularStimulus(magnitude_stimulus,
                                                                          0.05,
                                                                          1000,
                                                                          10.0));

        boost::shared_ptr<EulerIvpOdeSolver> p_solver(new EulerIvpOdeSolver); //define the solver
        HeartConfig::Instance()->SetOdeTimeStep(0.001);
        CellMahajan2008FromCellML rabbit_ode_system(p_solver, p_stimulus);

        //default values for scale factors. They are tested separately in the nightly build.
        rabbit_ode_system.SetParameter("ScaleFactorGks",1.0);
        rabbit_ode_system.SetParameter("ScaleFactorIto",1.0);
        rabbit_ode_system.SetParameter("ScaleFactorGkr",1.0);

        //Test the GetIIonic method against one hardcoded value.
        TS_ASSERT_DELTA(rabbit_ode_system.GetIIonic(), 0.0027, 1e-3);

        // Solve and write to file
        RunOdeSolverWithIonicModel(&rabbit_ode_system,
                                   80,/*end time, in milliseconds for this model*/
                                   "Mahajan2008",
                                   1000);
        // Check against validated data
        // (the code for the Mahajan model was generated from a CellML code known to be valid)
        CheckCellModelResults("Mahajan2008");
    }

    void TestBackwardEulerMahajan(void)
    {
        clock_t ck_start, ck_end;
        // Set stimulus
        double magnitude = -1800;
        double duration  = 0.05 ;  // ms
        double when = 5.0; // ms

        boost::shared_ptr<SimpleStimulus> p_stimulus(new SimpleStimulus(magnitude, duration, when));
        double end_time = 50.0;

        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.005, 0.05, 0.05);

        // Define solver passed in both constructor but used only by forward Euler
        boost::shared_ptr<EulerIvpOdeSolver> p_solver(new EulerIvpOdeSolver);

        // Solve using backward euler
        CellMahajan2008FromCellMLBackwardEuler mahajan_backward_euler(p_solver, p_stimulus);

        ck_start = clock();
        RunOdeSolverWithIonicModel(&mahajan_backward_euler,
                                   end_time,
                                   "MahajanBackwardEuler");
        ck_end = clock();
        double backward = (double)(ck_end - ck_start)/CLOCKS_PER_SEC;

        TS_ASSERT_DELTA( mahajan_backward_euler.GetIIonic(), 0.1244, 1e-3);

        // Solve using forward euler
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.005, 0.05, 0.05);

        CellMahajan2008FromCellML mahajan_ode_system(p_solver, p_stimulus);
        ck_start = clock();
        RunOdeSolverWithIonicModel(&mahajan_ode_system,
                                   end_time,
                                   "MahajanForward");
        ck_end = clock();
        double forward = (double)(ck_end - ck_start)/CLOCKS_PER_SEC;

        // Compare results
        CompareCellModelResults("MahajanForward", "MahajanBackwardEuler", 0.03, true);

        std::cout << "Run times:\n\tForward: " << forward << "\n\tBackward: "
                  << backward << std::endl;
    }

    void TestMaleckar(void)
    {
        boost::shared_ptr<RegularStimulus> p_stimulus(new RegularStimulus(-5.6, // Changed because now it is in the right units.
                                                                          6,
                                                                          1000,
                                                                          4.0));

        boost::shared_ptr<EulerIvpOdeSolver> p_solver(new EulerIvpOdeSolver); //define the solver
        HeartConfig::Instance()->SetOdeTimeStep(0.001);
        CellMaleckar2008FromCellML atrial_ode_system(p_solver, p_stimulus);

        //default values, other values tested in the nightly build
        atrial_ode_system.SetParameter("ScaleFactorGks",1.0);
        atrial_ode_system.SetParameter("ScaleFactorIto",1.0);
        atrial_ode_system.SetParameter("ScaleFactorGkr",1.0);
        atrial_ode_system.SetParameter("ScaleFactorGna",1.0);
        atrial_ode_system.SetParameter("ScaleFactorAch",1e-24);
        atrial_ode_system.SetParameter("ScaleFactorGNaK",1.0);
        atrial_ode_system.SetParameter("ScaleFactorGNaCa",1.0);
        atrial_ode_system.SetParameter("ScaleFactorGKur",1.0);
        atrial_ode_system.SetParameter("ScaleFactorGK1",1.0);
        atrial_ode_system.SetParameter("ScaleFactorGCaL",1.0);
        atrial_ode_system.SetParameter("ScaleFactorAZD",0.0);

        // Solve and write to file for a short time
        RunOdeSolverWithIonicModel(&atrial_ode_system,
                                   25,/*end time*/
                                   "Maleckar2009",
                                   1000);

        // Check against validated data(The full AP is attached to ticket 1194)
        CheckCellModelResults("Maleckar2009");

        //Test the GetIIonic method against one hardcoded value.
        TS_ASSERT_DELTA(atrial_ode_system.GetIIonic(), 1.4426, 1e-3);

    }
//    Uncomment the includes for the models too
//
//    void TestSolverForN98WithSimpleStimulus(void)
//    {
//        clock_t ck_start, ck_end;
//
//        // Set stimulus
//        double magnitude_stimulus = 0.0;  // uA/cm2
//        double duration_stimulus = 0.5;  // ms
//        double start_stimulus = 10.0;   // ms
//        SimpleStimulus stimulus(magnitude_stimulus,
//                                 duration_stimulus,
//                                 start_stimulus);
//
//        // Solve forward
//        HeartConfig::Instance()->SetOdeTimeStep(0.0005);
//        boost::shared_ptr<EulerIvpOdeSolver> p_solver(new EulerIvpOdeSolver);
//        CellNobleVargheseKohlNoble1998aFromCellMLOpt n98_ode_system(p_solver, p_stimulus);
//
//        std::vector<double> dY(22);
//
//        n98_ode_system.EvaluateYDerivatives(0, n98_ode_system.rGetStateVariables(), dY);
//
//        for (unsigned i = 0; i < dY.size(); ++i)
//        {
//            std::cout << dY[i] << "\n";
//        }
//
//        ck_start = clock();
//        RunOdeSolverWithIonicModel(&n98_ode_system,
//                                   150.0,
//                                   "N98RegResult");
//        ck_end = clock();
//        double forward = (double)(ck_end - ck_start)/CLOCKS_PER_SEC;
//
//        TS_ASSERT( !std::isnan( n98_ode_system.GetIIonic() ) );
//
//        // Solve backward
//        HeartConfig::Instance()->SetOdeTimeStep(0.1);
//        CML_noble_varghese_kohl_noble_1998_basic_pe_lut_be n98_be_ode_system(p_stimulus);
//
//        ck_start = clock();
//        RunOdeSolverWithIonicModel(&n98_be_ode_system,
//                                   150.0,
//                                   "BeN98RegResult");
//        ck_end = clock();
//        double backward = (double)(ck_end - ck_start)/CLOCKS_PER_SEC;
//
//        TS_ASSERT( !std::isnan( n98_be_ode_system.GetIIonic() ) );
//
//        CompareCellModelResults("N98RegResult", "BeN98RegResult", 0.2);
//
//        std::cout << "Run times:\n\tForward: " << forward
//                  << "\n\tBackward: " << backward
//                  << std::endl;
//
//
//    }

    void TestLR1991AndN98WithSacArchiving(void)
    {
        //Archive
        OutputFileHandler handler("archive", false);
        handler.SetArchiveDirectory();
        std::string archive_filename = ArchiveLocationInfo::GetProcessUniqueFilePath("lr91.arch");

        // Save
        {
            // Set stimulus
            double magnitude_stimulus = -3;  // uA/cm2
            double duration_stimulus = 3;  // ms
            double start_stimulus = 10.0;   // ms
            boost::shared_ptr<SimpleStimulus> p_stimulus(new SimpleStimulus(magnitude_stimulus,
                                                                            duration_stimulus,
                                                                            start_stimulus));
            boost::shared_ptr<EulerIvpOdeSolver> p_solver(new EulerIvpOdeSolver);
            double time_step = 0.01;

            HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(time_step, time_step, time_step);

            // Check Standard
            AbstractCardiacCellInterface* const p_luo_rudy_cell = new CellLuoRudy1991FromCellML(p_solver, p_stimulus);
            AbstractCardiacCellInterface* const p_n98_with_sac = new CML_noble_varghese_kohl_noble_1998_basic_with_sac(p_solver, p_stimulus);

            p_n98_with_sac->SetStretch(1.1);

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            output_arch <<  p_luo_rudy_cell;
            output_arch <<  p_n98_with_sac;

            // These results are in the repository and should be replicated after the load below
//            RunOdeSolverWithIonicModel(p_luo_rudy_cell,
//                           50.0,
//                           "LRAfterArchiveValidData");

            delete p_luo_rudy_cell;
            delete p_n98_with_sac;
        }
        // Load
        {
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            AbstractCardiacCellInterface* p_luo_rudy_cell;
            input_arch >> p_luo_rudy_cell;

            AbstractCardiacCellInterface* p_n98_with_sac;
            input_arch >> p_n98_with_sac;

            TS_ASSERT_EQUALS( p_luo_rudy_cell->GetNumberOfStateVariables(), 8U );
            TS_ASSERT_EQUALS( p_n98_with_sac->GetNumberOfStateVariables(), 22U );

            CML_noble_varghese_kohl_noble_1998_basic_with_sac*   p_n98_with_sac_conc
               = dynamic_cast<CML_noble_varghese_kohl_noble_1998_basic_with_sac*>(p_n98_with_sac);
            TS_ASSERT_DELTA( p_n98_with_sac_conc->GetStretch(), 1.1, 1e-5 );

            RunOdeSolverWithIonicModel(p_luo_rudy_cell,
                                       50.0,
                                       "LRAfterArchive");

            CheckCellModelResults("LRAfterArchive");

            delete p_luo_rudy_cell;
            delete p_n98_with_sac;
        }
    }

    void TestMaleckar2009Archiving(void)
    {
        //Archive
        OutputFileHandler handler("archive", false);
        handler.SetArchiveDirectory();
        std::string archive_filename = ArchiveLocationInfo::GetProcessUniqueFilePath("Maleckar.arch");

        // Save
        {
            // Set stimulus
            boost::shared_ptr<RegularStimulus> p_stimulus(new RegularStimulus(-5.6, // Now in consistent Chaste units.
                                                                          6,
                                                                          1000,
                                                                          4.0));
            boost::shared_ptr<EulerIvpOdeSolver> p_solver(new EulerIvpOdeSolver);
            double time_step = 0.01;

            HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(time_step, time_step, time_step);

            // Check Standard
            AbstractCardiacCellInterface* const p_maleckar_cell = new CellMaleckar2008FromCellML(p_solver, p_stimulus);

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            output_arch <<  p_maleckar_cell;

             //These results have been copied in the repository
             // after running this line and should be replicated after the load below
//            RunOdeSolverWithIonicModel(p_maleckar_cell,
//                           20.0,
//                           "MaleckarAfterArchiveValidData");

            delete p_maleckar_cell;
        }
        // Load
        {
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            AbstractCardiacCellInterface* p_maleckar_cell;
            input_arch >> p_maleckar_cell;

            TS_ASSERT_EQUALS( p_maleckar_cell->GetNumberOfStateVariables(), 30U );


            RunOdeSolverWithIonicModel(p_maleckar_cell,
                                       20.0,
                                       "MaleckarAfterArchive");

            CheckCellModelResults("MaleckarAfterArchive");

            delete p_maleckar_cell;
        }
    }

    void TestBackwardCellsArchiving(void)
    {
        //Archive
        OutputFileHandler handler("archive", false);
        handler.SetArchiveDirectory();
        std::string archive_filename = ArchiveLocationInfo::GetProcessUniqueFilePath("backward_cells.arch");

        double time_step = 0.01;

        // Save
        {
            // Set stimulus
            double magnitude_stimulus = -3;  // uA/cm2
            double magnitude_stimulus_noble = magnitude_stimulus/0.095;  // uA/cm2
            double duration_stimulus = 3;  // ms
            double start_stimulus = 10.0;   // ms

            boost::shared_ptr<SimpleStimulus> p_stimulus(new SimpleStimulus(magnitude_stimulus,
                                                                            duration_stimulus,
                                                                            start_stimulus));

            boost::shared_ptr<SimpleStimulus> p_noble_stimulus(new SimpleStimulus(magnitude_stimulus_noble,
                                                                                  duration_stimulus,
                                                                                  start_stimulus));

            HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(time_step, time_step, time_step);

            boost::shared_ptr<EulerIvpOdeSolver> p_solver(new EulerIvpOdeSolver);
            AbstractCardiacCellInterface* const p_backward_cell1 = new CellLuoRudy1991FromCellMLBackwardEuler(p_solver, p_stimulus);
            AbstractCardiacCellInterface* const p_backward_cell2 = new CellFoxModel2002FromCellMLBackwardEuler(p_solver, p_stimulus);
            AbstractCardiacCellInterface* const p_backward_cell3 = new CellNobleVargheseKohlNoble1998aFromCellMLBackwardEuler(p_solver, p_noble_stimulus);
            AbstractCardiacCellInterface* const p_backward_cell4 = new CellMahajan2008FromCellMLBackwardEuler(p_solver, p_stimulus);
            AbstractCardiacCellInterface* const p_backward_cell5 = new CellTenTusscher2006EpiFromCellMLBackwardEuler(p_solver, p_stimulus);

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            output_arch <<  p_backward_cell1;
            output_arch <<  p_backward_cell2;
            output_arch <<  p_backward_cell3;
            output_arch <<  p_backward_cell4;
            output_arch <<  p_backward_cell5;

            // These results are in the repository and should be replicated after the load below
//            RunOdeSolverWithIonicModel(p_backward_cell1,
//                                       50.0,
//                                       "Backward1AfterArchiveValidData");
//
//            RunOdeSolverWithIonicModel(p_backward_cell2,
//                                       50.0,
//                                       "Backward2AfterArchiveValidData");
//
//            RunOdeSolverWithIonicModel(p_backward_cell3,
//                                       50.0,
//                                       "Backward3AfterArchiveValidData");
            delete p_backward_cell1;
            delete p_backward_cell2;
            delete p_backward_cell3;
            delete p_backward_cell4;
            delete p_backward_cell5;
        }
        // Load
        {
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            AbstractCardiacCellInterface* p_backward_cell1;
            AbstractCardiacCellInterface* p_backward_cell2;
            AbstractCardiacCellInterface* p_backward_cell3;
            AbstractCardiacCellInterface* p_backward_cell4;
            AbstractCardiacCellInterface* p_backward_cell5;
            input_arch >> p_backward_cell1;
            input_arch >> p_backward_cell2;
            input_arch >> p_backward_cell3;
            input_arch >> p_backward_cell4;
            input_arch >> p_backward_cell5;

            TS_ASSERT_EQUALS( p_backward_cell1->GetNumberOfStateVariables(), 8U );
            TS_ASSERT_EQUALS( p_backward_cell2->GetNumberOfStateVariables(), 13U );
            TS_ASSERT_EQUALS( p_backward_cell3->GetNumberOfStateVariables(), 22U );
            TS_ASSERT_EQUALS( p_backward_cell4->GetNumberOfStateVariables(), 26U );
            TS_ASSERT_EQUALS( p_backward_cell5->GetNumberOfStateVariables(), 19U );

            RunOdeSolverWithIonicModel(p_backward_cell1,
                                       50.0,
                                       "Backward1AfterArchive");

            RunOdeSolverWithIonicModel(p_backward_cell2,
                                       50.0,
                                       "Backward2AfterArchive");

            RunOdeSolverWithIonicModel(p_backward_cell3,
                                       50.0,
                                       "Backward3AfterArchive");

            RunOdeSolverWithIonicModel(p_backward_cell4,
                                       50.0,
                                       "Backward4AfterArchive");

            RunOdeSolverWithIonicModel(p_backward_cell5,
                                       50.0,
                                       "Backward5AfterArchive");


            CheckCellModelResults("Backward1AfterArchive", "", 2e-3);
            CheckCellModelResults("Backward2AfterArchive");
            CheckCellModelResults("Backward3AfterArchive");
            CheckCellModelResults("Backward4AfterArchive");
            CheckCellModelResults("Backward5AfterArchive");

            delete p_backward_cell1;
            delete p_backward_cell2;
            delete p_backward_cell3;
            delete p_backward_cell4;
            delete p_backward_cell5;
        }
     }

    void TestPyCMLArchiving(void)
    {
        //Archive
        OutputFileHandler handler("archive", false);
        handler.SetArchiveDirectory();
        std::string archive_filename = ArchiveLocationInfo::GetProcessUniqueFilePath("noble98.arch");

        // Save
        {
            // Set stimulus
            double magnitude_stimulus = -3/0.095;  // uA/cm2
            double duration_stimulus = 3;  // ms
            double start_stimulus = 10.0;   // ms
            boost::shared_ptr<SimpleStimulus> p_stimulus(new SimpleStimulus(magnitude_stimulus,
                                                                            duration_stimulus,
                                                                            start_stimulus));
            boost::shared_ptr<EulerIvpOdeSolver> p_solver(new EulerIvpOdeSolver);
            double time_step = 0.01;

            HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(time_step, time_step, time_step);

            // Check Standard
            AbstractCardiacCellInterface* const p_n98_cell = new CellNobleVargheseKohlNoble1998aFromCellML(p_solver, p_stimulus);
            // and SAC
            AbstractCardiacCellInterface* const p_n98_sac = new CML_noble_varghese_kohl_noble_1998_basic_with_sac(p_solver, p_stimulus);
            // and "0d" backward Euler
            AbstractCardiacCellInterface* const p_hh52_be = new CellHodgkinHuxley1952FromCellMLBackwardEuler(p_solver, p_stimulus);

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            output_arch << p_n98_cell;
            output_arch << p_n98_sac;
            output_arch << p_hh52_be;

            // These results are in the repository and should be replicated after the load below
//            RunOdeSolverWithIonicModel(p_n98_cell,
//                                       50.0,
//                                       "N98AfterArchiveValidData");

            delete p_n98_cell;
            delete p_n98_sac;
            delete p_hh52_be;
        }
        // Load
        {
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            AbstractCardiacCellInterface* p_n98_cell;
            input_arch >> p_n98_cell;

            TS_ASSERT_EQUALS( p_n98_cell->GetNumberOfStateVariables(), 22U );
            TS_ASSERT(dynamic_cast<CellNobleVargheseKohlNoble1998aFromCellML*>(p_n98_cell));

            RunOdeSolverWithIonicModel(p_n98_cell,
                                       50.0,
                                       "N98AfterArchive");

            CheckCellModelResults("N98AfterArchive");
            AbstractCardiacCellInterface* p_n98_cell_sac;
            input_arch >> p_n98_cell_sac;
            TS_ASSERT_EQUALS( p_n98_cell_sac->GetNumberOfStateVariables(), 22U );
            TS_ASSERT(dynamic_cast<CML_noble_varghese_kohl_noble_1998_basic_with_sac*>(p_n98_cell_sac));

            AbstractCardiacCellInterface* p_hh52_be;
            input_arch >> p_hh52_be;
            TS_ASSERT_EQUALS(p_hh52_be->GetNumberOfStateVariables(), 4u);
            TS_ASSERT(dynamic_cast<CellHodgkinHuxley1952FromCellMLBackwardEuler*>(p_hh52_be));

            delete p_n98_cell;
            delete p_n98_cell_sac;
            delete p_hh52_be;
        }
     }

    void TestBackwardEulerDifficultCase()
    {

        //These data come from a failing human heart mesh test, but have been rounded
        double dodgy_state_vars_array[19] = {-6.15475,0.00808679,0.284434,0.00633525,0.994096,0.0321343,0.402544,0.730188,0.856068,0.959682,0.998295,0.912007,0.02408,0.000115671,3.63196,0.00175538,0.937932,8.62141,136.891};
        std::vector<double> dodgy_state_vars(dodgy_state_vars_array,dodgy_state_vars_array+19);

        double step=0.1; //Time step in full 3D test was 0.05


        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(step, step, step);

        boost::shared_ptr<ZeroStimulus> p_stimulus(new ZeroStimulus);
        boost::shared_ptr<EulerIvpOdeSolver> p_solver(new EulerIvpOdeSolver);
        CellTenTusscher2006EpiFromCellMLBackwardEuler tt06_backward_euler(p_solver, p_stimulus);

        tt06_backward_euler.rGetStateVariables() = dodgy_state_vars;
        tt06_backward_euler.ComputeExceptVoltage(0.0, 3*step);
    }


private:
    void TryTestLr91WithVoltageDrop(unsigned ratio) //
    {
        double end_time = 10;        // ms
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.01/ratio, 0.01, 0.01);

        boost::shared_ptr<ZeroStimulus> p_zero_stimulus(new ZeroStimulus);
        boost::shared_ptr<EulerIvpOdeSolver> p_solver(new EulerIvpOdeSolver);
        CellLuoRudy1991FromCellML lr91_ode_system(p_solver, p_zero_stimulus);
        double time=0.0;
        double start_voltage=-83.853;
        double end_voltage=-100;
        while (time<end_time)
        {
            double next_time=time + HeartConfig::Instance()->GetPdeTimeStep();
            lr91_ode_system.SetVoltage( start_voltage + (end_voltage-start_voltage)*time/end_time );
            lr91_ode_system.ComputeExceptVoltage(time, next_time);
            time=next_time;
        }
    }
};


#endif //_TESTIONICMODELS_HPP_
