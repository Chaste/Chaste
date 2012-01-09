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


#ifndef TESTNHSMODELWITHBACKWARDSOLVER_HPP_
#define TESTNHSMODELWITHBACKWARDSOLVER_HPP_

#include <cxxtest/TestSuite.h>

#include "NhsModelWithBackwardSolver.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "LuoRudy1991.hpp"
#include "ZeroStimulus.hpp"
#include "SimpleDataWriter.hpp"
#include "Timer.hpp"



class TestNhsModelWithBackwardSolver : public CxxTest::TestSuite
{
private:
    double GetSampleCaIValue()
    {
        boost::shared_ptr<EulerIvpOdeSolver> p_euler_solver(new EulerIvpOdeSolver);
        boost::shared_ptr<ZeroStimulus> p_zero_stimulus(new ZeroStimulus);
        CellLuoRudy1991FromCellML lr91(p_euler_solver, p_zero_stimulus);
        return lr91.rGetStateVariables()[lr91.GetStateVariableIndex("cytosolic_calcium_concentration")];
    }

public:
    void TestSolverSingleTimestep()
    {
        NhsModelWithBackwardSolver nhs_backward;

        // lam=const, dlamdt not zero doesn't make much sense, just for testing purposes
        nhs_backward.SetStretchAndStretchRate(0.5, 0.1);
        double Ca_I = GetSampleCaIValue();
        nhs_backward.SetIntracellularCalciumConcentration(Ca_I);

        // solve system (but don't update state vars yet)
        nhs_backward.RunDoNotUpdate(0, 0.1, 0.1); // one timestep

        NhsContractionModel nhs_forward;

        unsigned num_vars = nhs_backward.GetNumberOfStateVariables();
        for(unsigned i=0; i<num_vars; i++)
        {
            // both should be the same (ie initial values)
            TS_ASSERT_DELTA(nhs_backward.rGetStateVariables()[i],
                            nhs_forward.rGetStateVariables()[i],
                            1e-12);
        }

        // solve system with euler
        nhs_forward.SetStretchAndStretchRate(0.5, 0.1);
        nhs_forward.SetIntracellularCalciumConcentration(Ca_I);
        EulerIvpOdeSolver euler_solver;
        euler_solver.SolveAndUpdateStateVariable(&nhs_forward, 0, 0.1, 0.1);  // one timestep

        // update state vars on implicit system
        nhs_backward.UpdateStateVariables();

        for(unsigned i=0; i<num_vars; i++)
        {
            std::cout << nhs_backward.rGetStateVariables()[i] << " "
                      << nhs_forward.rGetStateVariables()[i] << "\n";

            // we want these within 10% of each other. Note we expect the implicit
            // solver to be more accurate than the explicit solver, and the timestep
            // is quite large (as we want non-zero solutions), so can't expect them
            // to be too close.
            TS_ASSERT_DELTA(nhs_backward.rGetStateVariables()[i],
                            nhs_forward.rGetStateVariables()[i],
                            std::max(fabs(nhs_backward.rGetStateVariables()[i]*1e-1),1e-10));
        }
    }


    void TestSolverManyTimestepsCompareWithEuler()
    {
        clock_t ck_start, ck_end;

        NhsModelWithBackwardSolver nhs_backward;

        // lam=const, dlamdt not zero doesn't make much sense, just for testing purposes
        nhs_backward.SetStretchAndStretchRate(0.5, 0.1);
        double Ca_I = GetSampleCaIValue();
        // bigger Ca_I, so we get some active tension (and so the a few iterations are
        // needed when solving for T_a and z
        nhs_backward.SetIntracellularCalciumConcentration(10*Ca_I);

        // solve system and update
        ck_start = clock();
        nhs_backward.RunAndUpdate(0, 100, 0.01);
        ck_end = clock();
        double implicit_solve_time = (double)(ck_end - ck_start)/CLOCKS_PER_SEC;

        //nhs_backward.UpdateStateVariables();

        // GetNextActiveTension should now be equal to baseclass::GetActiveTension(),
        // as the state vars have been updated
        TS_ASSERT_DELTA(nhs_backward.GetNextActiveTension(),
                        nhs_backward.GetActiveTension(),
                        1e-12);

        // solve system with euler
        NhsContractionModel nhs_forward;
        nhs_forward.SetStretchAndStretchRate(0.5, 0.1);
        nhs_forward.SetIntracellularCalciumConcentration(10*Ca_I);

        ck_start = clock();
        nhs_forward.RunAndUpdate(0, 100, 0.01);
        ck_end = clock();
        double explicit_solve_time = (double)(ck_end - ck_start)/CLOCKS_PER_SEC;

        unsigned num_vars = nhs_backward.GetNumberOfStateVariables();
        for(unsigned i=0; i<num_vars; i++)
        {
            std::cout << nhs_backward.rGetStateVariables()[i] << " "
                      << nhs_forward.rGetStateVariables()[i] << "\n";

            // small timesteps, want these to be very close (1%). and they are
            TS_ASSERT_DELTA(nhs_backward.rGetStateVariables()[i],
                            nhs_forward.rGetStateVariables()[i],
                            fabs(nhs_backward.rGetStateVariables()[i]*1e-2));
        }
        std::cout << "TIMES: " << implicit_solve_time << " " << explicit_solve_time << "\n\n";
    }

//// test how large a timestep the implicit solver can get away with. needs more study
//    void TestImplicitSolverWithLargeTimeSteps()
//    {
//        NhsModelWithBackwardSolver nhs_backward;
//        nhs_backward.SetStretchAndStretchRate(0.5, 0.1);
//        nhs_backward.SetIntracellularCalciumConcentration(10*GetSampleCaIValue());
//
//        nhs_backward.RunDoNotUpdate(0, 100, 0.01);
//        nhs_backward.UpdateStateVariables();
//
//        NhsModelWithBackwardSolver nhs_backward2;
//        nhs_backward2.SetStretchAndStretchRate(0.5, 0.1);
//        nhs_backward2.SetIntracellularCalciumConcentration(10*GetSampleCaIValue());
//
//        nhs_backward2.RunDoNotUpdate(0, 100, 1);
//        nhs_backward2.UpdateStateVariables();
//
//        unsigned num_vars = nhs_backward.GetNumberOfStateVariables();
//        for(unsigned i=0; i<num_vars; i++)
//        {
//            std::cout << nhs_backward.rGetStateVariables()[i] << " "
//                      << nhs_backward2.rGetStateVariables()[i] << "\n";
//
//            TS_ASSERT_DELTA(nhs_backward.rGetStateVariables()[i],
//                            nhs_backward2.rGetStateVariables()[i],
//                            fabs(nhs_backward.rGetStateVariables()[i]*1e-2));
//        }
//    }

    // test that checks RunDoNotUpdate does not do anything permanent on the class,
    // by checking by doing things repeatedly, and changing the order, makes no
    // difference
    void TestRunDoesNotUpdate()
    {
        NhsModelWithBackwardSolver system;

        double Ca_I = GetSampleCaIValue();
        system.SetIntracellularCalciumConcentration(Ca_I);

        // get initial active tension
        double init_Ta = system.GetActiveTension();

        system.SetStretchAndStretchRate(0.6, 0.1);
        system.RunDoNotUpdate(0, 1, 0.01);

        double Ta1 = system.GetNextActiveTension();

        system.SetStretchAndStretchRate(0.6, 0.2);
        system.RunDoNotUpdate(0, 1, 0.01);

        double Ta2 = system.GetNextActiveTension();

        // note that lam/end time etc must be large enough for there
        // to be non-zero Ta at the next time
        TS_ASSERT_DIFFERS(init_Ta, Ta1);
        TS_ASSERT_DIFFERS(init_Ta, Ta2);

        system.SetStretchAndStretchRate(0.6, 0.2);
        system.RunDoNotUpdate(0, 1, 0.01);

        double should_be_Ta2 = system.GetNextActiveTension();

        system.SetStretchAndStretchRate(0.6, 0.1);
        system.RunDoNotUpdate(0, 1, 0.01);

        double should_be_Ta1 = system.GetNextActiveTension();

        TS_ASSERT_EQUALS(Ta1, should_be_Ta1);
        TS_ASSERT_EQUALS(Ta2, should_be_Ta2);

        system.SetStretchAndStretchRate(0.6, 0.1);
        system.RunDoNotUpdate(0, 1, 0.01);

        double should_also_be_Ta1 = system.GetNextActiveTension();

        system.SetStretchAndStretchRate(0.6, 0.2);
        system.RunDoNotUpdate(0, 1, 0.01);

        double should_also_be_Ta2 = system.GetNextActiveTension();

        TS_ASSERT_EQUALS(Ta1, should_also_be_Ta1);
        TS_ASSERT_EQUALS(Ta2, should_also_be_Ta2);
    }

    void TestGetActiveTension()
    {
        NhsModelWithBackwardSolver system;

        double Ca_I = GetSampleCaIValue();
        system.SetIntracellularCalciumConcentration(Ca_I);
        system.SetStretchAndStretchRate(0.6, 0.1);
        system.RunDoNotUpdate(0, 1, 0.01);

        double Ta_at_next_time_before_update = system.GetNextActiveTension();

        system.UpdateStateVariables();

        double Ta = system.GetActiveTension();

        TS_ASSERT_DELTA(Ta, Ta_at_next_time_before_update, 1e-12);
    }
};

#endif /*TESTNHSMODELWITHBACKWARDSOLVER_HPP_*/
