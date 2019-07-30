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

#ifndef TESTVANLEEUWEN2009WNTSWATCELLCYCLEODESYSTEM_HPP_
#define TESTVANLEEUWEN2009WNTSWATCELLCYCLEODESYSTEM_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

// Must be included before any other cell_based or crypt headers
#include "CellBasedSimulationArchiver.hpp"

#include <cstdio>

#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "CryptSimulation2d.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "SloughingCellKiller.hpp"
#include "CryptCellsGenerator.hpp"
#include "VanLeeuwen2009WntSwatCellCycleModelHypothesisOne.hpp"
#include "VanLeeuwen2009WntSwatCellCycleModelHypothesisTwo.hpp"
#include "CylindricalHoneycombMeshGenerator.hpp"
#include "RungeKutta4IvpOdeSolver.hpp"
#include "RungeKuttaFehlbergIvpOdeSolver.hpp"
#include "BackwardEulerIvpOdeSolver.hpp"
#include "ApcOneHitCellMutationState.hpp"
#include "ApcTwoHitCellMutationState.hpp"
#include "BetaCateninOneHitCellMutationState.hpp"
#include "WildTypeCellMutationState.hpp"
#include "SmartPointers.hpp"
#include "CellMutationStatesCountWriter.hpp"
#include "FakePetscSetup.hpp"

class TestVanLeeuwen2009WntSwatCellCycleOdeSystem : public CxxTest::TestSuite
{
public:

    void TestVanLeeuwen2009WntSwatCellCycleEquationsHighWnt()
    {
        double wnt_level = 1.0;

        TS_ASSERT_THROWS_THIS(VanLeeuwen2009WntSwatCellCycleOdeSystem wnt_cell_cycle_system(3,wnt_level),
                              "You must set up this cell cycle ODE system with hypothesis one or two.");

        VanLeeuwen2009WntSwatCellCycleOdeSystem wnt_cell_cycle_system(1,wnt_level);

        double time = 0.0;
        std::vector<double> initial_conditions = wnt_cell_cycle_system.GetInitialConditions();
        TS_ASSERT_EQUALS(initial_conditions.size(), 22u);
        std::vector<double> derivs(initial_conditions.size());

        // Test initial conditions are being set up properly
        TS_ASSERT_DELTA(initial_conditions[0], 7.357000000000000e-01, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[1], 1.713000000000000e-01, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[2], 6.900000000000001e-02, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[3], 3.333333333333334e-03, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[4], 1.000000000000000e-04, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[5], 1.428571428571428e-01, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[6], 2.857142857142857e-02, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[7], 2.120643654085212e-01, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[8], 1.439678172957394e+01, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[9], 0, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[10], 0, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[11], 0, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[12], 1.000000000000000e+01, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[13], 1.028341552112424e+02, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[14], 0, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[15], 2.500000000000000e+01, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[16], 1.439678172957394e+01, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[17], 0, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[18], 0, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[19], 0, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[20], 2.235636835087720e+00, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[21], 1.000000000000000e+00, 1e-7);

        boost::shared_ptr<AbstractCellMutationState> p_state(new WildTypeCellMutationState);
        wnt_cell_cycle_system.SetMutationState(p_state);    // for coverage
        wnt_cell_cycle_system.EvaluateYDerivatives(time, initial_conditions, derivs);

        // Test derivatives are correct at t=0 for these initial conditions

        // Swat's
        TS_ASSERT_DELTA(derivs[0],-1.586627673253325e-02, 1e-5);
        TS_ASSERT_DELTA(derivs[1],-5.532201118824132e-05, 1e-5);
        TS_ASSERT_DELTA(derivs[2],5.676110199115955e+00, 1e-5); // changes due to new beta-catenin levels
        TS_ASSERT_DELTA(derivs[3],-7.449833887043188e-03, 1e-5);
        TS_ASSERT_DELTA(derivs[4],1.549680000000000e-02, 1e-5);

        // Inge's
        for (unsigned i=5; i<initial_conditions.size(); i++)
        {
            TS_ASSERT_DELTA(derivs[i],0.0, 1e-9);
        }
    }

    void TestVanLeeuwen2009WntSwatCellCycleEquationsHighWntHypothesisTwo()
    {
        double wnt_level = 1.0;
        VanLeeuwen2009WntSwatCellCycleOdeSystem wnt_cell_cycle_system(2,wnt_level);

        double time = 0.0;

        std::vector<double> initial_conditions = wnt_cell_cycle_system.GetInitialConditions();
        TS_ASSERT_EQUALS(initial_conditions.size(), 22u);
        std::vector<double> derivs(initial_conditions.size());

        // Test initial conditions are being set up properly
        TS_ASSERT_DELTA(initial_conditions[0], 7.357000000000000e-01, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[1], 1.713000000000000e-01, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[2], 6.900000000000001e-02, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[3], 3.333333333333334e-03, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[4], 1.000000000000000e-04, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[5], 1.428571428571428e-01, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[6], 2.857142857142857e-02, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[7], 2.120643654085212e-01, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[8], 9.391557285347905e-01, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[9], 1.345762600103915e+01, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[10], 0, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[11], 0, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[12], 1.000000000000000e+01, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[13], 6.708255203819932e+00, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[14], 0, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[15], 2.500000000000000e+01, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[16], 9.391557285347906e-01, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[17], 1.345762600103915e+01, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[18], 0, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[19], 0, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[20], 2.235636835087720e+00, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[21], 1.000000000000000e+00, 1e-7);

        wnt_cell_cycle_system.EvaluateYDerivatives(time, initial_conditions, derivs);

        // Test derivatives are correct at t=0 for these initial conditions

        // Swat's
        TS_ASSERT_DELTA(derivs[0],-1.586627673253325e-02, 1e-5);
        TS_ASSERT_DELTA(derivs[1],-5.532201118824132e-05, 1e-5);
        TS_ASSERT_DELTA(derivs[2],5.676110199115955e+00, 1e-5); // changes due to new beta-catenin levels
        TS_ASSERT_DELTA(derivs[3],-7.449833887043188e-03, 1e-5);
        TS_ASSERT_DELTA(derivs[4],1.549680000000000e-02, 1e-5);

        // Inge's
        for (unsigned i=5; i<derivs.size(); i++)
        {
            TS_ASSERT_DELTA(derivs[i],0.0, 1e-9);
        }
    }

    /**
     * And the same for a healthy cell at a low Wnt level
     */
    void TestVanLeeuwen2009WntSwatCellCycleEquationsLowWnt()
    {
        double time = 0.0;
        double wnt_level = 0.0;
        boost::shared_ptr<AbstractCellMutationState> p_wt_state(new WildTypeCellMutationState);
        VanLeeuwen2009WntSwatCellCycleOdeSystem wnt_cell_cycle_system2(1, wnt_level, p_wt_state);
        std::vector<double> initial_conditions = wnt_cell_cycle_system2.GetInitialConditions();
        std::vector<double> derivs(initial_conditions.size());

        // Test initial conditions are being set up properly
        TS_ASSERT_DELTA(initial_conditions[0], 7.357000000000000e-01, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[1], 1.713000000000000e-01, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[2], 6.900000000000001e-02, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[3], 3.333333333333334e-03, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[4], 1.000000000000000e-04, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[5], 6.666666666666666e-01, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[6], 6.666666666666667e-02, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[7], 4.491941767913325e-01, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[8], 2.540291160433374e+00, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[9], -4.440892098500626e-16, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[10], 0, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[11], 0, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[12], 1.000000000000000e+01, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[13], 1.814493686023838e+01, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[14], 0, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[15], 2.500000000000000e+01, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[16], 2.540291160433374e+00, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[17], 0, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[18], 0, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[19], 0, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[20], 4.834939252004746e-01, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[21], 0.000000000000000e+00, 1e-7);

        wnt_cell_cycle_system2.EvaluateYDerivatives(time, initial_conditions, derivs);

        // Test derivatives are correct at t=0 for these initial conditions
        // (figures from Matlab code)
        TS_ASSERT_DELTA(derivs[0],-1.586627673253325e-02, 1e-5);
        TS_ASSERT_DELTA(derivs[1],-5.532201118824132e-05, 1e-5);
        TS_ASSERT_DELTA(derivs[2],9.335139714597300e-01, 1e-5);
        TS_ASSERT_DELTA(derivs[3],-7.449833887043188e-03, 1e-5);
        TS_ASSERT_DELTA(derivs[4],1.549680000000000e-02, 1e-5);

        for (unsigned i=5; i<initial_conditions.size(); i++)
        {
            TS_ASSERT_DELTA(derivs[i],0.0, 1e-9);
        }
    }

    /**
     * A test for the case mutation = 1
     * (An APC +/- mutation)
     */
    void TestVanLeeuwen2009WntSwatCellCycleEquationsAPCOneHit()
    {
        double time = 0.0;
        double wnt_level = 0.5;

        boost::shared_ptr<AbstractCellMutationState> p_apc1(new ApcOneHitCellMutationState);

        VanLeeuwen2009WntSwatCellCycleOdeSystem wnt_cell_cycle_system3(1, wnt_level, p_apc1);

        std::vector<double> initial_conditions = wnt_cell_cycle_system3.GetInitialConditions();

        std::vector<double> derivs(initial_conditions.size());

        TS_ASSERT_DELTA(initial_conditions[0], 7.357000000000000e-01, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[1], 1.713000000000000e-01, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[2], 6.900000000000001e-02, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[3], 3.333333333333334e-03, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[4], 1.000000000000000e-04, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[5], 1.481481481481481e-01, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[6], 4.444444444444445e-02, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[7], 2.186081417430197e-01, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[8], 1.406959291284901e+01, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[9],                     0, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[10],                     0, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[11],                     0, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[12], 1.000000000000000e+01, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[13], 1.004970922346358e+02, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[14],                     0, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[15], 2.500000000000000e+01, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[16], 1.406959291284901e+01, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[17],                     0, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[18],                     0, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[19],                     0, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[20], 2.195986001032853e+00, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[21], 5.000000000000000e-01, 1e-7);

        wnt_cell_cycle_system3.EvaluateYDerivatives(time, initial_conditions, derivs);

        // Test derivatives are correct at t=0 for these initial conditions
        // (figures from Matlab code)
        TS_ASSERT_DELTA(derivs[0],-1.586627673253325e-02, 1e-5);
        TS_ASSERT_DELTA(derivs[1],-5.532201118824132e-05, 1e-5);
        TS_ASSERT_DELTA(derivs[2], 5.545234672425986e+00, 1e-4);
        TS_ASSERT_DELTA(derivs[3],-7.449833887043188e-03, 1e-5);
        TS_ASSERT_DELTA(derivs[4], 1.549680000000000e-02, 1e-5);
        for (unsigned i=5; i<initial_conditions.size(); i++)
        {
            TS_ASSERT_DELTA(derivs[i],0.0, 1e-9);
        }
    }

    /**
     * A test for the case mutation = 2
     * (A beta-cat delta45 mutation)
     */
    void TestVanLeeuwen2009WntSwatCellCycleEquationsBetaCatOneHit()
    {
        double time = 0.0;
        double wnt_level = 1.0;

        boost::shared_ptr<AbstractCellMutationState> p_bcat1(new BetaCateninOneHitCellMutationState);

        VanLeeuwen2009WntSwatCellCycleOdeSystem wnt_cell_cycle_system4(1, wnt_level, p_bcat1);

        std::vector<double> initial_conditions = wnt_cell_cycle_system4.GetInitialConditions();

        std::vector<double> derivs(initial_conditions.size());

        TS_ASSERT_DELTA(initial_conditions[0], 7.357000000000000e-01, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[1], 1.713000000000000e-01, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[2], 6.900000000000001e-02, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[3], 3.333333333333334e-03, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[4], 1.000000000000000e-04, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[5], 1.428571428571428e-01, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[6], 2.857142857142857e-02, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[7], 2.120643654085212e-01, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[8], 1.439678172957394e+01, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[9],                     0, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[10],                     0, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[11],                     0, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[12], 1.000000000000000e+01, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[13], 1.028341552112424e+02, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[14],                     0, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[15], 2.500000000000000e+01, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[16], 1.439678172957394e+01, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[17],                     0, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[18],                     0, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[19],                     0, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[20], 2.235636835087720e+00, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[21], 1, 1e-7);

        wnt_cell_cycle_system4.EvaluateYDerivatives(time, initial_conditions, derivs);

        // Test derivatives are correct at t=0 for these initial conditions
        // (figures from Matlab code)
        TS_ASSERT_DELTA(derivs[0],-1.586627673253325e-02, 1e-5);
        TS_ASSERT_DELTA(derivs[1],-5.532201118824132e-05, 1e-5);
        TS_ASSERT_DELTA(derivs[2],5.676110199115955e+00, 1e-4);
        TS_ASSERT_DELTA(derivs[3],-7.449833887043188e-03, 1e-5);
        TS_ASSERT_DELTA(derivs[4],1.549680000000000e-02, 1e-5);
        TS_ASSERT_DELTA(derivs[5], 0, 1e-5);
        TS_ASSERT_DELTA(derivs[6], 0, 1e-5);
        TS_ASSERT_DELTA(derivs[7], 0, 1e-5);
        TS_ASSERT_DELTA(derivs[8], -12.5, 1e-5);
        TS_ASSERT_DELTA(derivs[9], 0, 1e-5);
        TS_ASSERT_DELTA(derivs[10], 12.5, 1e-5);

        for (unsigned i=11; i<initial_conditions.size(); i++)
        {
            TS_ASSERT_DELTA(derivs[i],0.0, 1e-9);
        }
    }
    /**
     * A test for the case mutation = 3
     * (An APC -/- mutation)
     */
    void TestVanLeeuwen2009WntSwatCellCycleEquationsAPCTwoHit()
    {
        double time = 0.0;
        double wnt_level = 1.0;

        boost::shared_ptr<AbstractCellMutationState> p_apc2(new ApcTwoHitCellMutationState);

        VanLeeuwen2009WntSwatCellCycleOdeSystem wnt_cell_cycle_system5(1, wnt_level, p_apc2);

        std::vector<double> initial_conditions = wnt_cell_cycle_system5.GetInitialConditions();

        TS_ASSERT_DELTA(initial_conditions[0], 7.357000000000000e-01, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[1], 1.713000000000000e-01, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[2], 6.900000000000001e-02, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[3], 3.333333333333334e-03, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[4], 1.000000000000000e-04, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[5],                     0, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[6], 3.333333333333333e-02, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[7],                      0, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[8],                    25, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[9],                      0, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[10],                     0, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[11],                     0, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[12],                     10, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[13], 1.785714285714286e+02, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[14],                     0, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[15],                   25, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[16],                   25, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[17],                     0, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[18],                     0, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[19],                     0, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[20], 3.333333333333333e+00, 1e-7);
        TS_ASSERT_DELTA(initial_conditions[21], 1, 1e-7);

        std::vector<double> derivs(initial_conditions.size());

        wnt_cell_cycle_system5.EvaluateYDerivatives(time, initial_conditions, derivs);

        // Test derivatives are correct at t=0 for these initial conditions
        // (figures from Matlab code)
        TS_ASSERT_DELTA(derivs[0],-1.586627673253325e-02, 1e-5);
        TS_ASSERT_DELTA(derivs[1],-5.532201118824132e-05, 1e-5);
        TS_ASSERT_DELTA(derivs[2],9.917397507286381e+00, 1e-4);
        TS_ASSERT_DELTA(derivs[3],-7.449833887043188e-03, 1e-5);
        TS_ASSERT_DELTA(derivs[4],1.549680000000000e-02, 1e-5);

        for (unsigned i=5; i<initial_conditions.size(); i++)
        {
            TS_ASSERT_DELTA(derivs[i],0.0, 1e-9);
        }
    }

    void TestVanLeeuwen2009WntSwatCellCycleSolver()
    {
        double wnt_level = 1.0;
        boost::shared_ptr<AbstractCellMutationState> p_wt_state(new WildTypeCellMutationState);
        VanLeeuwen2009WntSwatCellCycleOdeSystem wnt_system(1, wnt_level, p_wt_state);

        // Solve system using rk4 solver

        // Matlab's strictest bit uses 0.01 below and relaxes it on flatter bits

        double h_value = 0.0001;

        RungeKutta4IvpOdeSolver rk4_solver;
        RungeKuttaFehlbergIvpOdeSolver rkf_solver;
        BackwardEulerIvpOdeSolver back_solver(9);

        OdeSolution solutions;

        std::vector<double> initial_conditions = wnt_system.GetInitialConditions();

        Timer::Reset();
        solutions = rk4_solver.Solve(&wnt_system, initial_conditions, 0.0, 100.0, h_value, h_value);
        Timer::Print("1. Runge-Kutta");

        // Test solutions are correct for a small time increase
        int end = solutions.rGetSolutions().size() - 1;

        // Test the simulation is ending at the right time (entering S phase at 5.971 hours)
        TS_ASSERT_DELTA(solutions.rGetTimes()[end], 6.193774457302713, 1e-2);

        // Proper values calculated using the Matlab stiff ODE solver ode15s. Note that
        // large tolerances are required for the tests to pass with both chaste and CVODE solvers
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][0], 2.937457584307182e-01, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][1], 1.000117721173146e+00, 1.01e-2);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][2], 2.400566063511684e+00, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][3], 1.392886957364369e+00, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][4], 1.358673176849681e-01, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][5], 1.428571428571429e-01, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][6], 2.857142857142857e-02, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][7], 2.120643654085166e-01, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][8], 1.439678172957273e+01, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][9], 0.0, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][10], 0.0, 1.01e-2);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][11], 0.0, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][12], 1.000000000000039e+01, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][13], 1.028341552112378e+02, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][14], 0.0, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][15], 2.500000000000000e+01, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][16], 1.439678172957305e+01, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][17], 0.0, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][18], 0.0, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][19], 0.0, 1.01e-2);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][20], 2.235636835087647e+00, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][21], 1.0, 1e-3);
    }

    void TestVanLeeuwen2009WntSwatCellCycleSolverWithAPCSingleHit()
    {
        double wnt_level = 1.0;
        boost::shared_ptr<AbstractCellMutationState> p_apc1(new ApcOneHitCellMutationState);
        VanLeeuwen2009WntSwatCellCycleOdeSystem wnt_system(1, wnt_level, p_apc1);

        // Solve system using rk4 solver

        // Matlab's strictest bit uses 0.01 below and relaxes it on flatter bits.

        double h_value = 0.0001;

        RungeKutta4IvpOdeSolver rk4_solver;
        OdeSolution solutions;

        std::vector<double> initial_conditions = wnt_system.GetInitialConditions();

        solutions = rk4_solver.Solve(&wnt_system, initial_conditions, 0.0, 100.0, h_value, h_value);

        // Test solutions are correct for a small time increase
        int end = solutions.rGetSolutions().size() - 1;

        // Test the simulation is ending at the right time (entering S phase at 3.94 hours)
        TS_ASSERT_DELTA(solutions.rGetTimes()[end], 4.722377242770206, 1e-2);

        // Proper values calculated using the Matlab stiff ODE solver ode15s. Note that
        // large tolerances are required for the tests to pass with both chaste and CVODE solvers
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][0], 2.449671985497571e-01, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][1], 1.00, 1.01e-2);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][2], 2.970519972449688, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][3], 1.962479928865283, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][4], 1.594820065531480e-01, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][5], 7.692307692307693e-02, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][6], 3.076923076923077e-02, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][7], 1.216821534917367e-01, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][8], 1.891589232541249e+01, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][9], 0.0, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][10], 0.0, 1.01e-2);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][11], 0.0, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][12], 1.000000000000039e+01, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][13], 1.351135166100946e+02, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][14], 0.0, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][15], 2.500000000000000e+01, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][16], 1.891589232541281e+01, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][17], 0.0, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][18], 0.0, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][19], 0.0, 1.01e-2);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][20], 2.744779424184835e+00, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][21], 1.0, 1e-3);
    }

    void TestVanLeeuwen2009WntSwatCellCycleSolverWithBetaCateninHit()
    {
        double wnt_level = 1.0;
        boost::shared_ptr<AbstractCellMutationState> p_bcat1(new BetaCateninOneHitCellMutationState);
        VanLeeuwen2009WntSwatCellCycleOdeSystem wnt_system(1, wnt_level, p_bcat1);

        // Solve system using rk4 solver

        // Matlab's strictest bit uses 0.01 below and relaxes it on flatter bits.

        double h_value = 0.0001;

        RungeKutta4IvpOdeSolver rk4_solver;
        OdeSolution solutions;

        std::vector<double> initial_conditions = wnt_system.GetInitialConditions();

        solutions = rk4_solver.Solve(&wnt_system, initial_conditions, 0.0, 100.0, h_value, h_value);

        // Test solutions are correct for a small time increase
        int end = solutions.rGetSolutions().size() - 1;

        // Test the simulation is ending at the right time (entering S phase at 7.81 hours)
        TS_ASSERT_DELTA(solutions.rGetTimes()[end], 6.109381124487460, 1e-2);

        // Proper values calculated using the Matlab stiff ODE solver ode15s. Note that
        // large tolerances are required for the tests to pass with both chaste and CVODE solvers
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][0], 2.885925504994788e-01, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][1], 1.0, 1.01e-2);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][2], 2.461454077577112e+00, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][3], 1.443615084885885e+00, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][4], 1.383019697069664e-01, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][5], 1.428571428571428e-01, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][6], 2.857142857142857e-02, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][7], 1.826340175985522e-01, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][8], 8.847117964569135e+00, 1e-2);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][9], 0, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][10], 6.199836467090596e+00, 1.01e-2);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][11], 0, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][12], 9.735946067044630e+00, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][13], 6.153733439469387e+01, 1e-1);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][14], 4.310129026467735e+01, 1e-1);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][15], 2.476909200705334e+01, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][16], 8.766188466378457e+00, 1e-2);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][17], 0, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][18], 6.141626362262133e+00, 1e-2);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][19], 0, 1.01e-2);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][20], 2.283368993736725e+00, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][21], 1, 1e-3);
    }

    void TestVanLeeuwen2009WntSwatCellCycleSolverWithAPCDoubleHit()
    {
        double wnt_level = 1.0;
        boost::shared_ptr<AbstractCellMutationState> p_apc2(new ApcTwoHitCellMutationState);
        VanLeeuwen2009WntSwatCellCycleOdeSystem wnt_system(1, wnt_level, p_apc2);

        // Solve system using rk4 solver

        // Matlab's strictest bit uses 0.01 below and relaxes it on flatter bits.

        double h_value = 0.0001;

        RungeKutta4IvpOdeSolver rk4_solver;
        OdeSolution solutions;

        std::vector<double> initial_conditions = wnt_system.GetInitialConditions();

        solutions = rk4_solver.Solve(&wnt_system, initial_conditions, 0.0, 100.0, h_value, h_value);

        // Test solutions are correct for a small time increase
        int end = solutions.rGetSolutions().size() - 1;

        // Test the simulation is ending at the right time (entering S phase at 3.94 hours)
        TS_ASSERT_DELTA(solutions.rGetTimes()[end], 3.912928619944499e+00, 1e-2);

        // Proper values calculated using the Matlab stiff ODE solver ode15s. Note that
        // large tolerances are required for the tests to pass with both chaste and CVODE solvers
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][0], 2.040531616988712e-01, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][1], 1.000000000000046e+00, 1.01e-2);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][2], 3.738096511672503e+00, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][3], 2.726370715594771e+00, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][4], 1.844725907694925e-01, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][5], 0, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][6],  3.333333333333333e-02, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][7], 0, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][8], 2.499999999999840e+01, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][9], 0, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][10], 0, 1.01e-2);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][11], 0, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][12], 1.000000000000064e+01, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][13], 1.785714285714286e+02, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][14], 0, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][15], 2.500000000000076e+01, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][16], 2.499999999999917e+01, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][17], 0, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][18], 0, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][19], 0, 1.01e-2);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][20], 3.333333333333310e+00, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][21], 1, 1e-3);
    }

    /**
     * Test to check that none of the protein concentrations in the
     * VanLeeuwen2009WntSwatCellCycleModel go negative in a crypt simulation
     * (see ticket #629).
     */
    void TestVanLeeuwen2009WntOdeSolutionDoesNotGoNegative()
    {
        double time_of_each_run = 0.01; // for each run

        unsigned cells_across = 4;
        unsigned cells_up = 25;
        double crypt_width = 4.1;
        double crypt_length = 22.0;
        unsigned thickness_of_ghost_layer = 3;

        CylindricalHoneycombMeshGenerator generator(cells_across, cells_up,thickness_of_ghost_layer, crypt_width/cells_across);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);

        // Set up cells
        std::vector<CellPtr> cells;
        CryptCellsGenerator<VanLeeuwen2009WntSwatCellCycleModelHypothesisOne> cells_generator;
        cells_generator.Generate(cells, p_mesh, location_indices, true);

        for (unsigned i=0; i<cells.size(); i++)
        {
            cells[i]->SetBirthTime(-1.1); // just to make the test run a bit quicker
        }

        MeshBasedCellPopulationWithGhostNodes<2> crypt(*p_mesh, cells, location_indices);

        // Set cell population to output cell types
        crypt.AddCellPopulationCountWriter<CellMutationStatesCountWriter>();

        WntConcentration<2>::Instance()->SetType(LINEAR);
        WntConcentration<2>::Instance()->SetWntConcentrationParameter(1.0/3.0);
        WntConcentration<2>::Instance()->SetCellPopulation(crypt);
        WntConcentration<2>::Instance()->SetCryptLength(crypt_length);

        CryptSimulation2d simulator(crypt);
        simulator.SetOutputDirectory("IngeCellsNiceCryptSim_long");

        // Set length of simulation here
        simulator.SetEndTime(time_of_each_run);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        simulator.AddForce(p_linear_force);

        // Add cell Killer
        MAKE_PTR_ARGS(SloughingCellKiller<2>, p_killer, (&simulator.rGetCellPopulation(), crypt_length));
        simulator.AddCellKiller(p_killer);

        TS_ASSERT_THROWS_NOTHING(simulator.Solve());

        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
        WntConcentration<2>::Destroy();
    }

    void TestArchiving()
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "vanleeuwen_ode.arch";

        {
            double wnt_level = 1.0;
            boost::shared_ptr<AbstractCellMutationState> p_bcat1(new BetaCateninOneHitCellMutationState);
            VanLeeuwen2009WntSwatCellCycleOdeSystem ode_system(1, wnt_level, p_bcat1);

            TS_ASSERT_DELTA(ode_system.GetWntLevel(), 1.0, 1e-6);
            TS_ASSERT_EQUALS(ode_system.GetMutationState()->IsType<BetaCateninOneHitCellMutationState>(), true);
            TS_ASSERT_EQUALS(ode_system.GetHypothesis(), 1u);

            std::vector<double> initial_conditions = ode_system.GetInitialConditions();
            TS_ASSERT_EQUALS(initial_conditions.size(), 22u);
            TS_ASSERT_DELTA(initial_conditions[0],    0.73570000, 1e-7);
            TS_ASSERT_DELTA(initial_conditions[1],    0.17130000, 1e-7);
            TS_ASSERT_DELTA(initial_conditions[2],    0.06900000, 1e-7);
            TS_ASSERT_DELTA(initial_conditions[3],    0.00333333, 1e-7);
            TS_ASSERT_DELTA(initial_conditions[4],    0.00010000, 1e-7);
            TS_ASSERT_DELTA(initial_conditions[5],    0.14285714, 1e-7);
            TS_ASSERT_DELTA(initial_conditions[6],    0.02857142, 1e-7);
            TS_ASSERT_DELTA(initial_conditions[7],    0.21206436, 1e-7);
            TS_ASSERT_DELTA(initial_conditions[8],   14.39678172, 1e-7);
            TS_ASSERT_DELTA(initial_conditions[9],    0.00000000, 1e-7);
            TS_ASSERT_DELTA(initial_conditions[10],   0.00000000, 1e-7);
            TS_ASSERT_DELTA(initial_conditions[11],   0.00000000, 1e-7);
            TS_ASSERT_DELTA(initial_conditions[12],  10.00000000, 1e-7);
            TS_ASSERT_DELTA(initial_conditions[13], 102.83415521, 1e-7);
            TS_ASSERT_DELTA(initial_conditions[14],   0.00000000, 1e-7);
            TS_ASSERT_DELTA(initial_conditions[15],  25.00000000, 1e-7);
            TS_ASSERT_DELTA(initial_conditions[16],  14.39678172, 1e-7);
            TS_ASSERT_DELTA(initial_conditions[17],   0.00000000, 1e-7);
            TS_ASSERT_DELTA(initial_conditions[18],   0.00000000, 1e-7);
            TS_ASSERT_DELTA(initial_conditions[19],   0.00000000, 1e-7);
            TS_ASSERT_DELTA(initial_conditions[20],   2.23563683, 1e-7);
            TS_ASSERT_DELTA(initial_conditions[21],   1.00000000, 1e-7);

            // Create an output archive
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Archive ODE system
            AbstractOdeSystem* const p_const_ode_system = &ode_system;
            output_arch << p_const_ode_system;
        }

        {
            AbstractOdeSystem* p_ode_system;

            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            // Restore from the archive
            input_arch >> p_ode_system;

            // Check that archiving worked correctly
            TS_ASSERT_DELTA(static_cast<VanLeeuwen2009WntSwatCellCycleOdeSystem*>(p_ode_system)->GetWntLevel(), 1.0, 1e-6);
            TS_ASSERT_EQUALS(static_cast<VanLeeuwen2009WntSwatCellCycleOdeSystem*>(p_ode_system)->GetMutationState()->IsType<BetaCateninOneHitCellMutationState>(), true);
            TS_ASSERT_EQUALS(static_cast<VanLeeuwen2009WntSwatCellCycleOdeSystem*>(p_ode_system)->GetHypothesis(), true);

            std::vector<double> initial_conditions = p_ode_system->GetInitialConditions();
            TS_ASSERT_EQUALS(initial_conditions.size(), 22u);
            TS_ASSERT_DELTA(initial_conditions[0],    0.73570000, 1e-7);
            TS_ASSERT_DELTA(initial_conditions[1],    0.17130000, 1e-7);
            TS_ASSERT_DELTA(initial_conditions[2],    0.06900000, 1e-7);
            TS_ASSERT_DELTA(initial_conditions[3],    0.00333333, 1e-7);
            TS_ASSERT_DELTA(initial_conditions[4],    0.00010000, 1e-7);
            TS_ASSERT_DELTA(initial_conditions[5],    0.14285714, 1e-7);
            TS_ASSERT_DELTA(initial_conditions[6],    0.02857142, 1e-7);
            TS_ASSERT_DELTA(initial_conditions[7],    0.21206436, 1e-7);
            TS_ASSERT_DELTA(initial_conditions[8],   14.39678172, 1e-7);
            TS_ASSERT_DELTA(initial_conditions[9],    0.00000000, 1e-7);
            TS_ASSERT_DELTA(initial_conditions[10],   0.00000000, 1e-7);
            TS_ASSERT_DELTA(initial_conditions[11],   0.00000000, 1e-7);
            TS_ASSERT_DELTA(initial_conditions[12],  10.00000000, 1e-7);
            TS_ASSERT_DELTA(initial_conditions[13], 102.83415521, 1e-7);
            TS_ASSERT_DELTA(initial_conditions[14],   0.00000000, 1e-7);
            TS_ASSERT_DELTA(initial_conditions[15],  25.00000000, 1e-7);
            TS_ASSERT_DELTA(initial_conditions[16],  14.39678172, 1e-7);
            TS_ASSERT_DELTA(initial_conditions[17],   0.00000000, 1e-7);
            TS_ASSERT_DELTA(initial_conditions[18],   0.00000000, 1e-7);
            TS_ASSERT_DELTA(initial_conditions[19],   0.00000000, 1e-7);
            TS_ASSERT_DELTA(initial_conditions[20],   2.23563683, 1e-7);
            TS_ASSERT_DELTA(initial_conditions[21],   1.00000000, 1e-7);

            // Tidy up
            delete p_ode_system;
        }
    }
};

#endif /*TESTVANLEEUWEN2009WNTSWATCELLCYCLEODESYSTEM_HPP_*/
