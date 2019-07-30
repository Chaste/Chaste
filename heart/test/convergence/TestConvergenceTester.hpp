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


#ifndef TESTCONVERGENCETESTER_HPP_
#define TESTCONVERGENCETESTER_HPP_

#include <cxxtest/TestSuite.h>
#include "BidomainProblem.hpp"
#include "MonodomainProblem.hpp"
#include <petscvec.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>

#include "LuoRudy1991BackwardEuler.hpp"
#include "PdeConvergenceTester.hpp"
#include "SpaceConvergenceTester.hpp"
#include "KspConvergenceTester.hpp"
#include "OdeConvergenceTester.hpp"
#include "OdePdeConvergenceTester.hpp"
#include "PetscSetupAndFinalize.hpp"


class TestConvergenceTester : public CxxTest::TestSuite
{
public:
//    time_t start;
//    void setUp()
//    {
//        start=time(NULL);
//    }
//    void tearDown()
//    {
//        unsigned secs = time(NULL) - start;
//        std::cout<<"REAL TIME Test took "<< secs/60<<" minutes, "<< secs%60<<" seconds of real time\n";
//    }
    void Test1DOdeTime()
    {
        OdeConvergenceTester<CellLuoRudy1991FromCellMLBackwardEuler, MonodomainProblem<1>, 1, 1> tester;
        tester.MeshNum=1;
        tester.AbsoluteStimulus = -5e6; // The default of -1e7 causes V to go out of range for lookup tables
        tester.Converge(__FUNCTION__);
        TS_ASSERT(tester.Converged);
        TS_ASSERT(tester.IsConverged());
        TS_ASSERT_EQUALS(tester.OdeTimeStep, 0.0025);
    }

    void Test1DOdeTimeWarning()
    {
        OdeConvergenceTester<CellLuoRudy1991FromCellMLBackwardEuler, MonodomainProblem<1>, 1, 1> tester;
        tester.MeshNum = 1;
        tester.AbsoluteStimulus = -5e8; // We want V to go out of range for lookup tables
        tester.Converge(__FUNCTION__);

        TS_ASSERT_EQUALS(Warnings::Instance()->GetNumWarnings(), 10u);
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNextWarningMessage(),
                         "This run threw an exception.  Check convergence results\n");
        Warnings::Instance()->QuietDestroy();
    }

    void Test1DPdeTime()
    {
        PdeConvergenceTester<CellLuoRudy1991FromCellMLBackwardEuler, MonodomainProblem<1>, 1, 1> tester;
        tester.MeshNum=1;
        tester.AbsoluteStimulus = -8e6; // The default of -1e7 causes V to go out of range for lookup tables
        tester.RelativeConvergenceCriterion=3e-2;
        tester.Converge(__FUNCTION__);
        TS_ASSERT(tester.Converged);
        TS_ASSERT_EQUALS(tester.PdeTimeStep, 0.02); // Was 0.01 before using lookup tables
    }

    void Test1DOdePdeTime()
    {
        OdePdeConvergenceTester<CellLuoRudy1991FromCellMLBackwardEuler, MonodomainProblem<1>, 1, 1> tester;
        tester.MeshNum=1;
        tester.AbsoluteStimulus = -5e6; // The default of -1e7 causes V to go out of range for lookup tables
        tester.RelativeConvergenceCriterion=0.026458;
        tester.Converge(__FUNCTION__);
        TS_ASSERT(tester.Converged);
        TS_ASSERT_EQUALS(tester.PdeTimeStep, 0.005);
    }

    void Test1DPdeTimeRegion()
    {
        PdeConvergenceTester<CellLuoRudy1991FromCellMLBackwardEuler, MonodomainProblem<1>, 1, 1> tester;
        tester.MeshNum=1;
        tester.AbsoluteStimulus = -5e6; // The default of -1e7 causes V to go out of range for lookup tables
        tester.Stimulus=QUARTER;
        tester.RelativeConvergenceCriterion=0.022361;
        tester.Converge(__FUNCTION__);
        TS_ASSERT(tester.Converged);
        TS_ASSERT_EQUALS(tester.PdeTimeStep, 0.01);
    }

    void Test1DPdeTimeNeumann()
    {
        PdeConvergenceTester<CellLuoRudy1991FromCellMLBackwardEuler, MonodomainProblem<1>, 1, 1> tester;
        tester.MeshNum=1;
        tester.Stimulus=NEUMANN;
        tester.RelativeConvergenceCriterion=0.022361;
        tester.Converge(__FUNCTION__);
        TS_ASSERT(tester.Converged);
        TS_ASSERT_EQUALS(tester.PdeTimeStep, 0.02);
    }

    void TestSpaceConvergenceMonoIn1DWithRelativeTolerance()
    {
        SpaceConvergenceTester<CellLuoRudy1991FromCellMLBackwardEuler, MonodomainProblem<1>, 1, 1> tester;
        HeartConfig::Instance()->SetUseRelativeTolerance(1e-4);
        tester.RelativeConvergenceCriterion=0.14142;
        tester.AbsoluteStimulus = -5e6; // The default of -1e7 causes V to go out of range for lookup tables
        tester.Converge(__FUNCTION__);
        TS_ASSERT(tester.IsConverged());
        TS_ASSERT_EQUALS(tester.GetMeshNum(), 2);
        TS_ASSERT_DELTA(tester.GetSpaceStep(), 0.0125, 1e-8);
        TS_ASSERT_LESS_THAN(tester.LastDifference, 0.070711);
        HeartConfig::Instance()->Reset();
    }

    void TestSpaceConvergenceBidomainIn1DWithAbsoluteTolerance()
    {
        // Zero pivot detected in Cholesky factorisation for mesh 1. This is not an error and it may always happen when using bjacobi with singular systems.
        HeartConfig::Instance()->SetKSPPreconditioner("jacobi");

        SpaceConvergenceTester<CellLuoRudy1991FromCellMLBackwardEuler, BidomainProblem<1>, 1, 2> tester;
        HeartConfig::Instance()->SetUseAbsoluteTolerance(1e-4);
        tester.RelativeConvergenceCriterion=0.14142;
        tester.AbsoluteStimulus = -5.5e6; // The default of -1e7 causes V to go out of range for lookup tables
        tester.Converge(__FUNCTION__);
        TS_ASSERT(tester.IsConverged());
        TS_ASSERT_EQUALS(tester.GetMeshNum(), 2);
        TS_ASSERT_DELTA(tester.GetSpaceStep(), 0.0125, 1e-8);
        TS_ASSERT_LESS_THAN(tester.LastDifference, 0.062450);
        HeartConfig::Instance()->Reset();
     }
};

#endif /*TESTCONVERGENCETESTER_HPP_*/
