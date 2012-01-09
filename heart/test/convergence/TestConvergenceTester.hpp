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
    void Test1DOdeTime() throw(Exception)
    {
        OdeConvergenceTester<CellLuoRudy1991FromCellMLBackwardEuler, MonodomainProblem<1>, 1, 1> tester;
        tester.MeshNum=1;
        tester.AbsoluteStimulus = -5e6; // The default of -1e7 causes V to go out of range for lookup tables
        tester.Converge(__FUNCTION__);
        TS_ASSERT(tester.Converged);
        TS_ASSERT(tester.IsConverged());
        TS_ASSERT_EQUALS(tester.OdeTimeStep, 0.0025);
    }

    void Test1DOdeTimeWarning() throw(Exception)
    {
        OdeConvergenceTester<CellLuoRudy1991FromCellMLBackwardEuler, MonodomainProblem<1>, 1, 1> tester;
        tester.MeshNum = 1;
        tester.AbsoluteStimulus = -5e8; // We want V to go out of range for lookup tables
        tester.Converge(__FUNCTION__);

        TS_ASSERT_EQUALS(Warnings::Instance()->GetNumWarnings(), 2u);
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNextWarningMessage(),
                         "This run threw an exception.  Check convergence results\n");
    }

    void Test1DPdeTime() throw(Exception)
    {
        PdeConvergenceTester<CellLuoRudy1991FromCellMLBackwardEuler, MonodomainProblem<1>, 1, 1> tester;
        tester.MeshNum=1;
        tester.AbsoluteStimulus = -8e6; // The default of -1e7 causes V to go out of range for lookup tables
        tester.RelativeConvergenceCriterion=9e-4;
        tester.Converge(__FUNCTION__);
        TS_ASSERT(tester.Converged);
        TS_ASSERT_EQUALS(tester.PdeTimeStep, 0.02); // Was 0.01 before using lookup tables
    }

    void Test1DOdePdeTime() throw(Exception)
    {
        OdePdeConvergenceTester<CellLuoRudy1991FromCellMLBackwardEuler, MonodomainProblem<1>, 1, 1> tester;
        tester.MeshNum=1;
        tester.AbsoluteStimulus = -5e6; // The default of -1e7 causes V to go out of range for lookup tables
        tester.RelativeConvergenceCriterion=7e-4;
        tester.Converge(__FUNCTION__);
        TS_ASSERT(tester.Converged);
        TS_ASSERT_EQUALS(tester.PdeTimeStep, 0.005);
    }

    void Test1DPdeTimeRegion() throw(Exception)
    {
        PdeConvergenceTester<CellLuoRudy1991FromCellMLBackwardEuler, MonodomainProblem<1>, 1, 1> tester;
        tester.MeshNum=1;
        tester.AbsoluteStimulus = -5e6; // The default of -1e7 causes V to go out of range for lookup tables
        tester.Stimulus=QUARTER;
        tester.RelativeConvergenceCriterion=5e-4;
        tester.Converge(__FUNCTION__);
        TS_ASSERT(tester.Converged);
        TS_ASSERT_EQUALS(tester.PdeTimeStep, 0.01);
    }

    void Test1DPdeTimeNeumann() throw(Exception)
    {
        PdeConvergenceTester<CellLuoRudy1991FromCellMLBackwardEuler, MonodomainProblem<1>, 1, 1> tester;
        tester.MeshNum=1;
        tester.Stimulus=NEUMANN;
        tester.RelativeConvergenceCriterion=5e-4;
        tester.Converge(__FUNCTION__);
        TS_ASSERT(tester.Converged);
        TS_ASSERT_EQUALS(tester.PdeTimeStep, 0.02);
    }



    void TestSpaceConvergenceMonoIn1DWithRelativeTolerance() throw(Exception)
    {
        SpaceConvergenceTester<CellLuoRudy1991FromCellMLBackwardEuler, MonodomainProblem<1>, 1, 1> tester;
        HeartConfig::Instance()->SetUseRelativeTolerance(1e-4);
        tester.RelativeConvergenceCriterion=2e-2;
        tester.AbsoluteStimulus = -5e6; // The default of -1e7 causes V to go out of range for lookup tables
        tester.Converge(__FUNCTION__);
        TS_ASSERT(tester.IsConverged());
        TS_ASSERT_EQUALS(tester.GetMeshNum(), 2);
        TS_ASSERT_DELTA(tester.GetSpaceStep(), 0.0125, 1e-8);
        TS_ASSERT_LESS_THAN(tester.LastDifference, 5e-3);
        HeartConfig::Instance()->Reset();
    }

    void xxTestSpaceConvergenceBidomainIn1DWithAbsoluteTolerance() throw(Exception)
    {
        // Zero pivot detected in Cholesky factorisation for mesh 1. This is not an error and it may always happen when using bjacobi with singular systems.
        HeartConfig::Instance()->SetKSPPreconditioner("jacobi");

        SpaceConvergenceTester<CellLuoRudy1991FromCellMLBackwardEuler, BidomainProblem<1>, 1, 2> tester;
        HeartConfig::Instance()->SetUseAbsoluteTolerance(1e-4);
        tester.RelativeConvergenceCriterion=2e-2;
        tester.AbsoluteStimulus = -5.5e6; // The default of -1e7 causes V to go out of range for lookup tables
        tester.Converge(__FUNCTION__);
        TS_ASSERT(tester.IsConverged());
        TS_ASSERT_EQUALS(tester.GetMeshNum(), 2);
        TS_ASSERT_DELTA(tester.GetSpaceStep(), 0.0125, 1e-8);
        TS_ASSERT_LESS_THAN(tester.LastDifference, 3.9e-3); // Was 2.9e-3 before using lookup tables
        HeartConfig::Instance()->Reset();
     }
};

#endif /*TESTCONVERGENCETESTER_HPP_*/
