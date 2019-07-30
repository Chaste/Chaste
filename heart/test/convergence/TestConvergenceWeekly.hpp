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


#ifndef TESTCONVERGENCEWEEKLY_HPP_
#define TESTCONVERGENCEWEEKLY_HPP_

#include <cxxtest/TestSuite.h>
#include "BidomainProblem.hpp"
#include "MonodomainProblem.hpp"
#include <petscvec.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>

#include "LuoRudy1991BackwardEuler.hpp"
#include "LuoRudy1991.hpp"
#include "PdeConvergenceTester.hpp"
#include "SpaceConvergenceTester.hpp"
#include "OdeConvergenceTester.hpp"
#include "PetscSetupAndFinalize.hpp"


class TestConvergenceWeekly : public CxxTest::TestSuite
{
public:

    void xxTest3DSpace()
    {

        SpaceConvergenceTester<CellLuoRudy1991FromCellMLBackwardEuler, BidomainProblem<3>, 3, 2> tester;
        HeartConfig::Instance()->SetUseRelativeTolerance(1e-8);
        //tester.SetKspRelativeTolerance(1e-8);
        tester.SetMeshWidth(0.15);//cm
        tester.Converge(__FUNCTION__);
        TS_ASSERT(tester.Converged);
        TS_ASSERT_EQUALS(tester.MeshNum, 4u); ///Just to prove the thing works
        HeartConfig::Instance()->Reset();
    }

    //Experiments with ksp_atol follow.
    //This first one has to be done with GMRES as 1D are known to be a bit flakey
    void TestSpaceConvergencein1DWithAtol()
    {
        HeartConfig::Instance()->SetKSPSolver("gmres");
        HeartConfig::Instance()->SetKSPPreconditioner("jacobi");
        SpaceConvergenceTester<CellLuoRudy1991FromCellMLBackwardEuler, BidomainProblem<1>, 1, 2> tester;
        HeartConfig::Instance()->SetUseAbsoluteTolerance(1e-5);
       //tester.SetKspAbsoluteTolerance(1e-5);
        tester.Converge(__FUNCTION__);
        TS_ASSERT(tester.Converged);
        TS_ASSERT_EQUALS(tester.MeshNum, 5u);
        TS_ASSERT_LESS_THAN(tester.LastDifference, 0.0041039);
        //Has to be at least as good as the 1D with Rtol=1e-7
        //Note the final line fails with ksp_atol=1e-4
        HeartConfig::Instance()->Reset();
    }


    void Test3DSpace10()
    {
        HeartConfig::Instance()->SetKSPSolver("symmlq");
        HeartConfig::Instance()->SetKSPPreconditioner("bjacobi");
        PetscTools::SetOption("-log_summary", "");
        SpaceConvergenceTester<CellLuoRudy1991FromCellMLBackwardEuler, BidomainProblem<3>, 3, 2> tester;
        HeartConfig::Instance()->SetUseAbsoluteTolerance(1e-7);
        //tester.SetKspAbsoluteTolerance(1e-7);
        tester.OdeTimeStep /= 2.0;
        tester.PdeTimeStep /= 2.0;
        tester.SetMeshWidth(0.10);//cm

        tester.Converge(__FUNCTION__);
        TS_ASSERT(tester.Converged);
        TS_ASSERT_EQUALS(tester.MeshNum, 3u);
        HeartConfig::Instance()->Reset();
    }

    void Test3DSpace10RampedQuarterStimulus()
    {
        HeartConfig::Instance()->SetKSPSolver("symmlq");
        HeartConfig::Instance()->SetKSPPreconditioner("bjacobi");
        PetscTools::SetOption("-log_summary", "");
        SpaceConvergenceTester<CellLuoRudy1991FromCellMLBackwardEuler, BidomainProblem<3>, 3, 2> tester;
        HeartConfig::Instance()->SetUseAbsoluteTolerance(1e-7);
        //tester.SetKspAbsoluteTolerance(1e-7);
        tester.OdeTimeStep /= 2.0;
        tester.PdeTimeStep /= 2.0;
        tester.SetMeshWidth(0.10);//cm
        tester.Stimulus = QUARTER;

        tester.Converge(__FUNCTION__);
        TS_ASSERT(tester.Converged);
        TS_ASSERT_EQUALS(tester.MeshNum, 3u);
        HeartConfig::Instance()->Reset();
    }

    //More experiments with ksp_atol follow.
    void TestSpaceConvergencein2DWithAtol()
    {
        HeartConfig::Instance()->SetKSPSolver("symmlq");
        HeartConfig::Instance()->SetKSPPreconditioner("bjacobi");
        SpaceConvergenceTester<CellLuoRudy1991FromCellMLBackwardEuler, BidomainProblem<2>, 2, 2> tester;
        //tester.SetKspAbsoluteTolerance(1e-5);
        HeartConfig::Instance()->SetUseAbsoluteTolerance(1e-5);
        tester.Converge(__FUNCTION__);
        TS_ASSERT(tester.Converged);
        TS_ASSERT_EQUALS(tester.MeshNum, 5u);
        TS_ASSERT_LESS_THAN(tester.LastDifference, 0.0081583);
        //Comes in at 1.17118e-5
        //Has to be at least as good as the 2D with Rtol=5e-8
        HeartConfig::Instance()->Reset();

    }

    //Copied from projects/jmpf since this converges on mesh4
    void Test3DSpaceRelaxWidthWithAtol()
    {
        HeartConfig::Instance()->SetKSPSolver("symmlq");
        HeartConfig::Instance()->SetKSPPreconditioner("bjacobi");
        SpaceConvergenceTester<CellLuoRudy1991FromCellMLBackwardEuler, BidomainProblem<3>, 3, 2> tester;
        //tester.SetKspAbsoluteTolerance(1e-3);
        HeartConfig::Instance()->SetUseAbsoluteTolerance(1e-3);
        tester.SetMeshWidth(0.15);//cm
        tester.Converge(__FUNCTION__);
        TS_ASSERT(tester.Converged);
        TS_ASSERT_EQUALS(tester.MeshNum, 4u);
        HeartConfig::Instance()->Reset();
    }
};

#endif /*TESTCONVERGENCEWEEKLY_HPP_*/
