/*

Copyright (C) University of Oxford, 2005-2011

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

    void xxTest3DSpace() throw(Exception)
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
    void TestSpaceConvergencein1DWithAtol() throw(Exception)
    {
        HeartConfig::Instance()->SetKSPSolver("gmres");
        HeartConfig::Instance()->SetKSPPreconditioner("jacobi");
        SpaceConvergenceTester<CellLuoRudy1991FromCellMLBackwardEuler, BidomainProblem<1>, 1, 2> tester;
        HeartConfig::Instance()->SetUseAbsoluteTolerance(1e-5);
       //tester.SetKspAbsoluteTolerance(1e-5);
        tester.Converge(__FUNCTION__);
        TS_ASSERT(tester.Converged);
        TS_ASSERT_EQUALS(tester.MeshNum, 5u);
        TS_ASSERT_LESS_THAN(tester.LastDifference, 1.68417e-05);
        //Has to be at least as good as the 1D with Rtol=1e-7
        //Note the final line fails with ksp_atol=1e-4
        HeartConfig::Instance()->Reset();
    }


    void Test3DSpace10() throw(Exception)
    {
        HeartConfig::Instance()->SetKSPSolver("symmlq");
        HeartConfig::Instance()->SetKSPPreconditioner("bjacobi");
        PetscOptionsSetValue("-log_summary", "");
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

    void Test3DSpace10RampedQuarterStimulus() throw(Exception)
    {
        HeartConfig::Instance()->SetKSPSolver("symmlq");
        HeartConfig::Instance()->SetKSPPreconditioner("bjacobi");
        PetscOptionsSetValue("-log_summary", "");
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
    void TestSpaceConvergencein2DWithAtol() throw(Exception)
    {
        HeartConfig::Instance()->SetKSPSolver("symmlq");
        HeartConfig::Instance()->SetKSPPreconditioner("bjacobi");
        SpaceConvergenceTester<CellLuoRudy1991FromCellMLBackwardEuler, BidomainProblem<2>, 2, 2> tester;
        //tester.SetKspAbsoluteTolerance(1e-5);
        HeartConfig::Instance()->SetUseAbsoluteTolerance(1e-5);
        tester.Converge(__FUNCTION__);
        TS_ASSERT(tester.Converged);
        TS_ASSERT_EQUALS(tester.MeshNum, 5u);
        TS_ASSERT_LESS_THAN(tester.LastDifference, 6.65582e-05);
        //Comes in at 1.17118e-5
        //Has to be at least as good as the 2D with Rtol=5e-8
        HeartConfig::Instance()->Reset();

    }

    //Copied from projects/jmpf since this converges on mesh4
    void Test3DSpaceRelaxWidthWithAtol() throw(Exception)
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
