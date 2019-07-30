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


#ifndef TESTCONVERGENCENIGHTLY_HPP_
#define TESTCONVERGENCENIGHTLY_HPP_

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
#include "KspConvergenceTester.hpp"
#include "OdeConvergenceTester.hpp"
#include "OdePdeConvergenceTester.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "NobleVargheseKohlNoble1998aBackwardEuler.hpp"
#include "NobleVargheseKohlNoble1998a.hpp"
#include "NobleVargheseKohlNoble1998aOpt.hpp"


class TestConvergenceNightly : public CxxTest::TestSuite
{

public:

    void RunConvergenceTester(AbstractUntemplatedConvergenceTester* pTester, StimulusType stimulusType)
    {
        pTester->Stimulus = stimulusType;
        HeartConfig::Instance()->SetUseAbsoluteTolerance(5e-4);

        pTester->Converge("Automated_test");
        TS_ASSERT(pTester->Converged);
    }

    void ConvergeInVarious(StimulusType stimulusType)
    {
        TS_ASSERT(strcmp(HeartConfig::Instance()->GetKSPPreconditioner(), "bjacobi")==0);
        TS_ASSERT(strcmp(HeartConfig::Instance()->GetKSPSolver(), "cg")==0);
        HeartConfig::Instance()->SetKSPPreconditioner("jacobi");
        HeartConfig::Instance()->SetKSPSolver("gmres");
        {
            std::cout << "PdeConvergenceTester<CellLuoRudy1991FromCellMLBackwardEuler, BidomainProblem<1>, 1, 2>\n";
            PdeConvergenceTester<CellLuoRudy1991FromCellMLBackwardEuler, BidomainProblem<1>, 1, 2> tester;
            RunConvergenceTester(&tester, stimulusType);
            TS_ASSERT_DELTA(tester.PdeTimeStep, 5.0e-3, 1e-10);
        }

        {
            std::cout << "SpaceConvergenceTester<CellLuoRudy1991FromCellMLBackwardEuler, BidomainProblem<1>, 1, 2>\n";
            //Block Jacobi with CG can detect zero pivots in a 1-D convergence test
            SpaceConvergenceTester<CellLuoRudy1991FromCellMLBackwardEuler, BidomainProblem<1>, 1, 2> tester;
            RunConvergenceTester(&tester, stimulusType);
            TS_ASSERT_EQUALS(tester.MeshNum, 5u);
        }

        {
            std::cout << "KspConvergenceTester<CellLuoRudy1991FromCellMLBackwardEuler, BidomainProblem<1>, 1, 2>\n";
            TS_ASSERT_EQUALS(HeartConfig::Instance()->GetAbsoluteTolerance(), 5.0e-4);
            KspConvergenceTester<CellLuoRudy1991FromCellMLBackwardEuler, BidomainProblem<1>, 1, 2> tester;
            RunConvergenceTester(&tester, stimulusType);
            //Result of KSP tester:
            TS_ASSERT_DELTA(HeartConfig::Instance()->GetAbsoluteTolerance(), 1e-3, 1e-10);

            //See above - we've fiddled with HeartConfig...
            HeartConfig::Instance()->SetUseAbsoluteTolerance(5.0e-4);
        }

        {
            std::cout << "OdeConvergenceTester<CellLuoRudy1991FromCellML, BidomainProblem<1>, 1, 2>\n";
            OdeConvergenceTester<CellLuoRudy1991FromCellML, BidomainProblem<1>, 1, 2> tester;
            RunConvergenceTester(&tester, stimulusType);
            TS_ASSERT_DELTA(tester.OdeTimeStep, 0.0025, 1e-10);
        }

        {
            std::cout << "OdeConvergenceTester<CellLuoRudy1991FromCellMLBackwardEuler, BidomainProblem<1>, 1, 2>\n";
            OdeConvergenceTester<CellLuoRudy1991FromCellMLBackwardEuler, BidomainProblem<1>, 1, 2> tester;
            tester.PdeTimeStep=0.01;
            RunConvergenceTester(&tester, stimulusType);
            TS_ASSERT_DELTA(tester.OdeTimeStep, 0.0025, 1e-10);
        }
        //Put the KSP defaults back (see above)
        TS_ASSERT(strcmp(HeartConfig::Instance()->GetKSPPreconditioner(), "jacobi")==0);
        TS_ASSERT(strcmp(HeartConfig::Instance()->GetKSPSolver(), "gmres")==0);
        HeartConfig::Instance()->SetKSPPreconditioner("bjacobi");
        HeartConfig::Instance()->SetKSPSolver("cg");
    }

public:

    void TestStimulatePlanein1D()
    {
        ConvergeInVarious(PLANE);
    }

    void TestStimulateRegionin1D()
    {
        ConvergeInVarious(QUARTER);
    }


    void TestFullActionPotential()
    {
        SpaceConvergenceTester<CellLuoRudy1991FromCellMLBackwardEuler, BidomainProblem<1>, 1, 2> tester;
        tester.SimulateFullActionPotential=true;
        //Time steps are okay for giving a sensible upstroke
        tester.PdeTimeStep=0.1;
        tester.OdeTimeStep=0.1;

        tester.Converge(__FUNCTION__);
        TS_ASSERT(tester.IsConverged());

        ///Note that long plateau phase will force convergence to happen earlier (mesh 3 instead of 4)
        TS_ASSERT_EQUALS(tester.MeshNum, 3u);

        TS_ASSERT_DELTA(329.5, tester.Apd90FirstQn, 1.5); //330.8
        TS_ASSERT_DELTA(329.5, tester.Apd90ThirdQn, 1.5); //329.3
        TS_ASSERT_DELTA(0.0588, tester.ConductionVelocity, 1e-3);
    }

    void TestFullActionPotentialWithRampedStimulus()
    {
        SpaceConvergenceTester<CellLuoRudy1991FromCellMLBackwardEuler, BidomainProblem<1>, 1, 2> tester;
        tester.SimulateFullActionPotential=true;
        //Time steps are okay for giving a sensible upstroke
        tester.PdeTimeStep=0.1;
        tester.OdeTimeStep=0.1;
        tester.Stimulus = QUARTER;

        tester.Converge(__FUNCTION__);
        TS_ASSERT(tester.IsConverged());

        ///Note that long plateau phase will force convergence to happen earlier (mesh 2 instead of mesh 4)
        TS_ASSERT_EQUALS(tester.MeshNum, 2u);

        TS_ASSERT_DELTA(329.5, tester.Apd90FirstQn, 1.5); //329.5
        TS_ASSERT_DELTA(329.5, tester.Apd90ThirdQn, 1.5); //328.7
        TS_ASSERT_DELTA(0.0588, tester.ConductionVelocity, 1e-3);
    }

    //Current test takes about 20 mins.
    //This is much longer (1 hour?) with default ksp
    void Test2DSpaceSymmLq()
    {
        HeartConfig::Instance()->SetKSPSolver("symmlq");
        HeartConfig::Instance()->SetKSPPreconditioner("bjacobi");
        SpaceConvergenceTester<CellLuoRudy1991FromCellMLBackwardEuler, BidomainProblem<2>, 2, 2> tester;
        //tester.SetKspAbsoluteTolerance(1e-3);
         HeartConfig::Instance()->SetUseAbsoluteTolerance(1e-3);

        tester.Converge(__FUNCTION__);
        TS_ASSERT(tester.Converged);
        TS_ASSERT_EQUALS(tester.MeshNum, 5u);
        TS_ASSERT(tester.IsConverged());
        //TS_ASSERT_EQUALS(tester.GetMeshNum(), 4);
        //TS_ASSERT_DELTA(tester.GetSpaceStep(), 0.0023 /*cm*/, 1e-4 /*Allowed error*/);
        HeartConfig::Instance()->Reset();
    }

    void Test2DSpaceWithRegionStimulus()
    {
        HeartConfig::Instance()->SetKSPSolver("symmlq");
        HeartConfig::Instance()->SetKSPPreconditioner("bjacobi");
        SpaceConvergenceTester<CellLuoRudy1991FromCellMLBackwardEuler, BidomainProblem<2>, 2, 2> tester;
        tester.Stimulus = QUARTER;
        HeartConfig::Instance()->SetUseAbsoluteTolerance(1e-3);

        tester.Converge(__FUNCTION__);
        TS_ASSERT(tester.Converged);
        TS_ASSERT_EQUALS(tester.MeshNum, 4u);
        HeartConfig::Instance()->Reset();
    }

    //Currently takes about 15 seconds to do mesh0, mesh1 and mesh2
    void Test3DSpace()
    {
        HeartConfig::Instance()->SetKSPSolver("symmlq");
        HeartConfig::Instance()->SetKSPPreconditioner("bjacobi");
        SpaceConvergenceTester<CellLuoRudy1991FromCellMLBackwardEuler, BidomainProblem<3>, 3, 2> tester;
        HeartConfig::Instance()->SetUseAbsoluteTolerance(1e-3);

        tester.RelativeConvergenceCriterion=2e-1;//Just to prove the thing works
        tester.Converge(__FUNCTION__);
        TS_ASSERT(tester.Converged);
        TS_ASSERT_EQUALS(tester.MeshNum, 2u); ///Just to prove the thing works
        HeartConfig::Instance()->Reset();
    }

    void TestSpaceConvergencein1DWithBackwardN98()
    {
        SpaceConvergenceTester<CellNobleVargheseKohlNoble1998aFromCellMLOpt,  MonodomainProblem<1>, 1, 1> tester;
        tester.AbsoluteStimulus = -5e6; // The default of -1e7 causes V to go out of range for lookup tables
        tester.Converge(__FUNCTION__);
        TS_ASSERT(tester.Converged);
        TS_ASSERT_EQUALS(tester.MeshNum, 5u);
        TS_ASSERT_LESS_THAN(tester.LastDifference, 0.0041039);
    }

    void TestOdeConvergencein1DWithBackwardN98()
    {
        OdeConvergenceTester<CellNobleVargheseKohlNoble1998aFromCellMLBackwardEuler,  MonodomainProblem<1>, 1, 1> tester;
        tester.AbsoluteStimulus = -5e6; // The default of -1e7 causes V to go out of range for lookup tables
        tester.Converge(__FUNCTION__);
        TS_ASSERT(tester.Converged);
        TS_ASSERT_DELTA(tester.OdeTimeStep, 0.0025, 1e-10);
    }

    void TestOdePdeConvergencein1DWithBackwardN98()
    {
        OdePdeConvergenceTester<CellNobleVargheseKohlNoble1998aFromCellMLBackwardEuler,  MonodomainProblem<1>, 1, 1> tester;
        tester.NeumannStimulus = 5000;
        tester.Stimulus = NEUMANN;
        tester.Converge(__FUNCTION__);
        TS_ASSERT(tester.Converged);
        TS_ASSERT_DELTA(tester.OdeTimeStep, 0.005, 1e-10);
        TS_ASSERT_DELTA(tester.PdeTimeStep, 0.005, 1e-10);
    }

    void TestOdePdeConvergencein1DWithForwardLookupN98()
    {
        OdePdeConvergenceTester<CellNobleVargheseKohlNoble1998aFromCellMLOpt,  MonodomainProblem<1>, 1, 1> tester;
        tester.NeumannStimulus = 5000;
        tester.Stimulus = NEUMANN;
        tester.Converge(__FUNCTION__);
        TS_ASSERT(tester.Converged);
        TS_ASSERT_DELTA(tester.OdeTimeStep, 0.0025, 1e-10);
        TS_ASSERT_DELTA(tester.PdeTimeStep, 0.0025, 1e-10);
    }
    void TestOdePdeConvergencein1DWithForwardBasicN98()
    {
        OdePdeConvergenceTester<CellNobleVargheseKohlNoble1998aFromCellML,  MonodomainProblem<1>, 1, 1> tester;
        tester.NeumannStimulus = 5000;
        tester.Stimulus = NEUMANN;
        tester.Converge(__FUNCTION__);
        TS_ASSERT(tester.Converged);
        TS_ASSERT_DELTA(tester.OdeTimeStep, 0.0025, 1e-10);
        TS_ASSERT_DELTA(tester.PdeTimeStep, 0.0025, 1e-10);
    }
};

#endif /*TESTCONVERGENCENIGHTLY_HPP_*/
