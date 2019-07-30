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


#ifndef TEST3DPERFORMANCE_HPP_
#define TEST3DPERFORMANCE_HPP_

#include <cxxtest/TestSuite.h>
#include "BidomainProblem.hpp"
#include "MonodomainProblem.hpp"
#include <petscvec.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <math.h>

#include "LuoRudy1991BackwardEuler.hpp"
#include "LuoRudy1991.hpp"
#include "PerformanceTester.hpp"

class Test3DPerformance : public CxxTest::TestSuite
{
public:
    void TestPerf()
    {
        HeartConfig::Instance()->SetKSPSolver("symmlq");
        HeartConfig::Instance()->SetKSPPreconditioner("bjacobi");
        PetscTools::SetOption("-log_summary", "");
        // write headings
        PerformanceTester<CellLuoRudy1991FromCellMLBackwardEuler, BidomainProblem<3>, 3>::DisplayHeadings();
        HeartEventHandler::Headings();

        // base line test
        PerformanceTester<CellLuoRudy1991FromCellMLBackwardEuler, BidomainProblem<3>, 3> tester("Test3DPerf");
        tester.MeshNum=1;
        tester.SimTime=4.0;
        tester.Run();
        HeartEventHandler::Report();

        // vary simulation time
        tester.SimTime=0.0025;
        tester.Run();
        HeartEventHandler::Report();

        tester.SimTime=8.0;
        tester.Run();
        HeartEventHandler::Report();

        tester.SimTime=4.0;

        // vary pde time step
        tester.PdeTimeStep=0.005;
        tester.Run();
        HeartEventHandler::Report();

        tester.PdeTimeStep=0.01;
        tester.Run();
        HeartEventHandler::Report();

        tester.PdeTimeStep=0.0025;
        // vary ode time step
        tester.OdeTimeStep=0.0025/2;
        tester.Run();
        HeartEventHandler::Report();

        tester.OdeTimeStep=0.0025/4;
        tester.Run();
        HeartEventHandler::Report();

        tester.OdeTimeStep=0.0025;

        // vary printing time step
        tester.PrintingTimeStep = 0.02;
        tester.Run();
        HeartEventHandler::Report();

        tester.PrintingTimeStep = 0.01;
        tester.Run();
        HeartEventHandler::Report();

        tester.PrintingTimeStep = 0.04;

        // vary mesh size
        tester.MeshNum++;
        tester.Run();
        HeartEventHandler::Report();

        tester.MeshNum++;
        tester.Run();
        HeartEventHandler::Report();

        tester.MeshNum++;
        tester.Run();
        HeartEventHandler::Report();
///\todo: Can't run mesh 5 yet: runs out of memory.
//        tester.MeshNum++;
//        tester.Run();
//        HeartEventHandler::Report();
    }
};

#endif /*TEST3DPERFORMANCE_HPP_*/
