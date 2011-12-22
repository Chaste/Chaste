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


#ifndef TESTPERFORMANCE_HPP_
#define TESTPERFORMANCE_HPP_

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

class TestPerformance : public CxxTest::TestSuite
{
public:
    void TestPerf() throw(Exception)
    {
        HeartConfig::Instance()->SetKSPSolver("symmlq");
        HeartConfig::Instance()->SetKSPPreconditioner("bjacobi");
        PetscOptionsSetValue("-log_summary", "");
        // write headings
        PerformanceTester<CellLuoRudy1991FromCellMLBackwardEuler, BidomainProblem<3>, 3>::DisplayHeadings();
        HeartEventHandler::Headings();

        // base line test
        PerformanceTester<CellLuoRudy1991FromCellMLBackwardEuler, BidomainProblem<3>, 3> tester;
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

#endif /*TESTPERFORMANCE_HPP_*/
