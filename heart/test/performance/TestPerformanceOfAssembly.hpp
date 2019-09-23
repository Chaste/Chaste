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


#ifndef TESTPERFORMANCEOFASSEMBLY_HPP_
#define TESTPERFORMANCEOFASSEMBLY_HPP_

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

class TestPerformanceOfAssembly : public CxxTest::TestSuite
{
public:

    void TestPerfAssembly()
    {
        // Write headings
        PerformanceTester<CellLuoRudy1991FromCellMLBackwardEuler, BidomainProblem<2>, 2>::DisplayHeadings();
        HeartEventHandler::Headings();

        // Base line
        PerformanceTester<CellLuoRudy1991FromCellMLBackwardEuler, BidomainProblem<2>, 2> tester("TestPerAssembly");
        tester.SimTime=0.0025;

        for (unsigned mesh_num=0; mesh_num<3; mesh_num++)
        {
            tester.MeshNum=mesh_num;
            tester.Run();
            HeartEventHandler::Report();
        }
    }
};

#endif /*TESTPERFORMANCE_HPP_*/
