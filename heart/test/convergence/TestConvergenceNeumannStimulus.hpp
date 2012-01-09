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


#ifndef TESTCONVERGENCENEUMANNSTIMULUS_HPP_
#define TESTCONVERGENCENEUMANNSTIMULUS_HPP_

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
//#include "StimulusConvergenceTester.hpp"

#include "PetscSetupAndFinalize.hpp"


class TestConvergenceNeumannStimulus : public CxxTest::TestSuite
{
public:

    void TestConvergenceMonodomain1d()
    {
        PdeConvergenceTester<CellLuoRudy1991FromCellMLBackwardEuler, MonodomainProblem<1>, 1, 1> tester;
        tester.Stimulus=NEUMANN;
        tester.Converge(__FUNCTION__);
        TS_ASSERT(tester.Converged);
        TS_ASSERT_DELTA(tester.PdeTimeStep, 5.0e-3, 1e-10);
    }

    void TestSpaceConvergence1d()
    {
        SpaceConvergenceTester<CellLuoRudy1991FromCellMLBackwardEuler, MonodomainProblem<1>, 1, 1> tester;
        tester.Stimulus=NEUMANN;
        tester.Converge(__FUNCTION__);
        TS_ASSERT(tester.Converged);
        TS_ASSERT_EQUALS(tester.MeshNum, 5u);
    }

    void TestSpaceConvergence2d()
    {
        SpaceConvergenceTester<CellLuoRudy1991FromCellMLBackwardEuler, MonodomainProblem<2>, 2, 1> tester;
        tester.Stimulus=NEUMANN;
        tester.Converge(__FUNCTION__);
        TS_ASSERT(tester.Converged);
        TS_ASSERT_EQUALS(tester.MeshNum, 5u);
    }


};

#endif /*TESTCONVERGENCENEUMANNSTIMULUS_HPP_*/
