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

#ifndef _TESTCONSTBOUNDARYCONDITION_HPP_
#define _TESTCONSTBOUNDARYCONDITION_HPP_

#include <cxxtest/TestSuite.h>
#include "ConstBoundaryCondition.hpp"

class TestConstDirichletBoundaryCondition : public CxxTest::TestSuite
{
public:
    void TestConstDirichletBoundaryConditionMethod()
    {
        ChastePoint<1> zero(0);
        ConstBoundaryCondition<1> const_dirich1(1.0);
        double value = const_dirich1.GetValue(zero);
        TS_ASSERT_DELTA(value,1.0,1e-12);
    }
};

#endif //_TESTCONSTBOUNDARYCONDITION_HPP_
