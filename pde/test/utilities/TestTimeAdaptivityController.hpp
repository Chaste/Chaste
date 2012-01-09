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

#ifndef TESTTIMEADAPTIVITYCONTROLLER_HPP_
#define TESTTIMEADAPTIVITYCONTROLLER_HPP_

#include <cxxtest/TestSuite.h>

#include "AbstractTimeAdaptivityController.hpp"

// Very simple toy time-adaptivity controller
class ToyController : public AbstractTimeAdaptivityController
{
    double ComputeTimeStep(double currentTime, Vec currentSolution)
    {
        if (currentTime < 1)
        {
            return 0.1;
        }
        else if (currentTime < 2)
        {
            return 0.5;
        }
        else
        {
            return 1.5;
        }
    }

public:
    ToyController(double minDt, double maxDt)
        : AbstractTimeAdaptivityController(minDt, maxDt)
    {
    }
};

class TestTimeAdaptivityController : public CxxTest::TestSuite
{
public:
    void TestToyController() throw(Exception)
    {
        ToyController controller(0.2,1.0);
        TS_ASSERT_EQUALS(controller.GetNextTimeStep(0.5,NULL), 0.2);
        TS_ASSERT_EQUALS(controller.GetNextTimeStep(1.5,NULL), 0.5);
        TS_ASSERT_EQUALS(controller.GetNextTimeStep(10 ,NULL), 1.0);
    }
};

#endif /*TESTTIMEADAPTIVITYCONTROLLER_HPP_*/
