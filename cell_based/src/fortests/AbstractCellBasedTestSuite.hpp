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

#ifndef ABSTRACTCELLBASEDTESTSUITE_HPP_
#define ABSTRACTCELLBASEDTESTSUITE_HPP_

#define CXXTEST_ABORT_TEST_ON_FAIL

#include "SimulationTime.hpp"
#include "RandomNumberGenerator.hpp"
#include "CellPropertyRegistry.hpp"

/**
 * This class provides setUp and tearDown methods that are common to
 * many cell_based test suites.  Such suites may inherit from this class
 * to avoid having to redefine them.
 */
class AbstractCellBasedTestSuite : public CxxTest::TestSuite
{
protected:

    /**
     * Overridden setUp() method. Initialises singleton classes.
     */
    void setUp()
    {
        // The following won't work: it returns from this setup method, but not the test suite
        //EXIT_IF_PARALLEL; // defined in PetscTools

        SimulationTime::Instance()->SetStartTime(0.0);
        RandomNumberGenerator::Instance()->Reseed(0);
        CellPropertyRegistry::Instance()->Clear();
    }

    /**
     * Overridden teardown() method. Clears up singleton classes.
     */
    void tearDown()
    {
        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
    }
};

#endif /*ABSTRACTCELLBASEDTESTSUITE_HPP_*/
