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

#ifndef TESTHEARTEVENTHANDLER_HPP_
#define TESTHEARTEVENTHANDLER_HPP_

#include "HeartEventHandler.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "PetscTools.hpp"

class TestHeartEventHandler : public CxxTest::TestSuite
{
public:

    void TestEvents() throw(Exception)
    {
        HeartEventHandler::BeginEvent(HeartEventHandler::EVERYTHING);
        HeartEventHandler::BeginEvent(HeartEventHandler::SOLVE_ODES);
        HeartEventHandler::MilliSleep(10);
        HeartEventHandler::EndEvent(HeartEventHandler::SOLVE_ODES);

        HeartEventHandler::BeginEvent(HeartEventHandler::READ_MESH);
        HeartEventHandler::MilliSleep(10);
        HeartEventHandler::EndEvent(HeartEventHandler::READ_MESH);

        HeartEventHandler::BeginEvent(HeartEventHandler::COMMUNICATION);
        HeartEventHandler::MilliSleep(10);

        HeartEventHandler::BeginEvent(HeartEventHandler::SOLVE_LINEAR_SYSTEM);
        HeartEventHandler::MilliSleep(10);
        HeartEventHandler::EndEvent(HeartEventHandler::SOLVE_LINEAR_SYSTEM);

        HeartEventHandler::EndEvent(HeartEventHandler::COMMUNICATION);
        HeartEventHandler::EndEvent(HeartEventHandler::EVERYTHING);

        HeartEventHandler::Headings();

        HeartEventHandler::Report();
    }

    void TestParallelPrinting() throw (Exception)
    {
        std::cout.flush();
        std::cerr.flush();
        PetscTools::Barrier("TestParallelPrinting");
        std::cout.flush();
        std::cerr.flush();

        HeartEventHandler::BeginEvent(HeartEventHandler::EVERYTHING);
        HeartEventHandler::BeginEvent(HeartEventHandler::READ_MESH);
        if (!PetscTools::AmTopMost())
        {
            HeartEventHandler::MilliSleep(50);
        }
        else
        {
            // Top process has smaller amount of work
            HeartEventHandler::MilliSleep(10);
        }
        HeartEventHandler::EndEvent(HeartEventHandler::READ_MESH);

        HeartEventHandler::EndEvent(HeartEventHandler::EVERYTHING);

        HeartEventHandler::Headings();

        HeartEventHandler::Report();

    }

    void TestEventExceptions() throw(Exception)
    {
        // Should not be able to end and event that has not yet begun
        TS_ASSERT_THROWS_THIS(HeartEventHandler::EndEvent(HeartEventHandler::EVERYTHING),
                "Error: The event associated with the counter for \'Total\' had not begun when EndEvent was called.");

        HeartEventHandler::BeginEvent(HeartEventHandler::EVERYTHING);

        // Beginning an event already begun should print an error message and disable the handler
        HeartEventHandler::BeginEvent(HeartEventHandler::EVERYTHING);
        TS_ASSERT(!HeartEventHandler::IsEnabled());
        // Report should then throw
        TS_ASSERT_THROWS_THIS(HeartEventHandler::Report(),
                "Asked to report on a disabled event handler.  Check for contributory errors above.");
    }
};

#endif /*TESTHEARTEVENTHANDLER_HPP_*/
