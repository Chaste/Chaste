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

#ifndef TESTHEARTEVENTHANDLER_HPP_
#define TESTHEARTEVENTHANDLER_HPP_

#include "HeartEventHandler.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "PetscTools.hpp"

class TestHeartEventHandler : public CxxTest::TestSuite
{
public:

    void TestEvents()
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

    void TestParallelPrinting()
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

    void TestEventExceptions()
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
