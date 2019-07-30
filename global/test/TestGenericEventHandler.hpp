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

#ifndef TESTGENERICEVENTHANDLER_HPP_
#define TESTGENERICEVENTHANDLER_HPP_

#include "GenericEventHandler.hpp"
#include "PetscSetupAndFinalize.hpp"

class AnEventHandler : public GenericEventHandler<3, AnEventHandler>
{
public:
    static const char* EventName[3];

    typedef enum
    {
        TEST1=0,
        TEST2,
        TEST3
    } EventType;
};

const char* AnEventHandler::EventName[] = { "Test1", "Test2", "Test3"};

class TestGenericEventHandler : public CxxTest::TestSuite
{
public:

    void TestEvents()
    {
        // Coverage
        AnEventHandler::Instance()->DisableImpl();
        TS_ASSERT_THROWS_THIS(AnEventHandler::Instance()->GetElapsedTime(AnEventHandler::TEST1),
                "Asked to report on a disabled event handler.  Check for contributory errors above.");
        AnEventHandler::Instance()->EnableImpl();

        AnEventHandler::BeginEvent(AnEventHandler::TEST1);
        // The first BeginEvent implicitly calls:
        // AnEventHandler::BeginEvent(AnEventHandler::TEST3);

        AnEventHandler::BeginEvent(AnEventHandler::TEST2);
        AnEventHandler::MilliSleep(10);
        AnEventHandler::EndEvent(AnEventHandler::TEST2);

        AnEventHandler::MilliSleep(10);
        AnEventHandler::EndEvent(AnEventHandler::TEST3);

        AnEventHandler::EndEvent(AnEventHandler::TEST1);

        AnEventHandler::Headings();

        AnEventHandler::Report();

        // No longer allowed to report twice
        TS_ASSERT_THROWS_THIS(AnEventHandler::Report(),"Asked to report on an event handler which is set to zero.");
    }

    void TestEventExceptions()
    {
        // Should not be able to end an event that has not yet begun
        TS_ASSERT_THROWS_THIS(AnEventHandler::EndEvent(AnEventHandler::TEST1),
                "Error: The event associated with the counter for \'Test1\' had not begun when EndEvent was called.");

        AnEventHandler::BeginEvent(AnEventHandler::TEST1);

        // Beginning an event already begun should print an error message and disable the handler
        AnEventHandler::BeginEvent(AnEventHandler::TEST1);
        TS_ASSERT(!AnEventHandler::IsEnabled());

        // Report should then throw
        TS_ASSERT_THROWS_THIS(AnEventHandler::Report(),
                "Asked to report on a disabled event handler.  Check for contributory errors above.");
    }

    void TestReset()
    {
        // Clear up from previous test
        AnEventHandler::Reset();
        AnEventHandler::BeginEvent(AnEventHandler::TEST1);
        AnEventHandler::BeginEvent(AnEventHandler::TEST2);
        AnEventHandler::Reset();

        // One can now being these events again because the state of the event handler was reset
        AnEventHandler::BeginEvent(AnEventHandler::TEST1);
        AnEventHandler::BeginEvent(AnEventHandler::TEST2);
    }

    void TestDisable()
    {
        AnEventHandler::Reset();
        AnEventHandler::Disable();
        AnEventHandler::BeginEvent(AnEventHandler::TEST1);
        AnEventHandler::BeginEvent(AnEventHandler::TEST1); // OK because event handling is disabled
        AnEventHandler::Enable();
    }

    void TestElapsedTime()
    {
        AnEventHandler::Reset();
        TS_ASSERT_EQUALS(AnEventHandler::GetElapsedTime(AnEventHandler::TEST1), 0.0);
        TS_ASSERT_EQUALS(AnEventHandler::GetElapsedTime(AnEventHandler::TEST2), 0.0);
        TS_ASSERT_EQUALS(AnEventHandler::GetElapsedTime(AnEventHandler::TEST3), 0.0);

        AnEventHandler::BeginEvent(AnEventHandler::TEST1);
        long dummy = 1;
        for (unsigned i=0; i<1e9; i++)
        {
            dummy += 2;
        }
        std::cout << "Printing variable to slow test down!: dummy = " << dummy << std::endl;
        TS_ASSERT_EQUALS(2000000001l, dummy); // try to avoid the loop being optimised away
        TS_ASSERT_LESS_THAN(0.0, AnEventHandler::GetElapsedTime(AnEventHandler::TEST1));
        AnEventHandler::EndEvent(AnEventHandler::TEST1);
        TS_ASSERT_LESS_THAN(0.0, AnEventHandler::GetElapsedTime(AnEventHandler::TEST1));

        AnEventHandler::BeginEvent(AnEventHandler::TEST2);
        AnEventHandler::MilliSleep(11);
        dummy = 0; // Separate the sleep from the end of the event
        AnEventHandler::EndEvent(AnEventHandler::TEST2);

        // Test in milliseconds (at least 10 and not too much)
        TS_ASSERT_LESS_THAN_EQUALS(10.0, AnEventHandler::GetElapsedTime(AnEventHandler::TEST2));
        TS_ASSERT_LESS_THAN_EQUALS(AnEventHandler::GetElapsedTime(AnEventHandler::TEST2), 60.0);
    }

    void TestSilentlyCloseEvent()
    {
        AnEventHandler::Headings();
        AnEventHandler::Reset();
        AnEventHandler::Enable();
        AnEventHandler::BeginEvent(AnEventHandler::TEST1);
        AnEventHandler::MilliSleep(10);

        AnEventHandler::Report();
    }

    void TestReportPrecision()
    {
        AnEventHandler::Headings();
        AnEventHandler::Reset();
        AnEventHandler::Enable();
        AnEventHandler::BeginEvent(AnEventHandler::TEST3);      // total time
        AnEventHandler::Instance()->mWallTime[AnEventHandler::TEST3] += 1e3; // fake short run time
        AnEventHandler::EndEvent(AnEventHandler::TEST3);
        AnEventHandler::Report();

        AnEventHandler::BeginEvent(AnEventHandler::TEST3);
        AnEventHandler::Instance()->mWallTime[AnEventHandler::TEST3] += 1e5; // fake longer run time
        AnEventHandler::EndEvent(AnEventHandler::TEST3);
        AnEventHandler::Report();

        AnEventHandler::BeginEvent(AnEventHandler::TEST3);
        AnEventHandler::Instance()->mWallTime[AnEventHandler::TEST3] += 1e7; // fake long run time
        AnEventHandler::EndEvent(AnEventHandler::TEST3);
        AnEventHandler::Report();
    }
};

#endif /*TESTGENERICEVENTHANDLER_HPP_*/
