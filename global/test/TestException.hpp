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

#ifndef _TESTEXCEPTION_HPP_
#define _TESTEXCEPTION_HPP_

#include <cxxtest/TestSuite.h>
#include "Exception.hpp"
#include "PetscSetupAndFinalize.hpp"

class TestException : public CxxTest::TestSuite
{
public:

    void TestGetMessage()
    {
        std::string msg("This is an exception");

        try
        {
            EXCEPTION(msg);
        }
        catch (Exception& exc)
        {
            std::string e_msg = exc.GetMessage();
            std::string::size_type e_len = e_msg.length();
            std::string::size_type len = msg.length();
            TS_ASSERT_EQUALS(e_msg.substr(e_len - len), msg);
        }

        TS_ASSERT_THROWS_EQUALS(EXCEPTION("Hello. I'm an exception"), const Exception &err,
                                err.GetShortMessage(), "Hello. I'm an exception" );
        /*
         * Note: The following test will fail if the number of lines above is changed drastically
         * (that's why method GetShortMessage() was introduced).
         */
#ifndef _MSC_VER
        TS_ASSERT_THROWS_EQUALS(EXCEPTION("Hello. I'm an exception"), const Exception &err,
                                err.GetMessage().find("Hello. I\'m an exception",0), 51u); // This appears at position 51 in full message (a bit more robust?!)
#endif // _MSC_VER
    }


    void TestCheckStreamingException()
    {
        int rule=42;
        try
        {
            EXCEPTION("This is a rule "<<rule<<" exception");
        }
        catch (const Exception& exc)
        {
            TS_ASSERT_EQUALS(exc.CheckShortMessage("This is a rule 42 exception"), "");
        }
    }

    void TestCheckMethods()
    {
        std::string msg("This is our message");

        // If another process threw an exception, this is what this process will actually throw
        std::string parallel_msg("Another process threw an exception; bailing out.");

        try
        {
            EXCEPTION(msg);
        }
        catch (const Exception& e)
        {
            TS_ASSERT_EQUALS(e.CheckShortMessage(msg), "");
            TS_ASSERT_DIFFERS(e.CheckShortMessage("This is not our message"), "");
        }

        TS_ASSERT_THROWS_THIS(EXCEPTION(msg), msg);
        TS_ASSERT_THROWS_THIS(EXCEPTION(parallel_msg), msg);
        TS_ASSERT_THROWS_THIS(EXCEPTION(parallel_msg), "Anything you like will work here");

        try
        {
            EXCEPTION(parallel_msg);
        }
        catch (const Exception& e)
        {
            TS_ASSERT_EQUALS(e.CheckShortMessage(msg), "");
            // Of course, in this case we don't know what the real message would have been...
            TS_ASSERT_EQUALS(e.CheckShortMessage("This is not our message"), "");
        }

        // Now check for messages containing a string
        try
        {
            EXCEPTION(msg + " extra bit");
        }
        catch (const Exception& e)
        {
            TS_ASSERT_EQUALS(e.CheckShortMessageContains(msg), "");
            TS_ASSERT_DIFFERS(e.CheckShortMessageContains("This is not our message"), "");
        }

        TS_ASSERT_THROWS_CONTAINS(EXCEPTION(msg + " extra bit"), msg);
        TS_ASSERT_THROWS_CONTAINS(EXCEPTION("Prefix " + msg + " extra bit"), msg);
        TS_ASSERT_THROWS_CONTAINS(EXCEPTION(msg), msg);
        TS_ASSERT_THROWS_CONTAINS(EXCEPTION(parallel_msg), msg);
        TS_ASSERT_THROWS_CONTAINS(EXCEPTION(parallel_msg), "Anything you like will work here");

        try
        {
            EXCEPTION(parallel_msg);
        }
        catch (const Exception& e)
        {
            TS_ASSERT_EQUALS(e.CheckShortMessageContains(msg), "");
            // Of course, in this case we don't know what the real message would have been...
            TS_ASSERT_EQUALS(e.CheckShortMessageContains("This is not our message"), "");
        }

        TS_ASSERT_THROWS_THIS(EXCEPT_IF_NOT(false), "Assertion tripped: false");
        TS_ASSERT_THROWS_NOTHING(EXCEPT_IF_NOT(true));
    }
};

#endif //_TESTEXCEPTION_HPP_
