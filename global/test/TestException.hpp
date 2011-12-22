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

#ifndef _TESTEXCEPTION_HPP_
#define _TESTEXCEPTION_HPP_

#include <cxxtest/TestSuite.h>
#include "Exception.hpp"

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
        catch (Exception e)
        {
            std::string e_msg = e.GetMessage();
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
        TS_ASSERT_THROWS_EQUALS(EXCEPTION("Hello. I'm an exception"), const Exception &err,
                                err.GetMessage().find("Hello. I\'m an exception",0), 51u); // This appears at position 51 in full message (a bit more robust?!)
    }


    void TestCheckStreamingException()
    {
        int rule=42;
        try
        {
            EXCEPTION("This is a rule "<<rule<<" exception");
        }
        catch (const Exception& e)
        {
            TS_ASSERT_EQUALS(e.CheckShortMessage("This is a rule 42 exception"), "");
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
