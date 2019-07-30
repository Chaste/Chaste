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

#ifndef _TESTSTEPSIZEEXCEPTION_HPP_
#define _TESTSTEPSIZEEXCEPTION_HPP_

#include <cxxtest/TestSuite.h>
#include "StepSizeException.hpp"
#include "PetscSetupAndFinalize.hpp"

class TestStepSizeException : public CxxTest::TestSuite
{
public:

    void TestGetMessage()
    {
        double new_step = 0.015;
        std::string message("This is a step size exception");

        try
        {
            StepSizeException(new_step, message, true);
        }
        catch (StepSizeException* exc)
        {
            const char* msg_char = exc->what();
            std::string exc_msg = std::string(msg_char);
            std::string::size_type e_len = exc_msg.length();
            std::string::size_type len = message.length();

            TS_ASSERT_EQUALS(exc_msg.substr(e_len - len), message);
            TS_ASSERT_EQUALS(exc_msg, message);

            TS_ASSERT_EQUALS(exc->IsTerminal(), true);
            TS_ASSERT_DELTA(exc->GetSuggestedNewStep(), new_step, 1e-6);
        }
    }
};

#endif //_TESTSTEPSIZEEXCEPTION_HPP_
