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
#ifndef TESTCELLBASEDEVENTHANDLER_HPP_
#define TESTCELLBASEDEVENTHANDLER_HPP_

#include "PetscSetupAndFinalize.hpp"
#include "CellBasedEventHandler.hpp"

/**
 * This class consists of a single test for the CellBasedEventHandler
 * class.
 */
class TestCellBasedEventHandler : public CxxTest::TestSuite
{
public:

    void TestEvents()
    {
        CellBasedEventHandler::BeginEvent(CellBasedEventHandler::EVERYTHING);

        CellBasedEventHandler::BeginEvent(CellBasedEventHandler::SETUP);
        CellBasedEventHandler::MilliSleep(10);

        CellBasedEventHandler::EndEvent(CellBasedEventHandler::SETUP);

        CellBasedEventHandler::BeginEvent(CellBasedEventHandler::DEATH);
        CellBasedEventHandler::MilliSleep(20);
        CellBasedEventHandler::EndEvent(CellBasedEventHandler::DEATH);

        CellBasedEventHandler::BeginEvent(CellBasedEventHandler::BIRTH);
        CellBasedEventHandler::MilliSleep(30);
        CellBasedEventHandler::EndEvent(CellBasedEventHandler::BIRTH);
        CellBasedEventHandler::BeginEvent(CellBasedEventHandler::UPDATECELLPOPULATION);
        CellBasedEventHandler::MilliSleep(40);
        CellBasedEventHandler::EndEvent(CellBasedEventHandler::UPDATECELLPOPULATION);

        CellBasedEventHandler::BeginEvent(CellBasedEventHandler::TESSELLATION);
        CellBasedEventHandler::MilliSleep(50);
        CellBasedEventHandler::EndEvent(CellBasedEventHandler::TESSELLATION);

        CellBasedEventHandler::BeginEvent(CellBasedEventHandler::FORCE);
        CellBasedEventHandler::MilliSleep(60);
        CellBasedEventHandler::EndEvent(CellBasedEventHandler::FORCE);

        CellBasedEventHandler::BeginEvent(CellBasedEventHandler::POSITION);
        CellBasedEventHandler::MilliSleep(70);
        CellBasedEventHandler::EndEvent(CellBasedEventHandler::POSITION);

        CellBasedEventHandler::BeginEvent(CellBasedEventHandler::OUTPUT);
        CellBasedEventHandler::MilliSleep(80);
        CellBasedEventHandler::EndEvent(CellBasedEventHandler::OUTPUT);

        CellBasedEventHandler::BeginEvent(CellBasedEventHandler::PDE);
        CellBasedEventHandler::MilliSleep(90);
        CellBasedEventHandler::EndEvent(CellBasedEventHandler::PDE);

        CellBasedEventHandler::EndEvent(CellBasedEventHandler::EVERYTHING);

        CellBasedEventHandler::Headings();

        CellBasedEventHandler::Report();
    }
};

#endif /*TESTCELLBASEDEVENTHANDLER_HPP_*/
