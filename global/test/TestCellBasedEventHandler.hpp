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

    void TestEvents() throw(Exception)
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
