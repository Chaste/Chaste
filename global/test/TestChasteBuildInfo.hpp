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

#ifndef TESTCHASTEBUILDINFO_HPP_
#define TESTCHASTEBUILDINFO_HPP_

#include <cxxtest/TestSuite.h>

#include <iostream>

#include "Version.hpp"

/**
 * The ChasteBuildInfo class isn't really amenable to testing.
 * We just print out all the information available.
 */
class TestChasteBuildInfo : public CxxTest::TestSuite
{
public:
    void TestShowInfo()
    {
        std::cout << "Information on Chaste build." << std::endl;
        std::cout << "Chaste root directory: " << ChasteBuildInfo::GetRootDir() << std::endl;
        std::cout << "Chaste version (multiple methods): " << ChasteBuildInfo::GetMajorReleaseNumber()
                  << "." << ChasteBuildInfo::GetMinorReleaseNumber()
                  << "." << ChasteBuildInfo::GetRevisionNumber() << std::endl;
        std::cout << "Chaste version (single method): " << ChasteBuildInfo::GetVersionString() << std::endl;
        if (ChasteBuildInfo::IsWorkingCopyModified())
        {
            std::cout << "  Built from a modified working copy." << std::endl;
        }
        std::cout << "Built on " << ChasteBuildInfo::GetBuildTime() << std::endl;
        std::cout << "Current time: " << ChasteBuildInfo::GetCurrentTime() << std::endl;
        std::cout << "Build machine: " << ChasteBuildInfo::GetBuilderUnameInfo() << std::endl;
        std::cout << "Build information: " << ChasteBuildInfo::GetBuildInformation() << std::endl;
        std::cout << "Provenance information: " << ChasteBuildInfo::GetProvenanceString() << std::endl;
    }
};

#endif /* TESTCHASTEBUILDINFO_HPP_ */
