/*

Copyright (c) 2005-2012, University of Oxford.
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

#ifndef TESTCHASTEBUILDINFO_HPP_
#define TESTCHASTEBUILDINFO_HPP_

#include <cxxtest/TestSuite.h>

#include <iostream>
#include <boost/foreach.hpp>
#include <hdf5.h>
#include <parmetis.h>

#include "Exception.hpp"
#include "PetscTools.hpp"
#include "PetscException.hpp"
#include "Version.hpp"
#include "ChasteSerialization.hpp"

#ifdef CHASTE_VTK
#define _BACKWARD_BACKWARD_WARNING_H 1 //Cut out the strstream deprecated warning for now (gcc4.3)
#include <vtkVersion.h>
#endif

#ifdef CHASTE_CVODE
#include <sundials/sundials_config.h>
#endif

/**
 * The ChasteBuildInfo class isn't really amenable to testing.
 * We just print out all the information available.
 */
class TestChasteBuildInfo : public CxxTest::TestSuite
{
    typedef std::pair<std::string, std::string> StringPair;
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

        std::cout << "Checked out projects:" << std::endl;
        BOOST_FOREACH(const StringPair& r_project_version, ChasteBuildInfo::rGetProjectVersions())
        {
            std::cout << "\t" << r_project_version.first << ": " << r_project_version.second << std::endl;
        }

        // The following is stolen from the ExecutableSupport() file, to log this info in test results.
        std::cout << std::endl << "Library versions:" << std::endl;

        std::cout << "<CompiledIn>" << std::endl;
        std::cout << "\t<PETSc>" << PETSC_VERSION_MAJOR << "." << PETSC_VERSION_MINOR << "." << PETSC_VERSION_SUBMINOR << "</PETSc>" << std::endl;
        std::cout << "\t<Boost>" << BOOST_VERSION  / 100000 << "." << BOOST_VERSION / 100 % 1000 << "." << BOOST_VERSION % 100 << "</Boost>" << std::endl;
        std::cout << "\t<HDF5>" << H5_VERS_MAJOR <<  "." << H5_VERS_MINOR << "." << H5_VERS_RELEASE << "</HDF5>" << std::endl;
        std::cout << "\t<Parmetis>" << PARMETIS_MAJOR_VERSION << "." << PARMETIS_MINOR_VERSION;
    #ifdef PARMETIS_SUBMINOR_VERSION // they only added this in v4.? !!
        std::cout << "." << PARMETIS_SUBMINOR_VERSION;
    #endif
        std::cout << "</Parmetis>" << std::endl;
        std::cout << "</CompiledIn>" << std::endl;

        std::cout << "<Binaries>" << std::endl;
        std::cout << "\t<XSD>" <<  ChasteBuildInfo::GetXsdVersion() << "</XSD>" << std::endl;
        std::cout << "</Binaries>" << std::endl;

        std::cout << "<Optional>" << std::endl;
    #ifdef CHASTE_VTK
        std::cout << "\t<VTK>" << VTK_MAJOR_VERSION << "." << VTK_MINOR_VERSION << "</VTK>" << std::endl;
    #else
        std::cout << "\t<VTK>no</VTK>" << std::endl;
    #endif

    #ifdef CHASTE_CVODE
        std::cout << "\t<SUNDIALS>" << SUNDIALS_PACKAGE_VERSION << "</SUNDIALS> <!-- includes Cvode of a different version number! --> " << std::endl;
    #else
        std::cout << "\t<SUNDIALS>no</SUNDIALS>" << std::endl;
    #endif

    #ifdef CHASTE_ADAPTIVITY
        std::cout << "\t<Adaptivity>yes</Adaptivity>" << std::endl;
    #else
        std::cout << "\t<Adaptivity>no</Adaptivity>" << std::endl;
    #endif
        std::cout << "</Optional>" << std::endl << std::flush;

    }
};

#endif /* TESTCHASTEBUILDINFO_HPP_ */
