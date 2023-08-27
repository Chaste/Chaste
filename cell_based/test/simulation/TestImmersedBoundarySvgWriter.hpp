/*

Copyright (c) 2005-2023, University of Oxford.
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

#ifndef TESTIMMERSEDBOUNDARYSVGWRITER_HPP_
#define TESTIMMERSEDBOUNDARYSVGWRITER_HPP_

// Needed for test framework
#include <cxxtest/TestSuite.h>

// Must be included before other cell_based headers
#include "CellBasedSimulationArchiver.hpp"

// Includes from projects/ImmersedBoundary
#include "ImmersedBoundaryMesh.hpp"
#include "ImmersedBoundarySvgWriter.hpp"
#include "FileComparison.hpp"

// This test is never run in parallel
#include "FakePetscSetup.hpp"

class TestImmersedBoundarySvgWriter : public CxxTest::TestSuite
{
public:

    void TestContructor()
    {
      ImmersedBoundarySvgWriter<2> svgWriter;
      TS_ASSERT_EQUALS(svgWriter.GetSamplingMultiple(), 100u);
      TS_ASSERT_EQUALS(svgWriter.GetSvgSize(), 1600.0);
    }

    void TestGetSetMethods()
    {
      ImmersedBoundarySvgWriter<2> svgWriter;
      TS_ASSERT_EQUALS(svgWriter.GetSamplingMultiple(), 100u);
      TS_ASSERT_EQUALS(svgWriter.GetSvgSize(), 1600.0);
      
      svgWriter.SetSamplingMultiple(10u);
      TS_ASSERT_EQUALS(svgWriter.GetSamplingMultiple(), 10u);
      
      svgWriter.SetSvgSize(1200.0);
      TS_ASSERT_EQUALS(svgWriter.GetSvgSize(), 1200.0);
    }

    void TestWriting()
    {
      
    }

    void TestArchiving()
    {
        // Create a file for archiving
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "ImmersedBoundarySvgWriter.arch";

        // Separate scope to write the archive
        {
            // Initialise a growth modifier and set a non-standard mature target area
            ImmersedBoundarySvgWriter<2> writer;
            writer.SetSamplingMultiple(10u);
            writer.SetSvgSize(1200.0);

            // Create an output archive
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Serialize
            output_arch << writer;
        }

        // Separate scope to read the archive
        {
            ImmersedBoundarySvgWriter<2> writer;

            // Restore the modifier
            std::ifstream ifs(archive_filename.c_str());
            boost::archive::text_iarchive input_arch(ifs);

            input_arch >> writer;
            
            TS_ASSERT_EQUALS(writer.GetSamplingMultiple(), 10u);
            TS_ASSERT_EQUALS(writer.GetSvgSize(), 1200.0);

        }
    }

};

#endif /*TESTIMMERSEDBOUNDARYSVGWRITER_HPP_*/
