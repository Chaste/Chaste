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

#ifndef TESTHEARTFILEFINDER_HPP_
#define TESTHEARTFILEFINDER_HPP_

#include <cxxtest/TestSuite.h>
#include "HeartFileFinder.hpp"
#include "ChasteBuildRoot.hpp"
#include "OutputFileHandler.hpp"
#include "HeartConfig.hpp"
#include "PetscSetupAndFinalize.hpp"

class TestHeartFileFinder : public CxxTest::TestSuite
{
public:
    void TestHeartFileFinderOpening()
    {
        {
            // Can we find our own source file?
            std::string file_name = "heart/src/io/HeartFileFinder.hpp";
            cp::path_type path(file_name);
            path.relative_to(cp::relative_to_type::chaste_source_root);
            HeartFileFinder file_finder(path);
            TS_ASSERT(file_finder.Exists());
            // Check the path is as expected
            std::string abs_path = ChasteSourceRootDir() + file_name;
            TS_ASSERT_EQUALS(file_finder.GetAbsolutePath(), abs_path);

            // CWD should be the Chaste source root
            path.relative_to(cp::relative_to_type::cwd);
            HeartFileFinder file_finder2(path);
            TS_ASSERT(file_finder2.Exists());
            // Check the path is as expected
            TS_ASSERT_EQUALS(file_finder2.GetAbsolutePath(), abs_path);
        }

        {
            // Now check a file in the output directory
            std::string dir_name = "TestHeartFileFinder";
            OutputFileHandler handler(dir_name);
            std::string file_name = "TestFile";
            cp::path_type path(dir_name + "/" + file_name);
            path.relative_to(cp::relative_to_type::chaste_test_output);
            HeartFileFinder file_finder(path);
            TS_ASSERT(! file_finder.Exists());
            // Check the path is as expected
            std::string abs_path = handler.GetOutputDirectoryFullPath() + file_name;
            TS_ASSERT_EQUALS(file_finder.GetAbsolutePath(), abs_path);
            // Create the file
            PetscTools::Barrier("TestHeartFileFinderOpening-1");
            if (PetscTools::AmMaster())
            {
                out_stream fp = handler.OpenOutputFile(file_name);
                fp->close();
            }
            PetscTools::Barrier("TestHeartFileFinderOpening-2");
            TS_ASSERT(file_finder.Exists());

            // Check when providing an absolute path
            cp::path_type path_abs(abs_path);
            path_abs.relative_to(cp::relative_to_type::absolute);
            HeartFileFinder file_finder2(path_abs);
            TS_ASSERT(file_finder2.Exists());
            TS_ASSERT_EQUALS(file_finder2.GetAbsolutePath(), abs_path);
        }

        {
            // Check we can find a sibling to the XML parameters file
            HeartConfig::Instance()->SetParametersFile("ChasteParameters.xml");
            std::string file_name = "SConstruct";
            cp::path_type sibling_path(file_name);
            sibling_path.relative_to(cp::relative_to_type::this_file);
            HeartFileFinder sibling(sibling_path);
            TS_ASSERT(sibling.IsFile());
            TS_ASSERT_EQUALS(sibling.GetAbsolutePath(), ChasteSourceRootDir() + file_name);
        }
    }
};

#endif /*TESTHEARTFILEFINDER_HPP_*/


