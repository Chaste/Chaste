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

#ifndef TESTHEARTFILEFINDER_HPP_
#define TESTHEARTFILEFINDER_HPP_

#include <cxxtest/TestSuite.h>
#include "HeartFileFinder.hpp"
#include "ChasteBuildRoot.hpp"
#include "OutputFileHandler.hpp"
#include "HeartConfig.hpp"

class TestHeartFileFinder : public CxxTest::TestSuite
{
public:
    void TestHeartFileFinderOpening() throw(Exception)
    {
        {
            // Can we find our own source file?
            std::string file_name = "heart/src/io/HeartFileFinder.hpp";
            cp::path_type path(file_name);
            path.relative_to(cp::relative_to_type::chaste_source_root);
            HeartFileFinder file_finder(path);
            TS_ASSERT(file_finder.Exists());
            // Check the path is as expected
            std::string abs_path = ChasteBuildRootDir() + file_name;
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
            out_stream fp = handler.OpenOutputFile(file_name);
            fp->close();
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
            TS_ASSERT_EQUALS(sibling.GetAbsolutePath(), ChasteBuildRootDir() + file_name);
        }
    }
};

#endif /*TESTHEARTFILEFINDER_HPP_*/


