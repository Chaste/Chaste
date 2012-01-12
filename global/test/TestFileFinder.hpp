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

#ifndef TESTFILEFINDER_HPP_
#define TESTFILEFINDER_HPP_

#include <cxxtest/TestSuite.h>
#include "FileFinder.hpp"
#include "ChasteBuildRoot.hpp"
#include "OutputFileHandler.hpp"
#include "GetCurrentWorkingDirectory.hpp"

class TestFileFinder : public CxxTest::TestSuite
{
public:

    void TestFileFinderOpening() throw(Exception)
    {
        {
            // Can we find our own source file?
            std::string file_name = "global/src/FileFinder.hpp";
            FileFinder file_finder(file_name, RelativeTo::ChasteSourceRoot);
            TS_ASSERT(file_finder.Exists());
            TS_ASSERT(file_finder.IsFile());
            TS_ASSERT(!file_finder.IsDir());

            // Check the path is as expected
            std::string abs_path = ChasteBuildRootDir() + file_name;
            TS_ASSERT_EQUALS(file_finder.GetAbsolutePath(), abs_path);

            // CWD should be the Chaste source root
            FileFinder file_finder2(file_name, RelativeTo::CWD);
            TS_ASSERT(file_finder2.Exists());
            TS_ASSERT(file_finder2.IsFile());
            TS_ASSERT(!file_finder2.IsDir());

            // Check the path is as expected
            TS_ASSERT_EQUALS(file_finder2.GetAbsolutePath(), abs_path);

            // file_name is relative
            file_finder2.SetPath(file_name, RelativeTo::AbsoluteOrCwd);
            TS_ASSERT(file_finder2.Exists());
            TS_ASSERT(file_finder2.IsFile());
            TS_ASSERT(!file_finder2.IsDir());

            // Check the path is as expected
            TS_ASSERT_EQUALS(file_finder2.GetAbsolutePath(), abs_path);

            // Check we can extract the leaf name
            TS_ASSERT_EQUALS(file_finder.GetLeafName(), "FileFinder.hpp");
            TS_ASSERT_EQUALS(file_finder.GetLeafNameNoExtension(), "FileFinder");

            // And the parent folder name
            FileFinder parent("global/src", RelativeTo::ChasteSourceRoot);
            TS_ASSERT_EQUALS(file_finder.GetParent().GetAbsolutePath(), parent.GetAbsolutePath());
        }

        {
            // Now check a file in the output directory
            std::string dir_name = "TestFileFinder";
            OutputFileHandler handler(dir_name);
            std::string file_name = "TestFile";
            FileFinder file_finder(dir_name + "/" + file_name, RelativeTo::ChasteTestOutput);
            TS_ASSERT(!file_finder.Exists());
            TS_ASSERT(!file_finder.IsFile());
            TS_ASSERT(!file_finder.IsDir());

            // Check the path is as expected
            std::string abs_path = handler.GetOutputDirectoryFullPath() + file_name;
            TS_ASSERT_EQUALS(file_finder.GetAbsolutePath(), abs_path);

            // Check finding a sibling fails
            TS_ASSERT_THROWS_THIS(FileFinder("sibling", file_finder),
                                  "Reference path '" + abs_path + "' does not exist.");

            // Create the file
            out_stream fp = handler.OpenOutputFile(file_name);
            fp->close();
            TS_ASSERT(file_finder.Exists());

            // Finding a sibling is now ok
            FileFinder sibling("sibling", file_finder);
            TS_ASSERT_EQUALS(sibling.GetAbsolutePath(), handler.GetOutputDirectoryFullPath() + "sibling");

            // Check when providing an absolute path
            FileFinder file_finder2(abs_path, RelativeTo::Absolute);
            TS_ASSERT(file_finder2.Exists());
            TS_ASSERT_EQUALS(file_finder2.GetAbsolutePath(), abs_path);

            FileFinder file_finder3(abs_path, RelativeTo::AbsoluteOrCwd);
            TS_ASSERT(file_finder3.Exists());
            TS_ASSERT_EQUALS(file_finder3.GetAbsolutePath(), abs_path);
        }
    }

    void TestNewer()
    {
        FileFinder file("global/src/FileFinder.hpp", RelativeTo::ChasteSourceRoot);

        // A file can't be newer than itself
        TS_ASSERT(!file.IsNewerThan(file));

        // A newly created file better be newer than ourself!
        OutputFileHandler handler("TestFileFinder");
        out_stream fp = handler.OpenOutputFile("new_file");
        fp->close();
        FileFinder new_file("TestFileFinder/new_file", RelativeTo::ChasteTestOutput);
        TS_ASSERT(new_file.IsNewerThan(file));
        TS_ASSERT(!file.IsNewerThan(new_file));
    }

    void TestIsAbsolutePath()
    {
        TS_ASSERT(!FileFinder::IsAbsolutePath("global/src/FileFinder.hpp"));
        TS_ASSERT(FileFinder::IsAbsolutePath("/root"));
    }

    void TestDirFinder()
    {
        FileFinder dir("global", RelativeTo::ChasteSourceRoot);
        TS_ASSERT(dir.Exists());
        TS_ASSERT(dir.IsDir());
        TS_ASSERT(!dir.IsFile());
        std::string abs_path = std::string(ChasteBuildRootDir()) + "global/";
        TS_ASSERT_EQUALS(dir.GetAbsolutePath(), abs_path);

        FileFinder dir2("global", RelativeTo::CWD); // CWD should be the same as ChasteSourceRoot for tests
        TS_ASSERT(dir2.Exists());
        TS_ASSERT(dir2.IsDir());
        TS_ASSERT(!dir2.IsFile());
        TS_ASSERT_EQUALS(dir2.GetAbsolutePath(), dir.GetAbsolutePath());

        dir2.SetPath(dir.GetAbsolutePath(), RelativeTo::Absolute);
        TS_ASSERT(dir2.Exists());
        TS_ASSERT(dir2.IsDir());
        TS_ASSERT(!dir2.IsFile());
        TS_ASSERT_EQUALS(dir2.GetAbsolutePath(), dir.GetAbsolutePath());

        OutputFileHandler handler("TestFileFinder");
        FileFinder new_dir("TestFileFinder", RelativeTo::ChasteTestOutput);
        TS_ASSERT(new_dir.Exists());
        TS_ASSERT(new_dir.IsDir());
        TS_ASSERT(!new_dir.IsFile());
        TS_ASSERT_EQUALS(new_dir.GetAbsolutePath(), handler.GetOutputDirectoryFullPath());

        FileFinder missing_dir("TestFileFinder/SubDir", RelativeTo::ChasteTestOutput);
        TS_ASSERT(!missing_dir.Exists());
        TS_ASSERT(!missing_dir.IsDir());
        TS_ASSERT(!missing_dir.IsFile());
        TS_ASSERT_EQUALS(missing_dir.GetAbsolutePath(), handler.GetOutputDirectoryFullPath() + "SubDir");

        // Check the parent folder works for directories too
        TS_ASSERT_EQUALS(missing_dir.GetParent().GetAbsolutePath(), new_dir.GetAbsolutePath());

        // Can we give an empty local path?
        FileFinder source_root("", RelativeTo::ChasteSourceRoot);
        TS_ASSERT_EQUALS(dir.GetParent().GetAbsolutePath(), source_root.GetAbsolutePath());

        // Finding children also works
        FileFinder child("src", dir);
        TS_ASSERT_EQUALS(child.GetAbsolutePath(), abs_path + "src/");
    }

    void TestHandyFilenameOperations()
    {
        std::string path_with_spaces = "a path/with spaces.txt";

        FileFinder::ReplaceSpacesWithUnderscores(path_with_spaces);

        TS_ASSERT_EQUALS(path_with_spaces,"a_path/with_spaces.txt");

        FileFinder::ReplaceUnderscoresWithSpaces(path_with_spaces);

        TS_ASSERT_EQUALS(path_with_spaces,"a path/with spaces.txt");
    }

    void TestFaking()
    {
        FileFinder::FakePath(RelativeTo::ChasteSourceRoot, "test");
        FileFinder path("file", RelativeTo::ChasteSourceRoot);
        TS_ASSERT_EQUALS(path.GetAbsolutePath(), "test/file");

        FileFinder::FakePath(RelativeTo::CWD, "test1");
        FileFinder::FakePath(RelativeTo::ChasteSourceRoot, "test2");
        path.SetPath("file", RelativeTo::CWD);
        TS_ASSERT_EQUALS(path.GetAbsolutePath(), GetCurrentWorkingDirectory() + "/file");
        path.SetPath("file", RelativeTo::ChasteSourceRoot);
        TS_ASSERT_EQUALS(path.GetAbsolutePath(), "test2/file");

        FileFinder::StopFaking();
        path.SetPath("file", RelativeTo::ChasteSourceRoot);
        TS_ASSERT_EQUALS(path.GetAbsolutePath(), std::string(ChasteBuildRootDir()) + "file");
    }
};

#endif /*TESTFILEFINDER_HPP_*/
