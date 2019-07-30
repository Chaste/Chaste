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

#ifndef TESTFILEFINDER_HPP_
#define TESTFILEFINDER_HPP_

#include <cxxtest/TestSuite.h>
#include "FileFinder.hpp"
#include "BoostFilesystem.hpp"
#include "ChasteBuildRoot.hpp"
#include "OutputFileHandler.hpp"
#include "GetCurrentWorkingDirectory.hpp"
#include "Warnings.hpp"
//This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

class TestFileFinder : public CxxTest::TestSuite
{
public:

    void TestFileFinderOpening()
    {
        {
            // Can we find our own source file?
            std::string file_name = "global/src/FileFinder.hpp";
            FileFinder file_finder(file_name, RelativeTo::ChasteSourceRoot);
            TS_ASSERT(file_finder.Exists());
            TS_ASSERT(file_finder.IsFile());
            TS_ASSERT(!file_finder.IsDir());
            TS_ASSERT(!file_finder.IsEmpty());
            TS_ASSERT(file_finder.IsPathSet());

            // Check the path is as expected
            std::string abs_path = ChasteSourceRootDir() + file_name;
            TS_ASSERT_EQUALS(file_finder.GetAbsolutePath(), abs_path);

            // CWD should be the Chaste source root
            FileFinder file_finder2(file_name, RelativeTo::CWD);
            std::cout << file_finder2.GetAbsolutePath();
            TS_ASSERT(file_finder2.Exists());
            TS_ASSERT(file_finder2.IsFile());
            TS_ASSERT(!file_finder2.IsDir());

            // Check the path is as expected
            TS_ASSERT( fs::equivalent( fs::path(file_finder2.GetAbsolutePath()) , fs::path(abs_path) ) );

            // file_name is relative
            file_finder2.SetPath(file_name, RelativeTo::AbsoluteOrCwd);
            TS_ASSERT(file_finder2.Exists());
            TS_ASSERT(file_finder2.IsFile());
            TS_ASSERT(!file_finder2.IsDir());

            // Check the path is as expected
            TS_ASSERT( fs::equivalent( fs::path(file_finder2.GetAbsolutePath()), fs::path(abs_path) ) );

            // Check we can extract the leaf name
            TS_ASSERT_EQUALS(file_finder.GetLeafName(), "FileFinder.hpp");
            TS_ASSERT_EQUALS(file_finder.GetLeafNameNoExtension(), "FileFinder");
            TS_ASSERT_EQUALS(file_finder.GetExtension(), ".hpp");

            // And the parent folder name
            FileFinder parent("global/src", RelativeTo::ChasteSourceRoot);
            TS_ASSERT_EQUALS(file_finder.GetParent().GetAbsolutePath(), parent.GetAbsolutePath());
            TS_ASSERT_EQUALS(parent.GetLeafName(), "src");
            TS_ASSERT_EQUALS(parent.GetLeafNameNoExtension(), "src");
            TS_ASSERT_EQUALS(parent.GetExtension(), "");

            // Check we can construct from a Boost path or a string
            TS_ASSERT(fs::equivalent(fs::path(FileFinder(fs::path(file_name)).GetAbsolutePath()),fs::path(abs_path)));
            TS_ASSERT(fs::equivalent(fs::path(FileFinder(file_name).GetAbsolutePath()),fs::path(abs_path)));
            TS_ASSERT(fs::equivalent(fs::path(FileFinder(abs_path).GetAbsolutePath()),fs::path(abs_path)));
            TS_ASSERT(fs::equivalent(fs::path(FileFinder("global/src/FileFinder.hpp").GetAbsolutePath()),fs::path(abs_path)));
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

            // Check finding a sibling and testing for emptiness fail
            TS_ASSERT_THROWS_THIS(file_finder.IsEmpty(), "The path '" + abs_path + "' does not exist.");
            TS_ASSERT_THROWS_THIS(FileFinder("sibling", file_finder),
                                  "Reference path '" + abs_path + "' does not exist.");

            // Create the file
            out_stream fp = handler.OpenOutputFile(file_name);
            fp->close();
            TS_ASSERT(file_finder.Exists());
            TS_ASSERT(file_finder.IsEmpty());

            // Finding a sibling is now ok
            FileFinder sibling("sibling", file_finder);
            TS_ASSERT_EQUALS(sibling.GetAbsolutePath(), handler.GetOutputDirectoryFullPath() + "sibling");
            TS_ASSERT(sibling.IsPathSet());

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
#ifndef _MSC_VER
        //This would always fail on Windows: it's not a single-root OS
        TS_ASSERT(FileFinder::IsAbsolutePath("/root"));
#endif
    }

    void TestDirFinder()
    {
        FileFinder dir("global", RelativeTo::ChasteSourceRoot);
        TS_ASSERT(dir.Exists());
        TS_ASSERT(dir.IsDir());
        TS_ASSERT(!dir.IsFile());
        TS_ASSERT(!dir.IsEmpty());
        TS_ASSERT(dir.IsPathSet());
        std::string abs_path = std::string(ChasteSourceRootDir()) + "global/";
        TS_ASSERT(fs::equivalent(fs::path(dir.GetAbsolutePath()),fs::path(abs_path)));

        FileFinder dir2("global", RelativeTo::CWD); // CWD should be the same as ChasteSourceRoot for tests
        TS_ASSERT(dir2.Exists());
        TS_ASSERT(dir2.IsDir());
        TS_ASSERT(!dir2.IsFile());
        TS_ASSERT(fs::equivalent(fs::path(dir2.GetAbsolutePath()),fs::path(dir.GetAbsolutePath())));

        dir2.SetPath(dir.GetAbsolutePath(), RelativeTo::Absolute);
        TS_ASSERT(dir2.Exists());
        TS_ASSERT(dir2.IsDir());
        TS_ASSERT(!dir2.IsFile());
        TS_ASSERT_EQUALS(dir2.GetAbsolutePath(), dir.GetAbsolutePath());

        OutputFileHandler handler("TestFileFinder");
        FileFinder new_dir("TestFileFinder", RelativeTo::ChasteTestOutput);
        TS_ASSERT(new_dir.Exists());
        TS_ASSERT(new_dir.IsDir());
        TS_ASSERT(new_dir.IsEmpty());
        TS_ASSERT(!new_dir.IsFile());
        TS_ASSERT_EQUALS(new_dir.GetAbsolutePath(), handler.GetOutputDirectoryFullPath());
        // Path for an existing dir will end in a '/'
        TS_ASSERT_EQUALS(*(new_dir.GetAbsolutePath().rbegin()), '/');

        FileFinder missing_dir("TestFileFinder/SubDir", RelativeTo::ChasteTestOutput);
        TS_ASSERT(!missing_dir.Exists());
        TS_ASSERT(!missing_dir.IsDir());
        TS_ASSERT(!missing_dir.IsFile());
        // Note no slash on the end of a missing dir's path
        TS_ASSERT_EQUALS(missing_dir.GetAbsolutePath(), handler.GetOutputDirectoryFullPath() + "SubDir");
        // But once we create the dir it will have one, even with the same finder
        OutputFileHandler sub_handler(missing_dir);
        TS_ASSERT_EQUALS(missing_dir.GetAbsolutePath(), handler.GetOutputDirectoryFullPath() + "SubDir/");

        // Check the parent folder works for directories too
        TS_ASSERT_EQUALS(missing_dir.GetParent().GetAbsolutePath(), new_dir.GetAbsolutePath());

        // Can we give an empty local path?
        FileFinder source_root("", RelativeTo::ChasteSourceRoot);
        TS_ASSERT_EQUALS(dir.GetParent().GetAbsolutePath(), source_root.GetAbsolutePath());

        // Finding children also works
        FileFinder child("src", dir);
        TS_ASSERT_EQUALS(child.GetAbsolutePath(), abs_path + "src/");

        // We can also compute relative paths
        TS_ASSERT_EQUALS(child.GetRelativePath(dir), "src/");
        TS_ASSERT_EQUALS(FileFinder("SConscript", dir).GetRelativePath(dir), "SConscript");
        TS_ASSERT_THROWS_CONTAINS(child.GetRelativePath(new_dir), "' is not relative to '");
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
        TS_ASSERT_EQUALS(path.GetAbsolutePath(), std::string(ChasteSourceRootDir()) + "file");
    }

    void TestRemove()
    {
        // We shouldn't be able to remove unsafe files, if possible
        FileFinder bad1a("global/src", RelativeTo::ChasteSourceRoot);
        FileFinder bad1b("/", RelativeTo::Absolute);
        TS_ASSERT_THROWS_CONTAINS(bad1a.Remove(), "is not located within the Chaste test output folder");
        TS_ASSERT_THROWS_CONTAINS(bad1b.Remove(), "is not located within the Chaste test output folder");
        FileFinder bad2("../..", RelativeTo::ChasteTestOutput);
        TS_ASSERT_THROWS_CONTAINS(bad2.Remove(), "contains a dangerous path component.");

        // We can delete individual files
        OutputFileHandler handler("TestFileFinder/TestRemove");
        handler.OpenOutputFile("delete_me");
        FileFinder file = handler.FindFile("delete_me");
        TS_ASSERT(file.Exists());
        file.Remove();
        TS_ASSERT(!file.Exists());

        // We can recursively delete folders
        FileFinder dir("TestFileFinder", RelativeTo::ChasteTestOutput);
        FileFinder subdir("TestFileFinder/TestRemove", RelativeTo::ChasteTestOutput);
        TS_ASSERT(subdir.Exists());
        dir.Remove();
        TS_ASSERT(!subdir.Exists());
        TS_ASSERT(!dir.Exists());

        // We can only delete (content in) folders created by an OutputFileHandler unless we force it
        fs::create_directory(dir.GetAbsolutePath());
        TS_ASSERT(dir.IsDir());
        TS_ASSERT_THROWS_CONTAINS(dir.Remove(),
                                  "because the signature file '.chaste_deletable_folder' is not present.");
        TS_ASSERT(dir.IsDir());
        file.SetPath("file", dir);
        { std::ofstream(file.GetAbsolutePath().c_str()); }
        TS_ASSERT(file.IsFile());
        TS_ASSERT_THROWS_CONTAINS(file.Remove(),
                                  "because the signature file '.chaste_deletable_folder' is not present.");
        TS_ASSERT(file.IsFile());
        file.DangerousRemove();
        TS_ASSERT(!file.Exists());
        dir.DangerousRemove();
        TS_ASSERT(!dir.Exists());

        FileFinder obscure_file("/TodayIsThe25thOfOctober2012AndWeLikedRafsCake.obscure", RelativeTo::Absolute);
        TS_ASSERT_THROWS_CONTAINS(obscure_file.Remove(),
                "as it is not located within the Chaste test output folder");
        TS_ASSERT_THROWS_CONTAINS(obscure_file.DangerousRemove(), ", the Chaste source folder");
    }

    void TestFindMatches()
    {
        std::string dirname("TestFileFinder_TestFindMatches");
        OutputFileHandler handler(dirname);
        FileFinder dir(dirname, RelativeTo::ChasteTestOutput);

        // Create some files to find
        const unsigned N = 5;
        for (unsigned i=0; i<N; ++i)
        {
            handler.OpenOutputFile("file", i, ".txt");
        }

        // ? matches a single character
        std::vector<FileFinder> matches = dir.FindMatches("file?.txt");
        TS_ASSERT_EQUALS(matches.size(), N);
        // Trailing * matches anything
        matches = dir.FindMatches("file*");
        TS_ASSERT_EQUALS(matches.size(), N);
        // Initial * matches anything
        matches = dir.FindMatches("*txt");
        TS_ASSERT_EQUALS(matches.size(), N);
        // Hidden files are ignored
        matches = dir.FindMatches("*");
        TS_ASSERT_EQUALS(matches.size(), N);
        // We can combine * & ? (to some extent)
        matches = dir.FindMatches("*le?.t??");
        TS_ASSERT_EQUALS(matches.size(), N);
        matches = dir.FindMatches("????1.*");
        TS_ASSERT_EQUALS(matches.size(), 1u);
        // Empty pattern matches nothing
        matches = dir.FindMatches("");
        TS_ASSERT_EQUALS(matches.size(), 0u);

        // Only initial or trailing * are supported
        matches = dir.FindMatches("file*txt");
        TS_ASSERT_EQUALS(matches.size(), 0u);
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNumWarnings(), 1u);
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNextWarningMessage(),
                         "A '*' only has special meaning at the start or end of a pattern.");
        Warnings::QuietDestroy();
        // This sort of pattern is hard too
        TS_ASSERT_THROWS_THIS(dir.FindMatches("*le?.t*"),
                              "The '*' wildcard may not be used at both the start and end of the pattern if the '?' wildcard is also used.");

        // Only works on folders
        FileFinder file("file0.txt", dir);
        TS_ASSERT_THROWS_CONTAINS(file.FindMatches("*"), "as it is not a directory.");
    }

    void TestCopying()
    {
        FileFinder source("global/test/TestFileFinder.hpp", RelativeTo::ChasteSourceRoot);
        std::string dest_dir_name("TestFileFinder_TestCopying");
        OutputFileHandler handler(dest_dir_name);

        // Copy to folder
        FileFinder dest_dir(dest_dir_name, RelativeTo::ChasteTestOutput);
        FileFinder expected_dest("TestFileFinder.hpp", dest_dir);
        TS_ASSERT(!expected_dest.Exists());
        FileFinder dest = source.CopyTo(dest_dir);
        TS_ASSERT(dest.IsFile());
        TS_ASSERT(fs::equivalent(fs::path(dest.GetAbsolutePath()),fs::path(expected_dest.GetAbsolutePath())));

        // Copy to existing file
        dest = source.CopyTo(dest_dir);
        TS_ASSERT(dest.IsFile());
        TS_ASSERT(fs::equivalent(fs::path(dest.GetAbsolutePath()),fs::path(expected_dest.GetAbsolutePath())));

        // Copy to new name
        expected_dest.SetPath("new_name", dest_dir);
        TS_ASSERT(!expected_dest.Exists());
        dest = source.CopyTo(expected_dest);
        TS_ASSERT(dest.IsFile());
        TS_ASSERT_EQUALS(dest.GetAbsolutePath(), expected_dest.GetAbsolutePath());

        // Recursive copy of a folder
        // Firstly we create some sub-folders to copy
        OutputFileHandler sub_handler1(dest_dir_name + "/sub1");
        OutputFileHandler sub_handler2(dest_dir_name + "/sub1/sub2");
        sub_handler2.OpenOutputFile("test.txt")->close();
        // Then we copy
        source = dest_dir;
        OutputFileHandler copy_handler("TestFileFinder_TestCopyingFolders");
        dest_dir = copy_handler.FindFile("");
        TS_ASSERT(dest_dir.IsDir());
        TS_ASSERT(!copy_handler.FindFile(dest_dir_name).Exists()); // We're copying into an empty folder
        dest = source.CopyTo(dest_dir);
        TS_ASSERT(copy_handler.FindFile(dest_dir_name).IsDir());
        TS_ASSERT(copy_handler.FindFile(dest_dir_name + "/sub1").IsDir());
        TS_ASSERT(copy_handler.FindFile(dest_dir_name + "/sub1/sub2").IsDir());
        TS_ASSERT(copy_handler.FindFile(dest_dir_name + "/sub1/sub2/test.txt").IsFile());
        TS_ASSERT(copy_handler.FindFile(dest_dir_name + "/TestFileFinder.hpp").IsFile());

        // Copy the tree to a name that doesn't exist
        dest_dir = copy_handler.FindFile("second_copy");
        TS_ASSERT(!dest_dir.Exists());
        dest = source.CopyTo(dest_dir);
        TS_ASSERT(dest.IsDir());
        TS_ASSERT_EQUALS(dest.GetLeafName(), "second_copy");
        TS_ASSERT(FileFinder("TestFileFinder.hpp", dest).IsFile());
        TS_ASSERT(FileFinder("sub1", dest).IsDir());

        // Error cases
        TS_ASSERT_THROWS_CONTAINS(FileFinder("global/no_file", RelativeTo::ChasteSourceRoot).CopyTo(dest),
                                  "as it does not exist.");
        // Trying to overwrite a file with a folder
        dest = handler.FindFile("new_name");
        TS_ASSERT(dest.IsFile());
        TS_ASSERT_THROWS_CONTAINS(source.CopyTo(dest), "as it would overwrite an existing file.");
        // Trying to copy a folder that has already been copied
        dest_dir = copy_handler.FindFile("second_copy");
        dest_dir = copy_handler.FindFile("");
        TS_ASSERT(dest_dir.IsDir());
        TS_ASSERT(FileFinder(source.GetLeafName(), dest_dir).IsDir());
        TS_ASSERT_THROWS_CONTAINS(source.CopyTo(dest_dir), "as it would overwrite an existing file.");
    }

    void TestDefaultConstructor()
    {
        FileFinder unset;
        TS_ASSERT(!unset.IsPathSet());
        unset.SetPath("", RelativeTo::ChasteSourceRoot);
        TS_ASSERT(unset.IsPathSet());
    }
};

#endif /*TESTFILEFINDER_HPP_*/
