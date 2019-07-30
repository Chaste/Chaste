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

#ifndef TESTOUTPUTFILEHANDLER_HPP_
#define TESTOUTPUTFILEHANDLER_HPP_

#include <cxxtest/TestSuite.h>
#include <cstdlib>
#include <string>
#include <fstream>
#include <sstream>
#include <petsc.h>

#include "OutputFileHandler.hpp"
#include "BoostFilesystem.hpp"
#include "FileFinder.hpp"
#include "PetscTools.hpp"
#include "ChasteSyscalls.hpp"
#include "PetscSetupAndFinalize.hpp"

class TestOutputFileHandler : public CxxTest::TestSuite
{
public:

    void TestHandler()
    {
        // Test that CHASTE_TEST_OUTPUT always has a trailing slash even before
        // a class object is instantiated
        const std::string chaste_test_output(OutputFileHandler::GetChasteTestOutputDirectory());
        TS_ASSERT_EQUALS( *(chaste_test_output.end()-1), '/');

        // Make a handler that points straight to the CHASTE_TEST_OUTPUT directory.
        OutputFileHandler handler("");
        const std::string handler_path(handler.GetOutputDirectoryFullPath());
        TS_ASSERT(handler_path.length() > 0);
        TS_ASSERT_EQUALS(handler_path, handler.GetChasteTestOutputDirectory());
        TS_ASSERT_EQUALS(handler.GetRelativePath(), "");

        // Test that CHASTE_TEST_OUTPUT always has a trailing slash
        TS_ASSERT_EQUALS( *(handler_path.end()-1), '/');

        // Make a handler that points to a sub-directory.
        std::string dir = "testhandler";
        OutputFileHandler handler2(dir);
        std::string full_dir = handler2.GetOutputDirectoryFullPath();
        TS_ASSERT_EQUALS(full_dir.substr(full_dir.length()-dir.length()-1), dir+"/");
        TS_ASSERT_EQUALS(full_dir.substr(0, full_dir.length()-dir.length()-1), handler_path);
        TS_ASSERT_EQUALS(handler2.GetRelativePath(), dir);

        // We can also create handlers from a FileFinder (provided it points to a location in CHASTE_TEST_OUTPUT)
        OutputFileHandler handler3(handler.FindFile("testhandler2"));
        full_dir = handler3.GetOutputDirectoryFullPath();
        TS_ASSERT_EQUALS(full_dir.substr(full_dir.length()-dir.length()-2), dir+"2/");
        TS_ASSERT_EQUALS(full_dir.substr(0, full_dir.length()-dir.length()-2), handler_path);
        TS_ASSERT_EQUALS(handler3.GetRelativePath(), "testhandler2");

        // Check that all three handlers can create files
        out_stream p_file_stream;
        p_file_stream = handler.OpenOutputFile("test_file", std::ios::out);
        TS_ASSERT(FileFinder(handler_path + "test_file").Exists());

        p_file_stream = handler.OpenOutputFile("test_file2");
        TS_ASSERT(FileFinder(handler_path + "test_file2").Exists());

        p_file_stream = handler2.OpenOutputFile("test_file");
        TS_ASSERT(FileFinder(handler2.GetOutputDirectoryFullPath() + "test_file").Exists());

        p_file_stream = handler2.OpenOutputFile("test_", 34, ".txt");
        TS_ASSERT(FileFinder(handler2.GetOutputDirectoryFullPath() + "test_34.txt").Exists());

        p_file_stream = handler3.OpenOutputFile("test_file");
        TS_ASSERT(FileFinder(handler3.GetOutputDirectoryFullPath() + "test_file").Exists());

        // This should try to write files to /, which isn't allowed (we hope!)
        TS_ASSERT_THROWS_CONTAINS(OutputFileHandler bad_handler("../../../../../../../../../../../../../../../", false),
                                  "due to it potentially being above, and cleaning, CHASTE_TEST_OUTPUT.");
        TS_ASSERT_THROWS_CONTAINS(OutputFileHandler bad_handler("/", false),
                                  "The constructor argument to OutputFileHandler must be a relative path");
        TS_ASSERT_THROWS_CONTAINS(OutputFileHandler bad_handler(FileFinder("/"), false),
                                  "The location provided to OutputFileHandler must be inside CHASTE_TEST_OUTPUT");

        // Check the CopyFileTo method
        FileFinder source_file("global/test/TestOutputFileHandler.hpp", RelativeTo::ChasteSourceRoot);
        TS_ASSERT(!handler2.FindFile("TestOutputFileHandler.hpp").Exists());
        PetscTools::Barrier("TestOutputFileHandler-0");
        FileFinder dest_file = handler2.CopyFileTo(source_file);
        TS_ASSERT(dest_file.Exists());
        FileFinder missing_file("global/no_file", RelativeTo::ChasteSourceRoot);
        TS_ASSERT_THROWS_CONTAINS(handler2.CopyFileTo(missing_file), "Can only copy single files");
        FileFinder global_dir("global", RelativeTo::ChasteSourceRoot);
        TS_ASSERT_THROWS_CONTAINS(handler2.CopyFileTo(global_dir), "Can only copy single files");

        // We don't want other people using CHASTE_TEST_OUTPUT whilst we are messing with it!
        PetscTools::Barrier("TestOutputFileHandler-1");

        // Test that the environment variable actually influences the location of files
        {

            setenv("CHASTE_TEST_OUTPUT", "", 1/*Overwrite*/);
            // Predict where Chaste puts output when CHASTE_TEST_OUTPUT is not set
            std::stringstream  tmp_directory;
            if (getenv("USER")!=NULL)
             {
                 tmp_directory << "/tmp/" << getenv("USER") << "/testoutput/NoEnvironmentForTestoutput";
             }
             else
             {
                 // No $USER in environment (which may be the case in Docker)
                 tmp_directory << "/tmp/chaste/testoutput/NoEnvironmentForTestoutput";
             }
            // Check this folder is not present
            FileFinder test_folder(tmp_directory.str(), RelativeTo::Absolute);
            TS_ASSERT(!test_folder.Exists());

            PetscTools::Barrier("TestOutputFileHandler-2");

            // Make a folder and erase it - NB only master can erase files and check it is successful!
            OutputFileHandler handler4("NoEnvironmentForTestoutput");

            TS_ASSERT(test_folder.Exists());
            PetscTools::Barrier("TestOutputFileHandler-2b");
            if (PetscTools::AmMaster())
            {
                test_folder.Remove();
            }
            PetscTools::Barrier("TestOutputFileHandler-2c");
        }

        {
            setenv("CHASTE_TEST_OUTPUT", "config__cyborg__T800__cooper", 1/*Overwrite*/);
            // Test that CHASTE_TEST_OUTPUT always has a trailing slash even before
            // a class object is instantiated and when the directory does not exist
            const std::string nonexistent_test_path(OutputFileHandler::GetChasteTestOutputDirectory());
            TS_ASSERT_EQUALS( *(nonexistent_test_path.end()-1), '/');
        }

        {
            // Check this folder is not present
            std::string test_folder("somewhere_without_trailing_forward_slash");
            TS_ASSERT(!FileFinder(test_folder, RelativeTo::CWD).Exists());
            PetscTools::Barrier("TestOutputFileHandler-3");

            setenv("CHASTE_TEST_OUTPUT", test_folder.c_str(), 1/*Overwrite*/);

            // Make a folder using a FileFinder, for coverage of the case where the root output folder doesn't exist
            FileFinder sub_folder("test_folder", RelativeTo::ChasteTestOutput);
            TS_ASSERT(!sub_folder.Exists());
            PetscTools::Barrier("TestOutputFileHandler-3a");
            OutputFileHandler creating_handler(sub_folder);
            TS_ASSERT(sub_folder.Exists());

            // Make a folder
            OutputFileHandler handler5("whatever");
            TS_ASSERT(FileFinder(test_folder, RelativeTo::CWD).Exists());
            PetscTools::Barrier("TestOutputFileHandler-3b");

            // Erase it
            if (PetscTools::AmMaster())
            {
                FileFinder(test_folder).DangerousRemove();
            }
        }

        // Reset the location of CHASTE_TEST_OUTPUT
        setenv("CHASTE_TEST_OUTPUT", chaste_test_output.c_str(), 1/*Overwrite*/);

        // We don't want other people using CHASTE_TEST_OUTPUT while we are messing with it!
        PetscTools::Barrier("TestOutputFileHandler-4");

        // Coverage of the case where we can't open a file for writing
        OutputFileHandler handler6("no_write_access");
        if (PetscTools::AmMaster())
        {
            std::string dir_path =  handler6.GetOutputDirectoryFullPath();
#ifndef _MSC_VER
            // This test can never pass on modern Windows OS! See: http://support.microsoft.com/kb/326549
            // You can't change DIRECTORY attributes
            chmod(dir_path.c_str(), CHASTE_READONLY);
            TS_ASSERT_THROWS_CONTAINS(p_file_stream = handler6.OpenOutputFile("test_file"),
                                      "Could not open file");
#endif
            chmod(dir_path.c_str(), CHASTE_READ_WRITE_EXECUTE);
            fs::remove(dir_path + ".chaste_deletable_folder");
            fs::remove(dir_path);
        }

        // Check behaviour of FindFile("")
        FileFinder handler_self = handler.FindFile("");
        TS_ASSERT_EQUALS(handler_self.GetAbsolutePath(), handler.GetOutputDirectoryFullPath());
    }

    void TestWeCanOnlyDeleteFoldersWeHaveMadeOurselves()
    {
        std::string test_folder = "cannot_delete_me";
        if (PetscTools::AmMaster())
        {
            ABORT_IF_THROWS(fs::create_directories(OutputFileHandler::GetChasteTestOutputDirectory() + test_folder));
        }
        // Wait until directory has been created, and check it exists
        PetscTools::Barrier("TestWeCanOnlyDeleteFoldersWeHaveMadeOurselves-1");
        FileFinder cannot_delete(test_folder, RelativeTo::ChasteTestOutput);
        TS_ASSERT(cannot_delete.IsDir());

        // Try to use it as an output folder
        TS_ASSERT_THROWS_CONTAINS(OutputFileHandler bad_handler(test_folder),
                                  "because signature file \".chaste_deletable_folder\" is not present");

        // Tidy up
        if (PetscTools::AmMaster())
        {
            TS_ASSERT(cannot_delete.Exists());
            cannot_delete.DangerousRemove();
            TS_ASSERT(!cannot_delete.Exists());
        }

        // Now create a folder the proper way
        test_folder = "can_delete_me";
        OutputFileHandler handler(test_folder);
        out_stream p_file_stream = handler.OpenOutputFile("test_file");
        p_file_stream->close(); // Windows does not like deleting open files

        // Test file is present
        FileFinder test_file = handler.FindFile("test_file");
        TS_ASSERT(test_file.Exists());
        PetscTools::Barrier("TestWeCanOnlyDeleteFoldersWeHaveMadeOurselves-2");

        OutputFileHandler handler2(test_folder, false /* don't clean */);

        // Test file is still present
        TS_ASSERT(test_file.Exists());
        PetscTools::Barrier("TestWeCanOnlyDeleteFoldersWeHaveMadeOurselves-3");

        OutputFileHandler handler3(test_folder, true /* do clean */);

        // Test file is deleted
        TS_ASSERT(!test_file.Exists());
        PetscTools::Barrier("TestWeCanOnlyDeleteFoldersWeHaveMadeOurselves-4");

        // Check we can delete the test_folder too
        if (PetscTools::AmMaster())
        {
            FileFinder folder = handler.FindFile("");
            TS_ASSERT(folder.Exists());
            folder.Remove();
            TS_ASSERT(!folder.Exists());
        }

        // Test we can make a directory of folders and delete them all
        OutputFileHandler handler4("what_about_me/and_me/and_me/and_da_da_da", true);

        // Check we have made a subdirectory
        FileFinder sub_folder("what_about_me/and_me", RelativeTo::ChasteTestOutput);
        TS_ASSERT(sub_folder.IsDir());
        PetscTools::Barrier("TestWeCanOnlyDeleteFoldersWeHaveMadeOurselves-5");

        OutputFileHandler handler5("what_about_me", true);

        // Check we have wiped the sub-directories
        TS_ASSERT(!sub_folder.Exists());
        PetscTools::Barrier("TestWeCanOnlyDeleteFoldersWeHaveMadeOurselves-6");

        // Check we can delete the main directory too
        if (PetscTools::AmMaster())
        {
            FileFinder folder = handler5.FindFile("");
            TS_ASSERT(folder.Exists());
            folder.Remove();
            TS_ASSERT(!folder.Exists());
        }
    }
};

#endif /*TESTOUTPUTFILEHANDLER_HPP_*/
