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

#ifndef TESTOUTPUTFILEHANDLER_HPP_
#define TESTOUTPUTFILEHANDLER_HPP_

#include <cxxtest/TestSuite.h>
#include <string>
#include <fstream>
#include <sstream>
#include <unistd.h> //For rmdir()
#include <petsc.h>
#include <sys/stat.h> //For chmod()
#include "OutputFileHandler.hpp"
#include "FileFinder.hpp"
#include "PetscTools.hpp"
#include "PetscSetupAndFinalize.hpp"

class TestOutputFileHandler : public CxxTest::TestSuite
{
public:

    void TestHandler() throw(Exception)
    {
        // Make a handler that points straight to the CHASTE_TEST_OUTPUT directory.
        OutputFileHandler handler("");
        TS_ASSERT(handler.GetOutputDirectoryFullPath().length() > 0);
        TS_ASSERT_EQUALS(handler.GetOutputDirectoryFullPath(),handler.GetChasteTestOutputDirectory());

        // Make a handler that points straight to a sub-directory.
        std::string dir = "testhandler";
        OutputFileHandler handler2(dir);
        std::string full_dir = handler2.GetOutputDirectoryFullPath();
        TS_ASSERT_EQUALS(full_dir.substr(full_dir.length()-dir.length()-1), dir+"/");
        TS_ASSERT_EQUALS(full_dir, handler2.GetOutputDirectoryFullPath());

        // Check that both can create files
        out_stream p_file_stream;
        p_file_stream = handler.OpenOutputFile("test_file", std::ios::out);
        EXPECT0(system, "test -e " + handler.GetOutputDirectoryFullPath() + "test_file");

        p_file_stream = handler.OpenOutputFile("test_file");
        EXPECT0(system, "test -e " + handler.GetOutputDirectoryFullPath() + "test_file");

        p_file_stream = handler2.OpenOutputFile("test_file");
        EXPECT0(system, "test -e " + handler2.GetOutputDirectoryFullPath() + "test_file");

        p_file_stream = handler2.OpenOutputFile("test_",34,".txt");
        EXPECT0(system, "test -e " + handler2.GetOutputDirectoryFullPath() + "test_34.txt");

        // This should try to write files to /, which isn't allowed (we hope!)
        TS_ASSERT_THROWS_CONTAINS(OutputFileHandler handler3("../../../../../../../../../../../../../../../",false),
                "due to it potentially being above, and cleaning, CHASTE_TEST_OUTPUT.");

        // Check the CopyFileTo method
        FileFinder source_file("global/test/TestOutputFileHandler.hpp", RelativeTo::ChasteSourceRoot);
        FileFinder dest_file = handler2.CopyFileTo(source_file);
        TS_ASSERT(dest_file.Exists());
        FileFinder missing_file("global/no_file", RelativeTo::ChasteSourceRoot);
        TS_ASSERT_THROWS_CONTAINS(handler2.CopyFileTo(missing_file), "Can only copy single files");
        FileFinder global_dir("global", RelativeTo::ChasteSourceRoot);
        TS_ASSERT_THROWS_CONTAINS(handler2.CopyFileTo(global_dir), "Can only copy single files");

        // We don't want other people using CHASTE_TEST_OUTPUT whilst we are messing with it!
        PetscTools::Barrier();

        // Test that the Chaste directory actually influences the location of files
        char *chaste_test_output = getenv("CHASTE_TEST_OUTPUT");

        setenv("CHASTE_TEST_OUTPUT", "", 1/*Overwrite*/);

        // Check this folder is not present
        std::string command = "test -d testoutput/whatever";
        int return_value = system(command.c_str());
        TS_ASSERT_DIFFERS(return_value, 0);
        PetscTools::Barrier();

        // Make a folder and erase it - NB only master can erase files and check it is successful!
        OutputFileHandler handler4("whatever");
        if (PetscTools::AmMaster())
        {
            ABORT_IF_NON0(system, "rm -rf testoutput/whatever");
        }

        // Check this folder is not present
        command = "test -d somewhere_without_trailing_forward_slash";
        return_value = system(command.c_str());
        TS_ASSERT_DIFFERS(return_value, 0);
        PetscTools::Barrier();

        setenv("CHASTE_TEST_OUTPUT", "somewhere_without_trailing_forward_slash", 1/*Overwrite*/);

        // Make a folder
        OutputFileHandler handler5("whatever");

        // Erase it
        if (PetscTools::AmMaster())
        {
            ABORT_IF_NON0(system, "rm -rf somewhere_without_trailing_forward_slash");
        }

        // Reset the location of CHASTE_TEST_OUTPUT
        setenv("CHASTE_TEST_OUTPUT", chaste_test_output, 1/*Overwrite*/);

        // We don't want other people using CHASTE_TEST_OUTPUT whilst we are messing with it!
        PetscTools::Barrier();

        // Coverage of the case where we can't open a file for writing
        OutputFileHandler handler6("no_write_access");
        if (PetscTools::AmMaster())
        {
            command = handler6.GetOutputDirectoryFullPath();
            chmod(command.c_str(),0444);
            TS_ASSERT_THROWS_CONTAINS(p_file_stream = handler6.OpenOutputFile("test_file"),
                    "Could not open file");
            chmod(command.c_str(),0755);
            rmdir(command.c_str());
        }
    }

    void TestWeCanOnlyDeleteFoldersWeHaveMadeOurselves() throw(Exception)
    {
        std::string command;
        if (PetscTools::AmMaster())
        {
            command = "mkdir -p " + OutputFileHandler::GetChasteTestOutputDirectory() + "cannot_delete_me";
            system(command.c_str());
        }

        // Wait until directory has been created
        PetscTools::Barrier();

        command = "test -d " + OutputFileHandler::GetChasteTestOutputDirectory() + "cannot_delete_me";

        // Check this folder has been created...
        TS_ASSERT_EQUALS(system(command.c_str()), 0);

        // Try to use it as an output folder
        TS_ASSERT_THROWS_CONTAINS(OutputFileHandler bad_handler("cannot_delete_me"),
                                  "because signature file \".chaste_deletable_folder\" is not present");

        OutputFileHandler handler("can_delete_me");
        out_stream p_file_stream;
        p_file_stream = handler.OpenOutputFile("test_file");

        // Test file is still present
        command = "test -e " + handler.GetOutputDirectoryFullPath() + "test_file";
        int return_value = system(command.c_str());
        TS_ASSERT_EQUALS(return_value, 0);
        PetscTools::Barrier();

        OutputFileHandler handler2("can_delete_me", false);

        // Test file is still present
        return_value = system(command.c_str());
        TS_ASSERT_EQUALS(return_value, 0);
        PetscTools::Barrier();

        OutputFileHandler handler3("can_delete_me", true);

        // Test file is deleted
        return_value = system(command.c_str());
        TS_ASSERT_DIFFERS(return_value, 0);
        PetscTools::Barrier();

        // Test we can make a directory of folders and delete them all
        OutputFileHandler handler4("what_about_me/and_me/and_me/and_da_da_da",true);

        // Check we have made a subdirectory
        command = "test -d " + OutputFileHandler::GetChasteTestOutputDirectory() + "what_about_me/and_me";
        return_value = system(command.c_str());
        TS_ASSERT_EQUALS(return_value, 0);
        PetscTools::Barrier();

        OutputFileHandler handler5("what_about_me", true);

        // Check we have wiped the sub-directories
        return_value = system(command.c_str());
        TS_ASSERT_DIFFERS(return_value, 0);
        PetscTools::Barrier();
    }
};

#endif /*TESTOUTPUTFILEHANDLER_HPP_*/
