/*

Copyright (c) 2005-2026, University of Oxford.
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

#ifndef TESTFILESYSTEMPERMISSIONS_HPP_
#define TESTFILESYSTEMPERMISSIONS_HPP_

#include <cxxtest/TestSuite.h>
#include <fstream>

#include "FilesystemPermissions.hpp"
#include "OutputFileHandler.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "PetscTools.hpp"

class TestFilesystemPermissions : public CxxTest::TestSuite
{
private:
    void AssertPermissions(const fs::path& rPath, fs::perms expectedPermissions)
    {
#ifndef _MSC_VER
        TS_ASSERT_EQUALS(static_cast<unsigned>(fs::status(rPath).permissions() & fs::perms::all),
                         static_cast<unsigned>(expectedPermissions));
#endif
    }

public:
    void TestSettingFileAndDirectoryPermissions()
    {
        OutputFileHandler handler("TestFilesystemPermissions");
        fs::path directory_path(handler.GetOutputDirectoryFullPath());
        directory_path /= "direct_directory";
        fs::path file_path(handler.GetOutputDirectoryFullPath());
        file_path /= "direct_file.txt";

        if (PetscTools::AmMaster())
        {
            TS_ASSERT(FilesystemPermissions::CreateDirectoryWithPermissions(directory_path));
            AssertPermissions(directory_path, FilesystemPermissions::GetDirectoryPermissions());

            std::ofstream file(file_path);
            file << "content";
            file.close();
            fs::permissions(file_path, fs::perms::owner_read, fs::perm_options::replace);
            FilesystemPermissions::SetFilePermissions(file_path);
            AssertPermissions(file_path, FilesystemPermissions::GetFilePermissions());

            fs::permissions(directory_path, fs::perms::owner_all, fs::perm_options::replace);
            FilesystemPermissions::SetDirectoryPermissions(directory_path);
            AssertPermissions(directory_path, FilesystemPermissions::GetDirectoryPermissions());
        }
        PetscTools::Barrier("TestSettingFileAndDirectoryPermissions");
    }

    void TestCreateDirectoriesSetsPermissionsOnEachCreatedDirectory()
    {
        OutputFileHandler handler("TestFilesystemPermissionsNested");
        fs::path first_directory(handler.GetOutputDirectoryFullPath());
        first_directory /= "first";
        fs::path second_directory(first_directory);
        second_directory /= "second";
        fs::path third_directory(second_directory);
        third_directory /= "third";

        if (PetscTools::AmMaster())
        {
            TS_ASSERT(FilesystemPermissions::CreateDirectoriesWithPermissions(third_directory));
            AssertPermissions(first_directory, FilesystemPermissions::GetDirectoryPermissions());
            AssertPermissions(second_directory, FilesystemPermissions::GetDirectoryPermissions());
            AssertPermissions(third_directory, FilesystemPermissions::GetDirectoryPermissions());
        }
        PetscTools::Barrier("TestCreateDirectoriesSetsPermissionsOnEachCreatedDirectory");
    }

    void TestCopyFileSetsPermissionsOnDestination()
    {
        OutputFileHandler handler("TestFilesystemPermissionsCopy");
        fs::path source_path(handler.GetOutputDirectoryFullPath());
        source_path /= "source.txt";
        fs::path destination_path(handler.GetOutputDirectoryFullPath());
        destination_path /= "destination.txt";

        if (PetscTools::AmMaster())
        {
            std::ofstream source_file(source_path);
            source_file << "content";
            source_file.close();
            fs::permissions(source_path, fs::perms::owner_all, fs::perm_options::replace);

            FilesystemPermissions::CopyFileWithPermissions(source_path, destination_path);

            TS_ASSERT(fs::is_regular_file(destination_path));
            AssertPermissions(destination_path, FilesystemPermissions::GetFilePermissions());
        }
        PetscTools::Barrier("TestCopyFileSetsPermissionsOnDestination");
    }
};

#endif /*TESTFILESYSTEMPERMISSIONS_HPP_*/
