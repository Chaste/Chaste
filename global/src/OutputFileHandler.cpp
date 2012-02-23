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

#include "OutputFileHandler.hpp"

#include <cstdlib>
#include <sys/stat.h>

#include "PetscTools.hpp"
#include "Exception.hpp"
#include "ArchiveLocationInfo.hpp"
#include "GetCurrentWorkingDirectory.hpp"
#include "FileFinder.hpp"

OutputFileHandler::OutputFileHandler(const std::string &rDirectory,
                                     bool cleanOutputDirectory)
{
    // Is it a valid request for a directory?
    if (rDirectory.find("..") != std::string::npos)
    {
        EXCEPTION("Will not create directory: " + rDirectory +
                  " due to it potentially being above, and cleaning, CHASTE_TEST_OUTPUT.");
    }

    mDirectory = MakeFoldersAndReturnFullPath(rDirectory);

    // Clean the directory (default)
    if (rDirectory != "" && cleanOutputDirectory) // Don't clean CHASTE_TEST_OUTPUT
    {
        std::string command = "test -e " + mDirectory + ".chaste_deletable_folder";
        int return_value = system(command.c_str());
        if (return_value != 0)
        {
            EXCEPTION("Cannot delete " + mDirectory + " because signature file \".chaste_deletable_folder\" is not present.");
        }

        // Are we the master process? Only the master should delete files
        if (PetscTools::AmMaster())
        {
            // Remove whatever was there before
            // Note that the /* part prevents removal of hidden files (.filename), which is useful in NFS systems
            ABORT_IF_NON0(system, "rm -rf " + mDirectory + "/*");
        }
        // Wait for master to finish before going on to use the directory.
        PetscTools::Barrier("OutputFileHandler");
    }
}

std::string OutputFileHandler::GetChasteTestOutputDirectory()
{
    char *chaste_test_output = getenv("CHASTE_TEST_OUTPUT");
    std::string directory_root;
    if (chaste_test_output == NULL || *chaste_test_output == 0)
    {
        // Default to 'testoutput' folder within the current directory
        directory_root = GetCurrentWorkingDirectory() + "/testoutput/";
    }
    else
    {
        if (*chaste_test_output != '/')
        {
            // It's a relative path; prepend with the CWD
            directory_root = GetCurrentWorkingDirectory() + "/" + chaste_test_output;
        }
        else
        {
#define COVERAGE_IGNORE
            // This branch is never taken on the build machines, because CHASTE_TEST_OUTPUT is set to a relative path.
            directory_root = std::string(chaste_test_output);
#undef COVERAGE_IGNORE
        }
    }
    AddTrailingSlash(directory_root);

    return directory_root;
}

std::string OutputFileHandler::MakeFoldersAndReturnFullPath(const std::string& rDirectory) const
{
    std::string directory_root = GetChasteTestOutputDirectory();
    std::string directories_to_add = rDirectory; // Get from a const to something we can mess with.
    AddTrailingSlash(directories_to_add);
    std::string directory = directory_root + directories_to_add;

    // Are we the master process? Only the master should make any new directories
    if (PetscTools::AmMaster())
    {
        // If necessary make the ChasteTestOutputDirectory - don't make it deleteable by Chaste
        std::string command = "test -d " + directory_root;
        int return_value = system(command.c_str());
        if (return_value!=0)
        {
            // We make as many folders as necessary here
            ABORT_IF_NON0(system, "mkdir -p " + directory_root);
        }

        // Now make all the sub-folders requested one-by-one and add the .chaste_deletable_folder file to them
        std::string remaining_directories = directories_to_add;
        std::string directory_to_add = "";

        // Create the output directory structure one folder at a time
        while (remaining_directories.find("/") != std::string::npos)
        {
            size_t found = remaining_directories.find_first_of("/");
            directory_to_add += remaining_directories.substr(0,found+1);
            remaining_directories = remaining_directories.substr(found+1);

            command = "test -d " + directory_root + directory_to_add;
            return_value = system(command.c_str());
            if (return_value!=0)
            {
                // We make only the next folder here
                ABORT_IF_NON0(system, "mkdir -p " + directory_root + directory_to_add);

                // Put the Chaste signature file into this folder
                ABORT_IF_NON0(system, "touch " + directory_root + directory_to_add + ".chaste_deletable_folder");
            }
        }
    }

    // Wait for master to finish before going on to use the directory.
    PetscTools::Barrier("OutputFileHandler::MakeFoldersAndReturnFullPath");

    return directory;
}

std::string OutputFileHandler::GetOutputDirectoryFullPath() const
{
    return mDirectory;
}

out_stream OutputFileHandler::OpenOutputFile(const std::string& rFileName,
                                             std::ios_base::openmode mode) const
{
    out_stream p_output_file(new std::ofstream((mDirectory+rFileName).c_str(), mode));
    if (!p_output_file->is_open())
    {
        EXCEPTION("Could not open file \"" + rFileName + "\" in " + mDirectory);
    }
    return p_output_file;
}

out_stream OutputFileHandler::OpenOutputFile(const std::string& rFileName,
                                             unsigned number,
                                             const std::string& rFileFormat,
                                             std::ios_base::openmode mode) const
{
    std::stringstream string_stream;
    string_stream << rFileName << number << rFileFormat;
    return OpenOutputFile(string_stream.str(), mode);
}

void OutputFileHandler::SetArchiveDirectory() const
{
    FileFinder dir(GetOutputDirectoryFullPath(), RelativeTo::Absolute);
    ArchiveLocationInfo::SetArchiveDirectory(dir);
}

void OutputFileHandler::AddTrailingSlash(std::string& rDirectory)
{
    // Add a trailing slash if not already there
    if (rDirectory!="" && !( *(rDirectory.end()-1) == '/'))
    {
        rDirectory = rDirectory + "/";
    }
}

FileFinder OutputFileHandler::CopyFileTo(const FileFinder& rSourceFile) const
{
    if (!rSourceFile.IsFile())
    {
        EXCEPTION("Can only copy single files:\n" << rSourceFile.GetAbsolutePath() << " is not a file.");
    }
    if (PetscTools::AmMaster())
    {
        ABORT_IF_NON0(system, "cp " + rSourceFile.GetAbsolutePath() + " " + GetOutputDirectoryFullPath());
    }
    PetscTools::Barrier("OutputFileHandler::CopyFileTo");
    return FindFile(rSourceFile.GetLeafName());
}

FileFinder OutputFileHandler::FindFile(std::string leafName) const
{
    return FileFinder(GetOutputDirectoryFullPath() + leafName, RelativeTo::Absolute);
}
