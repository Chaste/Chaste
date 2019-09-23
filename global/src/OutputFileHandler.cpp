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

#include "OutputFileHandler.hpp"

#include <cstdlib>
#include <fstream>
#include <sstream>

#include "ArchiveLocationInfo.hpp"
#include "BoostFilesystem.hpp"
#include "Exception.hpp"
#include "FileFinder.hpp"
#include "GetCurrentWorkingDirectory.hpp"
#include "PetscTools.hpp"


const std::string OutputFileHandler::SIG_FILE_NAME(".chaste_deletable_folder");

/**
 * Recursively remove the contents of the given folder, but leave any hidden
 * files present at the top level.
 *
 * @param rPath  path to the directory to clean
 * @param isTop  whether this is the top level directory
 */
void CleanFolder(const fs::path& rPath, bool isTop=true);

void CleanFolder(const fs::path& rPath, bool isTop)
{
    assert(fs::is_directory(rPath));
    fs::directory_iterator end_iter;
    // First recursively remove the children
    for (fs::directory_iterator dir_iter(rPath); dir_iter != end_iter; ++dir_iter)
    {
        if (fs::is_directory(dir_iter->status()))
        {
            CleanFolder(dir_iter->path(), false);
        }
        else
        {
            const fs::path& r_item_path(dir_iter->path());
            if (!isTop || PATH_LEAF_NAME(r_item_path)[0] != '.')
            {
                fs::remove(r_item_path);
            }
        }
    }
    // Now remove the folder itself, if not top
    if (!isTop)
    {
        fs::remove(rPath);
    }
}


OutputFileHandler::OutputFileHandler(const std::string& rDirectory,
                                     bool cleanOutputDirectory)
{
    CommonConstructor(rDirectory, cleanOutputDirectory);
}


OutputFileHandler::OutputFileHandler(const FileFinder& rDirectory,
                                     bool cleanOutputDirectory)
{
    FileFinder output_root("", RelativeTo::ChasteTestOutput);
    std::string relative_path;
    try
    {
        relative_path = rDirectory.GetRelativePath(output_root);
    }
    catch (const Exception&)
    {
        EXCEPTION("The location provided to OutputFileHandler must be inside CHASTE_TEST_OUTPUT; '"
                  << rDirectory.GetAbsolutePath() << "' is not under '"
                  << output_root.GetAbsolutePath() << "'.");
    }
    if (*output_root.GetAbsolutePath().rbegin() != '/' && !relative_path.empty())
    {
        assert(*relative_path.begin() == '/');
        relative_path.erase(0, 1); // Remove leading slash
    }
    CommonConstructor(relative_path, cleanOutputDirectory);
}


void OutputFileHandler::CommonConstructor(const std::string& rDirectory,
                                          bool cleanOutputDirectory)
{
    // Is it a valid request for a directory?
    if (rDirectory.find("..") != std::string::npos)
    {
        EXCEPTION("Will not create directory: " + rDirectory +
                  " due to it potentially being above, and cleaning, CHASTE_TEST_OUTPUT.");
        // Note: in Boost 1.48 and above we could use 'canonical' to check this better
    }
    //The notion of absolute path on Windows is rather different.
    //For example, / and /foo are not absolute paths.
    //However, fs::path.has_root_path() captures the intended semantics here as follows

    if (fs::path(rDirectory).has_root_path())
    {
        EXCEPTION("The constructor argument to OutputFileHandler must be a relative path; '"
                  << rDirectory << "' is absolute.");
    }

    mDirectory = MakeFoldersAndReturnFullPath(rDirectory);

    // Clean the directory (default)
    if (rDirectory != "" && cleanOutputDirectory) // Clean directory but don't ever clean CHASTE_TEST_OUTPUT at the top level
    {
        FileFinder signature_file(mDirectory + SIG_FILE_NAME, RelativeTo::Absolute);
        if (!signature_file.Exists())
        {
            EXCEPTION("Cannot delete " + mDirectory + " because signature file \"" + SIG_FILE_NAME + "\" is not present.");
        }

        // Are we the master process? Only the master should delete files
        if (PetscTools::AmMaster())
        {
            CleanFolder(mDirectory);
        }
        // Wait for master to finish before going on to use the directory.
        PetscTools::Barrier("OutputFileHandler");
    }
}

std::string OutputFileHandler::GetChasteTestOutputDirectory()
{
    char* chaste_test_output = getenv("CHASTE_TEST_OUTPUT");
    FileFinder directory_root;
    if (chaste_test_output == nullptr || *chaste_test_output == 0)
    {
        // Mimic the old SCons behaviour of setting CHASTE_TEST_OUTPUT: /tmp/'+os.environ['USER']+'/testoutput/
        std::stringstream  tmp_directory;
        if (getenv("USER")!=NULL)
        {
            tmp_directory << "/tmp/" << getenv("USER") << "/testoutput/";
        }
        else
        {
            // No $USER in environment (which may be the case in Docker)
            tmp_directory << "/tmp/chaste/testoutput/"; // LCOV_EXCL_LINE
        }
        directory_root.SetPath(tmp_directory.str(), RelativeTo::AbsoluteOrCwd);
        /* // Former behaviour: default to 'testoutput' folder within the current directory
         * directory_root.SetPath("testoutput", RelativeTo::CWD); */
    }
    else
    {
        directory_root.SetPath(chaste_test_output, RelativeTo::AbsoluteOrCwd);
    }
    // Note that FileFinder::GetAbsolutePath adds a trailing slash, but only
    // if the directory exists at the time of the call
    std::string chaste_test_output_directory = directory_root.GetAbsolutePath();
    AddTrailingSlash(chaste_test_output_directory);
    return chaste_test_output_directory;
}

std::string OutputFileHandler::MakeFoldersAndReturnFullPath(const std::string& rDirectory) const
{
    fs::path output_root(GetChasteTestOutputDirectory());
    fs::path rel_path(rDirectory);

    if (!rel_path.empty() && (*(--rel_path.end())) == ".")
    {
        // rDirectory has a trailing slash, which gives an unhelpful last component
        rel_path.remove_leaf();
    }

    // Make master wait (because other processes may be checking whether a directory exists)
    PetscTools::Barrier("OutputFileHandler::MakeFoldersAndReturnFullPathBeforeCreation");
    // Are we the master process? Only the master should make any new directories
    if (PetscTools::AmMaster())
    {
        try
        {
            // If necessary make the ChasteTestOutputDirectory - don't make it deleteable by Chaste
            fs::create_directories(output_root); // Note that this is a no-op if the folder exists already

            // Now make all the sub-folders requested one-by-one and add the .chaste_deletable_folder file to them
            fs::path next_folder(output_root);
            for (fs::path::iterator path_iter = rel_path.begin(); path_iter != rel_path.end(); ++path_iter)
            {
                next_folder /= *path_iter;
                bool created_dir = fs::create_directory(next_folder);
                if (created_dir)
                {
                    // Add the Chaste signature file
                    fs::ofstream sig_file(next_folder / SIG_FILE_NAME);
                    sig_file.close();
                }
            }
        }
        // LCOV_EXCL_START
        catch (const fs::filesystem_error& e)
        {
            TERMINATE("Error making test output folder: " << e.what());
        }
        // LCOV_EXCL_STOP
    }

    // Wait for master to finish before going on to use the directory.
    PetscTools::Barrier("OutputFileHandler::MakeFoldersAndReturnFullPath");

    std::string path_with_slash = (output_root / rel_path).string();
    AddTrailingSlash(path_with_slash);
    return path_with_slash;
}

std::string OutputFileHandler::GetOutputDirectoryFullPath() const
{
    return mDirectory;
}

std::string OutputFileHandler::GetRelativePath() const
{
    FileFinder output_root("", RelativeTo::ChasteTestOutput);
    std::string relative_path = FindFile("").GetRelativePath(output_root);
    if (!relative_path.empty() && *relative_path.rbegin() == '/')
    {
        relative_path.erase(--relative_path.end()); // Remove trailing slash
    }
    return relative_path;
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
    fs::path from_path(rSourceFile.GetAbsolutePath());
    fs::path to_path(GetOutputDirectoryFullPath());
    to_path /= from_path.leaf();
    if (PetscTools::AmMaster())
    {
        try
        {
            fs::copy_file(from_path, to_path);
        }
        // LCOV_EXCL_START
        catch (const fs::filesystem_error& e)
        {
            TERMINATE("Error copying file '" << rSourceFile.GetAbsolutePath() << "': " << e.what());
        }
        // LCOV_EXCL_STOP
    }
    PetscTools::Barrier("OutputFileHandler::CopyFileTo");
    return FileFinder(to_path.string(), RelativeTo::Absolute);
}

FileFinder OutputFileHandler::FindFile(std::string leafName) const
{
    return FileFinder(GetOutputDirectoryFullPath() + leafName, RelativeTo::Absolute);
}
