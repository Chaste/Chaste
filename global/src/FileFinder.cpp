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

#include "FileFinder.hpp"

#include <cassert>

#include "BoostFilesystem.hpp"
#include "ChasteBuildRoot.hpp"
#include "Exception.hpp"
#include "GetCurrentWorkingDirectory.hpp"
#include "OutputFileHandler.hpp"
#include "PosixPathFixer.hpp"
#include "Warnings.hpp"

bool FileFinder::msFaking = false;

RelativeTo::Value FileFinder::msFakeWhat = RelativeTo::Absolute;

std::string FileFinder::msFakePath = "";

#define UNSET_PATH "UNSET!"

/**
 * This macro converts boost filesystem exceptions into Chaste Exceptions.
 *
 * It should be put round any calls in this class' public methods to
 * either:
 *  * boost filesystem
 *  * or private methods of this class
 *
 * that are likely to throw.
 *
 * @param code some code that could throw a boost file system error
 */
#define CONVERT_ERROR(code)               \
    try                                   \
    {                                     \
        code;                             \
    }                                     \
    catch (const fs::filesystem_error& e) \
    {                                     \
        EXCEPTION(e.what());              \
    }

FileFinder::FileFinder()
        : mAbsPath(UNSET_PATH)
{
}

FileFinder::FileFinder(const std::string& rRelativePath, RelativeTo::Value relativeTo)
{
    SetPath(rRelativePath, relativeTo);
}

FileFinder::FileFinder(const std::string& rLeafName, const FileFinder& rParentOrSibling)
{
    SetPath(rLeafName, rParentOrSibling);
}

FileFinder::FileFinder(const fs::path& rPath)
{
    SetPath(fs::complete(rPath).string(), RelativeTo::Absolute);
}

FileFinder::~FileFinder()
{
}

void FileFinder::SetPath(const std::string& rRelativePath, RelativeTo::Value relativeTo)
{
    switch (relativeTo)
    {
        case RelativeTo::ChasteSourceRoot:
            mAbsPath = ChasteSourceRootDir() + rRelativePath;
            break;

        case RelativeTo::ChasteBuildRoot:
            mAbsPath = ChasteBuildRootDir() + rRelativePath;
            break;

        case RelativeTo::ChasteTestOutput:
            mAbsPath = OutputFileHandler::GetChasteTestOutputDirectory() + rRelativePath;
            break;

        case RelativeTo::CWD:
            mAbsPath = GetCurrentWorkingDirectory() + "/" + rRelativePath;
            break;

        case RelativeTo::Absolute:
            mAbsPath = rRelativePath;
            break;

        case RelativeTo::AbsoluteOrCwd:
            if (FileFinder::IsAbsolutePath(rRelativePath))
            {
                mAbsPath = rRelativePath;
            }
            else
            {
                mAbsPath = GetCurrentWorkingDirectory() + "/" + rRelativePath;
            }
            break;

        default:
            // Getting here is impossible
            NEVER_REACHED;
            break;
    }

    if (msFaking && msFakeWhat == relativeTo)
    {
        // Fake the resulting path
        mAbsPath = msFakePath + "/" + rRelativePath;
    }

    // Remove any trailing /
    std::string::iterator it = mAbsPath.end();
    while (it != mAbsPath.begin() && *(--it) == '/')
    {
        // Iterator was decremented in the while test
    }
    // it now points at the last non-slash character, if any
    if (it != mAbsPath.end() && (++it) != mAbsPath.end())
    {
        mAbsPath.erase(it, mAbsPath.end());
    }
}

void FileFinder::SetPath(const std::string& rLeafName, const FileFinder& rParentOrSibling)
{
    if (!rParentOrSibling.Exists())
    {
        EXCEPTION("Reference path '" << rParentOrSibling.GetAbsolutePath() << "' does not exist.");
    }
    if (rParentOrSibling.IsDir())
    {
        SetPath(rParentOrSibling.GetAbsolutePath() + rLeafName, RelativeTo::Absolute);
    }
    else
    {
        SetPath(rParentOrSibling.GetParent().GetAbsolutePath() + rLeafName, RelativeTo::Absolute);
    }
}

bool FileFinder::IsPathSet() const
{
    return mAbsPath != UNSET_PATH;
}

bool FileFinder::Exists() const
{
    return fs::exists(mAbsPath);
}

bool FileFinder::IsFile() const
{
    return fs::is_regular(mAbsPath);
}

bool FileFinder::IsDir() const
{
    return fs::is_directory(mAbsPath);
}

bool FileFinder::IsEmpty() const
{
    bool empty = true;
    if (IsFile())
    {
        empty = (fs::file_size(mAbsPath) == 0u);
    }
    else if (IsDir())
    {
        fs::directory_iterator end_iter;
        for (fs::directory_iterator dir_iter(mAbsPath); dir_iter != end_iter; ++dir_iter)
        {
            if (PATH_LEAF_NAME(dir_iter->path()).substr(0, 1) != ".")
            {
                empty = false;
                break;
            }
        }
    }
    else
    {
        EXCEPTION("The path '" << mAbsPath << "' does not exist.");
    }
    return empty;
}

std::string FileFinder::GetAbsolutePath() const
{
    if (IsDir())
    {
        return mAbsPath + '/';
    }
    return mAbsPath;
}

bool FileFinder::IsNewerThan(const FileFinder& rOtherEntity) const
{
    assert(Exists());
    assert(rOtherEntity.Exists());
    return fs::last_write_time(mAbsPath) > fs::last_write_time(rOtherEntity.mAbsPath);
}

std::string FileFinder::GetLeafName() const
{
    return PATH_LEAF_NAME(fs::path(mAbsPath));
}

std::string FileFinder::GetLeafNameNoExtension() const
{
    return fs::basename(mAbsPath);
}

std::string FileFinder::GetExtension() const
{
    return fs::extension(mAbsPath);
}

FileFinder FileFinder::GetParent() const
{
    fs::path our_path(mAbsPath);
    EXCEPT_IF_NOT(our_path.has_branch_path());
    return FileFinder(our_path.branch_path().string(),
                      RelativeTo::Absolute);
}

std::string FileFinder::GetRelativePath(const FileFinder& rBasePath) const
{
    const std::string base_path = rBasePath.GetAbsolutePath();
    const std::string our_path = GetAbsolutePath();
    if (our_path.substr(0, base_path.length()) != base_path)
    {
        EXCEPTION("The path '" << our_path << "' is not relative to '" << base_path << "'.");
    }
    return our_path.substr(base_path.length());
}

/**
 * Helper function for FileFinder::CopyTo - recursively copy the given path.
 * @param rFromPath
 * @param rToPath
 */
void RecursiveCopy(const fs::path& rFromPath, const fs::path& rToPath);

void RecursiveCopy(const fs::path& rFromPath, const fs::path& rToPath)
{
    fs::path dest = rToPath;
    // If rToPath is a folder, then we're copying to the source name *inside* this folder
    if (fs::is_directory(dest))
    {
        dest /= rFromPath.leaf();
    }
    // If the source is a folder, it's complicated
    if (fs::is_directory(rFromPath))
    {
        // Create the destination folder
        EXCEPT_IF_NOT(!fs::exists(dest));
        fs::create_directory(dest);
        // Recursively copy our contents
        fs::directory_iterator end_iter;
        for (fs::directory_iterator dir_iter(rFromPath); dir_iter != end_iter; ++dir_iter)
        {
            RecursiveCopy(dir_iter->path(), dest);
        }
    }
    else
    {
        fs::copy_file(rFromPath, dest); // Just copy!
    }
}

FileFinder FileFinder::CopyTo(const FileFinder& rDest) const
{
    if (!Exists())
    {
        EXCEPTION("Cannot copy '" << mAbsPath << "' as it does not exist.");
    }
    fs::path from_path(mAbsPath);
    fs::path to_path(rDest.mAbsPath);
    if (rDest.IsDir())
    {
        to_path /= from_path.leaf();
    }
    if (fs::exists(to_path))
    {
        if (IsFile())
        {
            CONVERT_ERROR(fs::remove(to_path));
        }
        else
        {
            EXCEPTION("Cannot copy '" << mAbsPath << "' to '" << to_path << "' as it would overwrite an existing file.");
        }
    }
    CONVERT_ERROR(RecursiveCopy(from_path, to_path));
    return FileFinder(to_path);
}

/**
 * Helper function for FileFinder::Remove - recursively remove the given path.
 * @param rPath
 */
void RemoveAll(const fs::path& rPath);

void RemoveAll(const fs::path& rPath)
{
    // First recursively remove any children
    if (fs::is_directory(rPath))
    {
        fs::directory_iterator end_iter;
        for (fs::directory_iterator dir_iter(rPath); dir_iter != end_iter; ++dir_iter)
        {
            RemoveAll(dir_iter->path());
        }
    }
    // Now remove the item itself
    fs::remove(rPath);
}

void FileFinder::PrivateRemove(bool dangerous) const
{
    // Test for bad paths
    const std::string test_output(OutputFileHandler::GetChasteTestOutputDirectory());
    const std::string test_output_path(ChastePosixPathFixer::ToPosix(fs::path(test_output)));
    const std::string absolute_path(ChastePosixPathFixer::ToPosix(fs::path(GetAbsolutePath())));
    bool in_testoutput = (absolute_path.substr(0, test_output_path.length()) == test_output_path);

    if (!in_testoutput)
    {
        if (dangerous)
        {
            const std::string source_folder(FileFinder("", RelativeTo::ChasteSourceRoot).GetAbsolutePath());
            const std::string source_folder_path = ChastePosixPathFixer::ToPosix(fs::path(source_folder));
            bool in_source = (absolute_path.substr(0, source_folder_path.length()) == source_folder_path);

            const std::string build_folder(FileFinder("", RelativeTo::ChasteBuildRoot).GetAbsolutePath());
            const std::string build_folder_path = ChastePosixPathFixer::ToPosix(fs::path(build_folder));
            bool in_build = (absolute_path.substr(0, build_folder_path.length()) == build_folder_path);

            if (!(in_source || in_build))
            {
                EXCEPTION("Cannot remove location '" << mAbsPath
                                                     << "' as it is not located within the Chaste test output folder ("
                                                     << test_output_path << "), the Chaste source folder ("
                                                     << source_folder_path << ") or the Chaste build folder ("
                                                     << build_folder_path << ").");
            }
        }
        else
        {
            EXCEPTION("Cannot remove location '" << mAbsPath
                                                 << "' as it is not located within the Chaste test output folder ("
                                                 << test_output_path << ").");
        }
    }

    if (mAbsPath.find("..") != std::string::npos)
    {
        EXCEPTION("Cannot remove location '" << mAbsPath
                                             << "' as it contains a dangerous path component.");
    }
    if (Exists())
    {
        if (!dangerous)
        {
            fs::path sig_file(mAbsPath);
            if (IsFile())
            {
                // We need to look for the signature file in the parent folder
                sig_file.remove_leaf();
            }
            sig_file /= OutputFileHandler::SIG_FILE_NAME;
            if (!fs::exists(sig_file))
            {
                EXCEPTION("Cannot remove location '" << mAbsPath << "' because the signature file '"
                                                     << OutputFileHandler::SIG_FILE_NAME << "' is not present.");
            }
        }
        // Do the removal
        CONVERT_ERROR(RemoveAll(mAbsPath));
    }
}

void FileFinder::Remove() const
{
    PrivateRemove();
}

void FileFinder::DangerousRemove() const
{
    PrivateRemove(true);
}

std::vector<FileFinder> FileFinder::FindMatches(const std::string& rPattern) const
{
    // Check for error/warning cases
    if (!IsDir())
    {
        EXCEPTION("Cannot search for matching files in '" << mAbsPath << "' as it is not a directory.");
    }
    size_t len = rPattern.length();
    size_t inner_star_pos = rPattern.find('*', 1);
    if (inner_star_pos != std::string::npos && inner_star_pos < len - 1)
    {
        WARNING("A '*' only has special meaning at the start or end of a pattern.");
    }

    // Note initial or trailing *, and use of ?
    std::string pattern(rPattern);
    bool star_fini = false;
    if (!pattern.empty() && *(pattern.rbegin()) == '*')
    {
        star_fini = true;
        pattern = pattern.substr(0, len - 1);
        len--;
    }
    bool star_init = false;
    if (!pattern.empty() && pattern[0] == '*')
    {
        star_init = true;
        pattern = pattern.substr(1);
        len--;
    }
    bool has_query = (pattern.find('?') != std::string::npos);
    // Disallow a harder case to match
    if (star_init && star_fini && has_query)
    {
        EXCEPTION("The '*' wildcard may not be used at both the start and end of the pattern if the '?' wildcard is also used.");
    }

    // Search the folder
    std::vector<FileFinder> results;
    if (!rPattern.empty())
    {
        fs::directory_iterator end_iter;
        fs::path our_path(mAbsPath);
        for (fs::directory_iterator dir_iter(our_path); dir_iter != end_iter; ++dir_iter)
        {
            std::string leafname = PATH_LEAF_NAME(dir_iter->path());
            size_t leaf_len = leafname.length();
            if (leafname[0] != '.' // Don't include hidden files
                && leaf_len >= len) // Ignore stuff that can't match
            {
                if (!has_query) // Easier case
                {
                    size_t pos = leafname.find(pattern);
                    if ((star_init || pos == 0) && (star_fini || pos + len == leaf_len))
                    {
                        results.push_back(FileFinder(our_path / leafname));
                    }
                }
                else
                {
                    std::string match;
                    if (star_init)
                    {
                        // Match against last len chars
                        match = leafname.substr(leaf_len - len);
                    }
                    else
                    {
                        // Match against first len chars
                        match = leafname.substr(0, len);
                    }
                    bool ok = true;
                    for (std::string::const_iterator it_p = pattern.begin(), it_m = match.begin();
                         it_p != pattern.end();
                         ++it_p, ++it_m)
                    {
                        if (*it_p != '?' && *it_p != *it_m)
                        {
                            ok = false;
                            break;
                        }
                    }
                    if (ok)
                    {
                        results.push_back(FileFinder(our_path / leafname));
                    }
                }
            }
        }
    }

    std::sort(results.begin(), results.end());
    return results;
}

bool FileFinder::IsAbsolutePath(const std::string& rPath)
{
    return fs::path(rPath).is_complete();
}

void FileFinder::ReplaceSpacesWithUnderscores(std::string& rPath)
{
    for (std::string::iterator it = rPath.begin(); it != rPath.end(); ++it)
    {
        if (*it == ' ')
        {
            *it = '_';
        }
    }
}

void FileFinder::ReplaceUnderscoresWithSpaces(std::string& rPath)
{
    for (std::string::iterator it = rPath.begin(); it != rPath.end(); ++it)
    {
        if (*it == '_')
        {
            *it = ' ';
        }
    }
}

bool FileFinder::operator<(const FileFinder& otherFinder) const
{
    return (mAbsPath < otherFinder.GetAbsolutePath());
}

void FileFinder::FakePath(RelativeTo::Value fakeWhat, const std::string& rFakePath)
{
    msFakeWhat = fakeWhat;
    msFakePath = rFakePath;
    msFaking = true;
}

void FileFinder::StopFaking()
{
    msFaking = false;
}
