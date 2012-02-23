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

#include "FileFinder.hpp"
#include "ChasteBuildRoot.hpp"
#include "OutputFileHandler.hpp"
#include "Exception.hpp"
#include "GetCurrentWorkingDirectory.hpp"
#include <fstream>
#include <cassert>
#include <sys/stat.h>

bool FileFinder::msFaking = false;

RelativeTo::Value FileFinder::msFakeWhat = RelativeTo::Absolute;

std::string FileFinder::msFakePath = "";


FileFinder::FileFinder()
    : mAbsPath("UNSET!")
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

void FileFinder::SetPath(const std::string& rRelativePath, RelativeTo::Value relativeTo)
{
    switch (relativeTo)
    {
        case RelativeTo::ChasteSourceRoot:
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

    if (IsDir())
    {
        if (*(mAbsPath.end()-1) != '/')
        {
            mAbsPath = mAbsPath + "/";
        }
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

bool FileFinder::Exists() const
{
    struct stat our_stats;
    int retcode = stat(GetAbsolutePath().c_str(), &our_stats);
    return (retcode == 0);
}

bool FileFinder::IsFile() const
{
    bool result;
    struct stat our_stats;
    int retcode = stat(GetAbsolutePath().c_str(), &our_stats);
    if (retcode == 0)
    {
        result = S_ISREG(our_stats.st_mode);
    }
    else
    {
        // If it doesn't exist, it isn't a file
        result = false;
    }
    return result;
}

bool FileFinder::IsDir() const
{
    bool result;
    struct stat our_stats;
    int retcode = stat(GetAbsolutePath().c_str(), &our_stats);
    if (retcode == 0)
    {
        result = S_ISDIR(our_stats.st_mode);
    }
    else
    {
        // If it doesn't exist, it isn't a directory
        result = false;
    }
    return result;
}

std::string FileFinder::GetAbsolutePath() const
{
    return mAbsPath;
}

bool FileFinder::IsNewerThan(const FileFinder& rOtherEntity) const
{
    assert(Exists());
    assert(rOtherEntity.Exists());
    struct stat our_stats, other_stats;
    stat(GetAbsolutePath().c_str(), &our_stats);
    stat(rOtherEntity.GetAbsolutePath().c_str(), &other_stats);
    return our_stats.st_mtime > other_stats.st_mtime;
}

std::string FileFinder::GetLeafName() const
{
    std::string full_name = GetAbsolutePath();
    size_t slash = full_name.rfind('/');
    EXCEPT_IF_NOT(slash != std::string::npos);
    return full_name.substr(slash+1);
}

std::string FileFinder::GetLeafNameNoExtension() const
{
    std::string leaf_name = GetLeafName();
    size_t dot = leaf_name.rfind('.');
    return leaf_name.substr(0,dot);
}

FileFinder FileFinder::GetParent() const
{
    std::string full_name = GetAbsolutePath();
    size_t limit = full_name.length() > 1 ? full_name.length() - 2 : std::string::npos;
    size_t slash = full_name.rfind('/', limit);
    EXCEPT_IF_NOT(slash != std::string::npos);
    return FileFinder(full_name.substr(0, slash), RelativeTo::Absolute);
}

bool FileFinder::IsAbsolutePath(const std::string& rPath)
{
    return rPath[0]=='/';
}

void FileFinder::ReplaceSpacesWithUnderscores(std::string& rPath)
{
    for(std::string::iterator it = rPath.begin(); it != rPath.end(); ++it)
    {
        if(*it == ' ')
        {
            *it = '_';
        }
    }
}

void FileFinder::ReplaceUnderscoresWithSpaces(std::string& rPath)
{
    for(std::string::iterator it = rPath.begin(); it != rPath.end(); ++it)
    {
        if(*it == '_')
        {
            *it = ' ';
        }
    }
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
