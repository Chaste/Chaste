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

#include "ArchiveLocationInfo.hpp"

#include <sstream>

#include "Exception.hpp"
#include "OutputFileHandler.hpp"

std::string ArchiveLocationInfo::mDirAbsPath = "";
std::string ArchiveLocationInfo::mMeshFilename = "mesh";

void ArchiveLocationInfo::SetMeshPathname(const FileFinder& rDirectory, const std::string& rFilename)
{
    SetArchiveDirectory(rDirectory);
    mMeshFilename = rFilename;
}

void ArchiveLocationInfo::SetMeshPathname(const std::string& rDirectory, const std::string& rFilename)
{
    bool is_absolute = FileFinder::IsAbsolutePath(rDirectory);
    RelativeTo::Value relative_to;
    if (!is_absolute)
    {
        relative_to = RelativeTo::ChasteTestOutput;
    }
    else
    {
        relative_to = RelativeTo::Absolute;
    }
    FileFinder dir(rDirectory, relative_to);
    SetArchiveDirectory(dir);
    mMeshFilename = rFilename;
}

void ArchiveLocationInfo::SetMeshFilename(const std::string& rFilename)
{
    mMeshFilename = rFilename;
}

std::string ArchiveLocationInfo::GetMeshFilename()
{
    if (mMeshFilename == "")
    {
        EXCEPTION("ArchiveLocationInfo::mMeshFilename has not been set");
    }
    return mMeshFilename;
}

std::string ArchiveLocationInfo::GetArchiveDirectory()
{
    if (mDirAbsPath == "")
    {
        EXCEPTION("ArchiveLocationInfo::mDirAbsPath has not been set");
    }
    return mDirAbsPath;
}

std::string ArchiveLocationInfo::GetProcessUniqueFilePath(
        const std::string& rFileName,
        unsigned procId)
{
    std::stringstream filepath;
    filepath << GetArchiveDirectory() << rFileName << "." << procId;
    return filepath.str();
}

void ArchiveLocationInfo::SetArchiveDirectory(const FileFinder& rDirectory)
{
    mDirAbsPath = rDirectory.GetAbsolutePath();
    if (!(*(mDirAbsPath.end()-1) == '/'))
    {
        mDirAbsPath = mDirAbsPath + "/";
    }
}

std::string ArchiveLocationInfo::GetArchiveRelativePath()
{
    std::string chaste_output = OutputFileHandler::GetChasteTestOutputDirectory();

    // Find occurrence of CHASTE_TEST_OUTPUT in string
    std::string::size_type pos = mDirAbsPath.find(chaste_output, 0);
    if (pos == 0)
    {
        return mDirAbsPath.substr(chaste_output.length());
    }
    else
    {
        return mDirAbsPath;
    }
}

bool ArchiveLocationInfo::GetIsDirRelativeToChasteTestOutput()
{
    std::string chaste_output = OutputFileHandler::GetChasteTestOutputDirectory();
    std::string::size_type pos = mDirAbsPath.find(chaste_output, 0);
    return (pos == 0);
}
